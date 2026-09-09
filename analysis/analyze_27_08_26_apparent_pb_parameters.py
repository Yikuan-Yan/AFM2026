#!/usr/bin/env python3
"""Brute-force apparent PB fits for the 27-08-26 pure-water palindrome maps.

This intentionally does *not* subtract a hydrodynamic term or correct the
observed acquisition-history relaxation.  It applies the previously used
equal-potential silica sphere-plane nonlinear-PB + van der Waals model to each
finite-speed curve and asks whether the resulting *apparent* Debye length and
surface-potential magnitude vary with clock time or approach speed.

The expensive nonlinear PB pressure tables are shared by all curves through a
common lambda-zeta grid.  This is a deliberately coarse/grid-based analysis,
not a claim that the fitted parameters are equilibrium material constants.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats
from scipy.constants import Boltzmann, elementary_charge, epsilon_0


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import fit_glycerol_surface_forces as base  # noqa: E402


UPSTREAM = ROOT / "analysis" / "palindrome_27_08_26_full_results"
RESULTS = ROOT / "analysis" / "palindrome_27_08_26_apparent_pb_results"
FIGURES = RESULTS / "figures"

MODEL = "nonlinear_pb_derjaguin"
EPSILON_R_WATER = 78.5
BASELINES = ("linear_drift_corrected", "far_constant_referenced")
FIT_VARIANTS = {
    "previous_primary_10_250nm": (10.0, 250.0),
    "measurement_window_20_200nm": (20.0, 200.0),
}
PRIMARY_VARIANT = "previous_primary_10_250nm"
DISTANCE_GRID_NM = np.arange(10.0, 255.0, 5.0, dtype=np.float64)
LAMBDA_GRID_NM = np.geomspace(1.0, 300.0, 181, dtype=np.float64)
ZETA_GRID_MV = np.linspace(1.0, 250.0, 167, dtype=np.float64)
BASELINE_BOUNDS_PN = (-500.0, 500.0)
MIN_SPATIAL_SIGMA_PN = 2.0
PARAMETERS = (
    ("lambda_D_nm", "Apparent Debye length", "nm"),
    ("zeta_magnitude_mV", "Apparent surface-potential magnitude", "mV"),
)
SPEED_COLORS = {0.05: "#2a9d8f", 0.1: "#e9b949", 0.2: "#e76f51"}
BLOCK_COLORS = {1: "#577590", 2: "#43aa8b", 3: "#f9c74f", 4: "#f9844a", 5: "#9b5de5"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def model_library() -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    """Evaluate the previous nonlinear-PB + vdW formula on one shared grid."""

    distance_nm = DISTANCE_GRID_NM
    distance_m = distance_nm * 1e-9
    lambda_m = LAMBDA_GRID_NM * 1e-9
    kappa = 1.0 / lambda_m
    dimensionless_gap = distance_nm[None, :] / LAMBDA_GRID_NM[:, None]
    absolute_permittivity = EPSILON_R_WATER * epsilon_0
    prefactor = (
        2.0
        * np.pi
        * base.PROBE_RADIUS_M
        * absolute_permittivity
        * (Boltzmann * base.TEMPERATURE_K / elementary_charge) ** 2
    )
    vdw_pN = (
        -base.HAMAKER_J
        * base.PROBE_RADIUS_M
        / (6.0 * distance_m**2)
        * 1e12
    )
    library = np.empty(
        (ZETA_GRID_MV.size, LAMBDA_GRID_NM.size, distance_nm.size),
        dtype=np.float64,
    )
    for zeta_index, zeta_mV in enumerate(ZETA_GRID_MV):
        dimensionless_surface = (
            elementary_charge
            * zeta_mV
            * 1e-3
            / (Boltzmann * base.TEMPERATURE_K)
        )
        dimensionless_integral = base.pb_dimensionless_G(
            dimensionless_gap, dimensionless_surface
        )
        edl_pN = (
            prefactor
            * kappa[:, None]
            * dimensionless_integral
            * 1e12
        )
        library[zeta_index] = edl_pN + vdw_pN[None, :]

    flat_library = library.reshape(-1, distance_nm.size)
    flat_lambda = np.tile(LAMBDA_GRID_NM, ZETA_GRID_MV.size)
    flat_zeta = np.repeat(ZETA_GRID_MV, LAMBDA_GRID_NM.size)

    absolute_errors: list[float] = []
    relative_errors: list[float] = []
    for zeta_index, lambda_index in ((0, 0), (41, 60), (83, 120), (166, 180)):
        reference = base.total_equilibrium_force_pN(
            distance_nm,
            float(LAMBDA_GRID_NM[lambda_index]),
            float(ZETA_GRID_MV[zeta_index]),
            0.0,
            EPSILON_R_WATER,
            MODEL,
        )
        candidate = library[zeta_index, lambda_index]
        absolute_errors.append(float(np.max(np.abs(candidate - reference))))
        relative_errors.append(
            float(
                np.max(
                    np.abs(candidate - reference)
                    / np.maximum(np.abs(reference), 1e-9)
                )
            )
        )
    check = {
        "shared_grid_formula_max_abs_error_pN": max(absolute_errors),
        "shared_grid_formula_max_relative_error": max(relative_errors),
    }
    if check["shared_grid_formula_max_abs_error_pN"] > 1e-8:
        raise AssertionError(f"shared PB grid does not match previous formula: {check}")
    if not np.all(np.isfinite(flat_library)):
        raise FloatingPointError("non-finite PB model library")
    return flat_library, flat_lambda, flat_zeta, check


def normalize_curve_points(rows: list[dict[str, str]], dataset_type: str) -> list[dict]:
    points: list[dict] = []
    for row in rows:
        if dataset_type == "map":
            force = float(row["force_median_pN"])
            q25 = float(row["force_q25_pN"])
            q75 = float(row["force_q75_pN"])
            count = int(row["pixel_count"])
        else:
            force = float(row["median_pN"])
            q25 = float(row["q25_pN"])
            q75 = float(row["q75_pN"])
            count = int(row["paired_pixel_count"])
        values = np.asarray([force, q25, q75], dtype=float)
        if not np.all(np.isfinite(values)):
            continue
        points.append(
            {
                "distance_nm": float(row["distance_nm"]),
                "force_pN": force,
                "q25_pN": q25,
                "q75_pN": q75,
                "count": count,
            }
        )
    return sorted(points, key=lambda item: item["distance_nm"])


def fit_curve_grid(
    points: list[dict],
    model_force: np.ndarray,
    flat_lambda: np.ndarray,
    flat_zeta: np.ndarray,
    variant: str,
) -> tuple[dict, list[dict]]:
    dmin_nm, dmax_nm = FIT_VARIANTS[variant]
    point_by_distance = {float(row["distance_nm"]): row for row in points}
    selected_indices = np.asarray(
        [
            index
            for index, distance in enumerate(DISTANCE_GRID_NM)
            if dmin_nm <= distance <= dmax_nm and distance in point_by_distance
        ],
        dtype=int,
    )
    if selected_indices.size < 20:
        raise RuntimeError(f"{variant}: only {selected_indices.size} usable points")
    selected_distances = DISTANCE_GRID_NM[selected_indices]
    selected = [point_by_distance[float(distance)] for distance in selected_distances]
    force = np.asarray([row["force_pN"] for row in selected], dtype=np.float64)
    spatial_sigma = np.asarray(
        [
            max(
                (row["q75_pN"] - row["q25_pN"]) / 1.349,
                MIN_SPATIAL_SIGMA_PN,
            )
            for row in selected
        ],
        dtype=np.float64,
    )
    if not np.all(np.isfinite(force)) or not np.all(np.isfinite(spatial_sigma)):
        raise FloatingPointError("non-finite force or fit scale")
    if np.any(spatial_sigma <= 0.0):
        raise FloatingPointError("non-positive fit scale")

    candidate = model_force[:, selected_indices]
    weight = 1.0 / spatial_sigma**2
    weight_sum = float(np.sum(weight))
    baseline = ((force[None, :] - candidate) * weight[None, :]).sum(axis=1)
    baseline /= weight_sum
    baseline = np.clip(baseline, *BASELINE_BOUNDS_PN)
    residual = candidate + baseline[:, None] - force[None, :]
    weighted_rss = np.einsum("ij,j,ij->i", residual, weight, residual)
    if not np.all(np.isfinite(weighted_rss)):
        raise FloatingPointError("non-finite weighted RSS")
    best_index = int(np.argmin(weighted_rss))
    best_prediction = candidate[best_index] + baseline[best_index]
    raw_residual = force - best_prediction
    rss = float(np.sum(raw_residual**2))
    total = float(np.sum((force - np.mean(force)) ** 2))
    r_squared = float("nan") if total <= 0.0 else 1.0 - rss / total
    degrees_of_freedom = int(force.size - 3)
    reduced_grid_objective = float(weighted_rss[best_index] / degrees_of_freedom)

    # A diagnostic profile region only.  The spatial IQR scale and adjacent
    # distance bins are not independent measurement errors, so this is not a
    # formal confidence interval.
    profile_delta = 2.30 * max(reduced_grid_objective, 1.0)
    profile = weighted_rss <= weighted_rss[best_index] + profile_delta
    profile_lambda = flat_lambda[profile]
    profile_zeta = flat_zeta[profile]

    lambda_value = float(flat_lambda[best_index])
    zeta_value = float(flat_zeta[best_index])
    baseline_value = float(baseline[best_index])
    lambda_boundary = bool(
        np.isclose(lambda_value, LAMBDA_GRID_NM[0])
        or np.isclose(lambda_value, LAMBDA_GRID_NM[-1])
    )
    zeta_boundary = bool(
        np.isclose(zeta_value, ZETA_GRID_MV[0])
        or np.isclose(zeta_value, ZETA_GRID_MV[-1])
    )
    baseline_boundary = bool(
        np.isclose(baseline_value, BASELINE_BOUNDS_PN[0])
        or np.isclose(baseline_value, BASELINE_BOUNDS_PN[1])
    )
    result = {
        "fit_variant": variant,
        "model": MODEL,
        "dmin_nm": dmin_nm,
        "dmax_nm": dmax_nm,
        "n_distance_bins": int(force.size),
        "lambda_D_nm": lambda_value,
        "lambda_profile_low_nm": float(np.min(profile_lambda)),
        "lambda_profile_high_nm": float(np.max(profile_lambda)),
        "zeta_magnitude_mV": zeta_value,
        "zeta_signed_silica_mV": -zeta_value,
        "zeta_profile_low_mV": float(np.min(profile_zeta)),
        "zeta_profile_high_mV": float(np.max(profile_zeta)),
        "baseline_pN": baseline_value,
        "weighted_grid_objective": float(weighted_rss[best_index]),
        "reduced_grid_objective": reduced_grid_objective,
        "rmse_pN": float(math.sqrt(np.mean(raw_residual**2))),
        "r2": r_squared,
        "lambda_grid_boundary": lambda_boundary,
        "zeta_grid_boundary": zeta_boundary,
        "baseline_boundary": baseline_boundary,
        "fit_status": "grid_boundary"
        if lambda_boundary or zeta_boundary or baseline_boundary
        else "interior_grid_solution",
        "weight_definition": "max(pixel_IQR/1.349,2_pN); spatial_scale_not_SE",
        "profile_interval_status": "diagnostic_correlated_bins_not_formal_CI",
    }
    predictions: list[dict] = []
    edl = base.edl_force_pN(
        selected_distances,
        lambda_value,
        zeta_value,
        EPSILON_R_WATER,
        MODEL,
    )
    distance_m = selected_distances * 1e-9
    vdw = -base.HAMAKER_J * base.PROBE_RADIUS_M / (6.0 * distance_m**2) * 1e12
    for distance, observed, sigma, edl_value, vdw_value, predicted in zip(
        selected_distances,
        force,
        spatial_sigma,
        edl,
        vdw,
        best_prediction,
        strict=True,
    ):
        predictions.append(
            {
                "distance_nm": float(distance),
                "observed_force_pN": float(observed),
                "spatial_sigma_pN": float(sigma),
                "fitted_edl_force_pN": float(edl_value),
                "fixed_vdw_force_pN": float(vdw_value),
                "fitted_baseline_pN": baseline_value,
                "fitted_total_force_pN": float(predicted),
                "residual_pN": float(observed - predicted),
            }
        )
    return result, predictions


def load_curves() -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    inventory = read_csv(UPSTREAM / "map_inventory_QC.csv")
    map_curve_rows = read_csv(UPSTREAM / "map_force_curves.csv")
    pair_curve_rows = read_csv(UPSTREAM / "palindrome_pair_curves.csv")
    inventory_by_source = {row["source"]: row for row in inventory}
    earliest_epoch = min(float(row["map_protocol_midpoint_epoch_s"]) for row in inventory)

    grouped_map: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in map_curve_rows:
        grouped_map[(row["source"], row["baseline_method"])].append(row)
    map_curves: list[dict] = []
    for (source, baseline), rows in grouped_map.items():
        meta = inventory_by_source[source]
        midpoint_epoch = float(meta["map_protocol_midpoint_epoch_s"])
        map_curves.append(
            {
                "dataset_type": "map",
                "curve_id": source,
                "source": source,
                "baseline_method": baseline,
                "acquisition_order": int(meta["acquisition_order"]),
                "block": int(meta["block"]),
                "nominal_speed_um_per_s": float(meta["nominal_speed_um_per_s"]),
                "actual_gap_speed_um_per_s": float(
                    meta["gap_speed_20_200nm_median_um_per_s"]
                ),
                "clock_epoch_s": midpoint_epoch,
                "elapsed_min": (midpoint_epoch - earliest_epoch) / 60.0,
                "points": normalize_curve_points(rows, "map"),
            }
        )

    grouped_pair: dict[tuple[int, float, str], list[dict[str, str]]] = defaultdict(list)
    for row in pair_curve_rows:
        if row["quantity"] != "F_sym":
            continue
        key = (
            int(row["block"]),
            float(row["nominal_speed_um_per_s"]),
            row["baseline_method"],
        )
        grouped_pair[key].append(row)
    pair_curves: list[dict] = []
    for (block, speed, baseline), rows in grouped_pair.items():
        first = rows[0]
        center_epoch = float(first["pair_center_time_epoch_s"])
        pair_curves.append(
            {
                "dataset_type": "palindrome_Fsym",
                "curve_id": f"block{block}_speed{speed:g}_{baseline}",
                "source": "",
                "baseline_method": baseline,
                "acquisition_order": float(first["pair_midpoint_acquisition_order"]),
                "block": block,
                "nominal_speed_um_per_s": speed,
                "actual_gap_speed_um_per_s": float(
                    first["actual_gap_speed_um_per_s"]
                ),
                "clock_epoch_s": center_epoch,
                "elapsed_min": (center_epoch - earliest_epoch) / 60.0,
                "pair_half_span_min": float(first["pair_half_span_min"]),
                "early_source": first["early_source"],
                "late_source": first["late_source"],
                "points": normalize_curve_points(rows, "pair"),
            }
        )
    return inventory, map_curve_rows, map_curves, pair_curves


def fit_all_curves(
    curves: list[dict],
    model_force: np.ndarray,
    flat_lambda: np.ndarray,
    flat_zeta: np.ndarray,
) -> tuple[list[dict], list[dict]]:
    fit_rows: list[dict] = []
    prediction_rows: list[dict] = []
    for curve in sorted(
        curves,
        key=lambda row: (
            row["dataset_type"],
            row["baseline_method"],
            row["acquisition_order"],
        ),
    ):
        for variant in FIT_VARIANTS:
            fit, predictions = fit_curve_grid(
                curve["points"], model_force, flat_lambda, flat_zeta, variant
            )
            metadata = {key: value for key, value in curve.items() if key != "points"}
            fit_rows.append({**metadata, **fit})
            if variant == PRIMARY_VARIANT:
                for prediction in predictions:
                    prediction_rows.append(
                        {
                            **metadata,
                            "fit_variant": variant,
                            **prediction,
                        }
                    )
    return fit_rows, prediction_rows


def ols_hc3(design: np.ndarray, response: np.ndarray) -> dict:
    design = np.asarray(design, dtype=np.float64)
    response = np.asarray(response, dtype=np.float64)
    if design.ndim != 2 or response.ndim != 1 or design.shape[0] != response.size:
        raise ValueError("OLS shape mismatch")
    coefficients, _, rank, _ = np.linalg.lstsq(design, response, rcond=None)
    fitted = design @ coefficients
    residual = response - fitted
    n, parameter_count = design.shape
    residual_df = n - int(rank)
    xtx_inverse = np.linalg.pinv(design.T @ design)
    leverage = np.einsum("ij,jk,ik->i", design, xtx_inverse, design)
    hc3_scale = residual**2 / np.maximum(1.0 - leverage, 1e-10) ** 2
    hc3_covariance = xtx_inverse @ (design.T @ (hc3_scale[:, None] * design)) @ xtx_inverse
    hc3_se = np.sqrt(np.maximum(0.0, np.diag(hc3_covariance)))
    rss = float(residual @ residual)
    total = float(np.sum((response - np.mean(response)) ** 2))
    condition_columns = np.linalg.norm(design, axis=0)
    scaled_design = design / np.where(condition_columns > 0.0, condition_columns, 1.0)
    return {
        "coefficients": coefficients,
        "hc3_se": hc3_se,
        "fitted": fitted,
        "residual": residual,
        "n": n,
        "parameter_count": parameter_count,
        "rank": int(rank),
        "residual_df": residual_df,
        "r2": float("nan") if total <= 0.0 else 1.0 - rss / total,
        "rmse": float(math.sqrt(rss / max(residual_df, 1))),
        "condition": float(np.linalg.cond(scaled_design)),
    }


def time_dependence_rows(map_fits: list[dict]) -> list[dict]:
    rows: list[dict] = []
    for baseline in BASELINES:
        for variant in FIT_VARIANTS:
            selected_all = [
                row
                for row in map_fits
                if row["baseline_method"] == baseline
                and row["fit_variant"] == variant
            ]
            for subset_name, selected in (
                ("all_29_maps", selected_all),
                (
                    "no_refresh_testB_blocks3_5",
                    [row for row in selected_all if int(row["block"]) in (3, 4, 5)],
                ),
            ):
                selected = sorted(selected, key=lambda row: row["elapsed_min"])
                elapsed_hours = np.asarray(
                    [float(row["elapsed_min"]) / 60.0 for row in selected]
                )
                centered_hours = elapsed_hours - np.mean(elapsed_hours)
                design = np.column_stack([np.ones(len(selected)), centered_hours])
                for parameter, _, _ in PARAMETERS:
                    response = np.log(
                        np.asarray([float(row[parameter]) for row in selected])
                    )
                    spearman = stats.spearmanr(elapsed_hours, response)
                    fit = ols_hc3(design, response)
                    beta = float(fit["coefficients"][1])
                    se = float(fit["hc3_se"][1])
                    tcrit = float(stats.t.ppf(0.975, fit["residual_df"]))
                    low = beta - tcrit * se
                    high = beta + tcrit * se
                    pvalue = (
                        float(2.0 * stats.t.sf(abs(beta / se), fit["residual_df"]))
                        if se > 0.0
                        else float("nan")
                    )
                    rows.append(
                        {
                            "analysis": "clock_time_dependence",
                            "dataset_subset": subset_name,
                            "baseline_method": baseline,
                            "fit_variant": variant,
                            "parameter": parameter,
                            "n": len(selected),
                            "spearman_rho": float(spearman.statistic),
                            "spearman_two_sided_p": float(spearman.pvalue),
                            "linear_beta_log_per_hour": beta,
                            "linear_HC3_SE_log_per_hour": se,
                            "linear_HC3_95CI_low_log_per_hour": low,
                            "linear_HC3_95CI_high_log_per_hour": high,
                            "linear_two_sided_p": pvalue,
                            "linear_percent_change_per_hour": 100.0 * math.expm1(beta),
                            "linear_percent_change_CI_low": 100.0 * math.expm1(low),
                            "linear_percent_change_CI_high": 100.0 * math.expm1(high),
                            "linear_R2": fit["r2"],
                            "design_condition": fit["condition"],
                            "claim_status": "descriptive_time_association",
                        }
                    )
    return rows


def time_adjusted_speed_rows(map_fits: list[dict]) -> list[dict]:
    rows: list[dict] = []
    for baseline in BASELINES:
        for variant in FIT_VARIANTS:
            selected = sorted(
                [
                    row
                    for row in map_fits
                    if row["baseline_method"] == baseline
                    and row["fit_variant"] == variant
                    and int(row["block"]) in (3, 4, 5)
                ],
                key=lambda row: row["elapsed_min"],
            )
            elapsed = np.asarray([float(row["elapsed_min"]) for row in selected])
            time_normalized = elapsed - np.mean(elapsed)
            time_normalized /= np.max(np.abs(time_normalized))
            speed = np.asarray(
                [float(row["actual_gap_speed_um_per_s"]) for row in selected]
            )
            speed_centered = speed - np.mean(speed)
            blocks = np.asarray([int(row["block"]) for row in selected])
            block_matrix = np.column_stack(
                [blocks == block for block in (3, 4, 5)]
            ).astype(float)
            within_time = np.zeros_like(elapsed)
            for block in (3, 4, 5):
                mask = blocks == block
                within_time[mask] = elapsed[mask] - np.mean(elapsed[mask])
            within_time /= np.max(np.abs(within_time))
            within_time_squared = within_time**2
            for block in (3, 4, 5):
                mask = blocks == block
                within_time_squared[mask] -= np.mean(within_time_squared[mask])

            designs: list[tuple[str, np.ndarray, int]] = []
            for degree in range(1, 4):
                time_terms = np.polynomial.legendre.legvander(
                    time_normalized, degree
                )[:, 1:]
                designs.append(
                    (
                        f"global_time_poly{degree}",
                        np.column_stack(
                            [np.ones(len(selected)), speed_centered, time_terms]
                        ),
                        1,
                    )
                )
            designs.extend(
                [
                    (
                        "block_fixed_linear_within",
                        np.column_stack([block_matrix, speed_centered, within_time]),
                        3,
                    ),
                    (
                        "block_fixed_quadratic_within",
                        np.column_stack(
                            [
                                block_matrix,
                                speed_centered,
                                within_time,
                                within_time_squared,
                            ]
                        ),
                        3,
                    ),
                ]
            )
            for parameter, _, _ in PARAMETERS:
                response = np.log(
                    np.asarray([float(row[parameter]) for row in selected])
                )
                for model_name, design, speed_column in designs:
                    fit = ols_hc3(design, response)
                    beta = float(fit["coefficients"][speed_column])
                    se = float(fit["hc3_se"][speed_column])
                    tcrit = float(stats.t.ppf(0.975, fit["residual_df"]))
                    low = beta - tcrit * se
                    high = beta + tcrit * se
                    pvalue = (
                        float(2.0 * stats.t.sf(abs(beta / se), fit["residual_df"]))
                        if se > 0.0
                        else float("nan")
                    )
                    delta_speed = 0.15
                    rows.append(
                        {
                            "analysis": "time_adjusted_speed_dependence",
                            "dataset_subset": "no_refresh_testB_blocks3_5",
                            "baseline_method": baseline,
                            "fit_variant": variant,
                            "parameter": parameter,
                            "time_model": model_name,
                            "n": len(selected),
                            "speed_beta_log_per_um_per_s": beta,
                            "speed_HC3_SE_log_per_um_per_s": se,
                            "speed_HC3_95CI_low_log_per_um_per_s": low,
                            "speed_HC3_95CI_high_log_per_um_per_s": high,
                            "speed_two_sided_p": pvalue,
                            "high_0p2_vs_low_0p05_percent": 100.0
                            * math.expm1(beta * delta_speed),
                            "high_vs_low_percent_CI_low": 100.0
                            * math.expm1(low * delta_speed),
                            "high_vs_low_percent_CI_high": 100.0
                            * math.expm1(high * delta_speed),
                            "model_R2": fit["r2"],
                            "design_rank": fit["rank"],
                            "parameter_count": fit["parameter_count"],
                            "residual_df": fit["residual_df"],
                            "design_condition": fit["condition"],
                            "claim_status": "sequential_time_model_sensitivity_not_randomized_causal_effect",
                        }
                    )
    return rows


def palindrome_speed_contrasts(
    pair_fits: list[dict],
) -> tuple[list[dict], list[dict]]:
    contrast_rows: list[dict] = []
    test_rows: list[dict] = []
    for baseline in BASELINES:
        for variant in FIT_VARIANTS:
            selected = [
                row
                for row in pair_fits
                if row["baseline_method"] == baseline
                and row["fit_variant"] == variant
                and int(row["block"]) in (3, 4, 5)
            ]
            lookup = {
                (int(row["block"]), float(row["nominal_speed_um_per_s"])): row
                for row in selected
            }
            for parameter, _, _ in PARAMETERS:
                differences: list[float] = []
                for block in (3, 4, 5):
                    low = lookup[(block, 0.05)]
                    high = lookup[(block, 0.2)]
                    low_value = float(low[parameter])
                    high_value = float(high[parameter])
                    difference_log = math.log(high_value / low_value)
                    differences.append(difference_log)
                    contrast_rows.append(
                        {
                            "baseline_method": baseline,
                            "fit_variant": variant,
                            "parameter": parameter,
                            "block": block,
                            "low_speed_um_per_s": float(
                                low["actual_gap_speed_um_per_s"]
                            ),
                            "high_speed_um_per_s": float(
                                high["actual_gap_speed_um_per_s"]
                            ),
                            "low_value": low_value,
                            "high_value": high_value,
                            "high_minus_low": high_value - low_value,
                            "log_high_over_low": difference_log,
                            "high_vs_low_percent": 100.0 * math.expm1(difference_log),
                        }
                    )
                values = np.asarray(differences, dtype=np.float64)
                ttest = stats.ttest_1samp(values, 0.0)
                positive = int(np.count_nonzero(values > 0.0))
                negative = int(np.count_nonzero(values < 0.0))
                nonzero = positive + negative
                sign_p = (
                    float(
                        stats.binomtest(
                            min(positive, negative), nonzero, 0.5, alternative="two-sided"
                        ).pvalue
                    )
                    if nonzero
                    else 1.0
                )
                test_rows.append(
                    {
                        "analysis": "palindrome_Fsym_0p2_vs_0p05",
                        "dataset_subset": "no_refresh_testB_blocks3_5",
                        "baseline_method": baseline,
                        "fit_variant": variant,
                        "parameter": parameter,
                        "block_n": 3,
                        "mean_log_high_over_low": float(np.mean(values)),
                        "geometric_high_vs_low_percent": 100.0
                        * math.expm1(float(np.mean(values))),
                        "paired_t_statistic": float(ttest.statistic),
                        "paired_t_two_sided_p": float(ttest.pvalue),
                        "positive_blocks": positive,
                        "negative_blocks": negative,
                        "exact_sign_two_sided_p": sign_p,
                        "claim_status": "n3_block_level_apparent_parameter_test",
                    }
                )
    return contrast_rows, test_rows


def parameter_summaries(fit_rows: list[dict]) -> list[dict]:
    summaries: list[dict] = []
    keys: dict[tuple, list[dict]] = defaultdict(list)
    for row in fit_rows:
        keys[
            (
                row["dataset_type"],
                row["baseline_method"],
                row["fit_variant"],
                float(row["nominal_speed_um_per_s"]),
            )
        ].append(row)
    for (dataset, baseline, variant, speed), rows in sorted(keys.items()):
        for parameter, _, unit in PARAMETERS:
            values = np.asarray([float(row[parameter]) for row in rows])
            summaries.append(
                {
                    "dataset_type": dataset,
                    "baseline_method": baseline,
                    "fit_variant": variant,
                    "nominal_speed_um_per_s": speed,
                    "parameter": parameter,
                    "unit": unit,
                    "n": len(values),
                    "median": float(np.median(values)),
                    "q25": float(np.quantile(values, 0.25)),
                    "q75": float(np.quantile(values, 0.75)),
                    "minimum": float(np.min(values)),
                    "maximum": float(np.max(values)),
                }
            )
    return summaries


def make_map_chronology_figure(map_fits: list[dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14.0, 9.5), sharex=True)
    for column, baseline in enumerate(BASELINES):
        selected = sorted(
            [
                row
                for row in map_fits
                if row["baseline_method"] == baseline
                and row["fit_variant"] == PRIMARY_VARIANT
            ],
            key=lambda row: row["elapsed_min"],
        )
        elapsed_hours = np.asarray([float(row["elapsed_min"]) / 60.0 for row in selected])
        for row_index, (parameter, label, unit) in enumerate(PARAMETERS):
            ax = axes[row_index, column]
            values = np.asarray([float(row[parameter]) for row in selected])
            ax.plot(elapsed_hours, values, color="0.72", lw=1.0, zorder=1)
            for speed in (0.05, 0.1, 0.2):
                mask = np.asarray(
                    [np.isclose(float(row["nominal_speed_um_per_s"]), speed) for row in selected]
                )
                ax.scatter(
                    elapsed_hours[mask],
                    values[mask],
                    s=46,
                    color=SPEED_COLORS[speed],
                    edgecolor="white",
                    linewidth=0.5,
                    label=f"{speed:g} µm/s",
                    zorder=3,
                )
            refresh_candidates = [
                float(row["elapsed_min"]) / 60.0
                for row in selected
                if int(row["block"]) == 3
            ]
            ax.axvline(min(refresh_candidates), color="0.25", ls="--", lw=0.9)
            ax.set_ylabel(f"{label} ({unit})")
            ax.set_title(baseline.replace("_", " "))
            ax.grid(alpha=0.2)
            if row_index == 0:
                ax.legend(frameon=False, fontsize=8, ncol=3)
    for ax in axes[-1]:
        ax.set_xlabel("Elapsed protocol-midpoint time (h)")
    fig.suptitle(
        "Per-map apparent PB parameters versus acquisition time\n"
        "dashed line: first map of shared block 3 after the refresh-affected stage"
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "map_apparent_parameters_chronology.png", dpi=220)
    plt.close(fig)


def make_palindrome_speed_figure(pair_fits: list[dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.5), sharex=True)
    for column, baseline in enumerate(BASELINES):
        selected = [
            row
            for row in pair_fits
            if row["baseline_method"] == baseline
            and row["fit_variant"] == PRIMARY_VARIANT
        ]
        for row_index, (parameter, label, unit) in enumerate(PARAMETERS):
            ax = axes[row_index, column]
            for block in range(1, 6):
                block_rows = sorted(
                    [row for row in selected if int(row["block"]) == block],
                    key=lambda row: row["actual_gap_speed_um_per_s"],
                )
                x = np.asarray(
                    [float(row["actual_gap_speed_um_per_s"]) for row in block_rows]
                )
                y = np.asarray([float(row[parameter]) for row in block_rows])
                ax.plot(
                    x,
                    y,
                    marker="o",
                    ms=5,
                    lw=1.5,
                    color=BLOCK_COLORS[block],
                    label=f"map{block}",
                )
            ax.set_ylabel(f"{label} ({unit})")
            ax.set_title(baseline.replace("_", " "))
            ax.grid(alpha=0.2)
            if row_index == 0:
                ax.legend(frameon=False, fontsize=8, ncol=3)
    for ax in axes[-1]:
        ax.set_xlabel("Actual approach gap speed (µm/s)")
    fig.suptitle("PB fits to palindrome mean Fsym curves: within-block speed patterns")
    fig.tight_layout()
    fig.savefig(FIGURES / "palindrome_apparent_parameters_by_speed.png", dpi=220)
    plt.close(fig)


def make_example_fit_figure(prediction_rows: list[dict]) -> None:
    orders = (1, 15, 29)
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.9), sharex=True, sharey=True)
    for ax, order in zip(axes, orders, strict=True):
        rows = sorted(
            [
                row
                for row in prediction_rows
                if row["dataset_type"] == "map"
                and row["baseline_method"] == "linear_drift_corrected"
                and int(row["acquisition_order"]) == order
            ],
            key=lambda row: row["distance_nm"],
        )
        distance = np.asarray([float(row["distance_nm"]) for row in rows])
        observed = np.asarray([float(row["observed_force_pN"]) for row in rows])
        prediction = np.asarray([float(row["fitted_total_force_pN"]) for row in rows])
        sigma = np.asarray([float(row["spatial_sigma_pN"]) for row in rows])
        ax.fill_between(
            distance,
            observed - sigma,
            observed + sigma,
            color="0.8",
            alpha=0.55,
            label="pixel IQR/1.349 scale",
        )
        ax.plot(distance, observed, color="#264653", lw=1.6, label="measured median")
        ax.plot(distance, prediction, color="#e76f51", lw=2.0, label="PB+vdW fit")
        speed = float(rows[0]["nominal_speed_um_per_s"])
        block = int(rows[0]["block"])
        ax.set_title(f"order {order}, map{block}, {speed:g} µm/s")
        ax.set_yscale("symlog", linthresh=10.0)
        ax.grid(alpha=0.2)
        ax.set_xlabel("Separation D (nm)")
    axes[0].set_ylabel("Force (pN, symlog)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Representative deliberately uncorrected apparent PB fits")
    fig.tight_layout()
    fig.savefig(FIGURES / "representative_apparent_pb_fits.png", dpi=220)
    plt.close(fig)


def make_time_model_speed_figure(speed_rows: list[dict]) -> None:
    models = [
        "global_time_poly1",
        "global_time_poly2",
        "global_time_poly3",
        "block_fixed_linear_within",
        "block_fixed_quadratic_within",
    ]
    fig, axes = plt.subplots(2, 2, figsize=(14.0, 9.5), sharex=True)
    for column, baseline in enumerate(BASELINES):
        for row_index, (parameter, label, _) in enumerate(PARAMETERS):
            ax = axes[row_index, column]
            rows = {
                row["time_model"]: row
                for row in speed_rows
                if row["baseline_method"] == baseline
                and row["fit_variant"] == PRIMARY_VARIANT
                and row["parameter"] == parameter
            }
            center = np.asarray(
                [float(rows[model]["high_0p2_vs_low_0p05_percent"]) for model in models]
            )
            low = np.asarray(
                [float(rows[model]["high_vs_low_percent_CI_low"]) for model in models]
            )
            high = np.asarray(
                [float(rows[model]["high_vs_low_percent_CI_high"]) for model in models]
            )
            x = np.arange(len(models))
            ax.errorbar(
                x,
                center,
                yerr=np.vstack([center - low, high - center]),
                fmt="o",
                color="#264653",
                ecolor="#8d99ae",
                capsize=3,
            )
            ax.axhline(0.0, color="0.35", lw=0.8)
            ax.set_ylabel(f"0.2 vs 0.05 change in {label} (%)")
            ax.set_title(baseline.replace("_", " "))
            ax.set_xticks(x, [model.replace("_", "\n") for model in models], rotation=20)
            ax.grid(alpha=0.2)
    fig.suptitle("Map-level speed coefficient after alternative time models (HC3 95% CI)")
    fig.tight_layout()
    fig.savefig(FIGURES / "time_adjusted_speed_effect_sensitivity.png", dpi=220)
    plt.close(fig)


def fmt(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}g}" if np.isfinite(value) else "n/a"


def build_report(
    fit_rows: list[dict],
    time_rows: list[dict],
    speed_rows: list[dict],
    contrast_rows: list[dict],
    contrast_tests: list[dict],
    model_check: dict[str, float],
) -> str:
    map_fits = [row for row in fit_rows if row["dataset_type"] == "map"]
    pair_fits = [row for row in fit_rows if row["dataset_type"] == "palindrome_Fsym"]
    primary_testb = sorted(
        [
            row
            for row in map_fits
            if row["baseline_method"] == "linear_drift_corrected"
            and row["fit_variant"] == PRIMARY_VARIANT
            and int(row["block"]) in (3, 4, 5)
        ],
        key=lambda row: row["acquisition_order"],
    )
    testb_first = primary_testb[0]
    testb_last = primary_testb[-1]
    testb_duration_min = float(testb_last["elapsed_min"]) - float(
        testb_first["elapsed_min"]
    )
    lambda_drop_percent = 100.0 * (
        float(testb_last["lambda_D_nm"]) / float(testb_first["lambda_D_nm"]) - 1.0
    )
    zeta_drop_percent = 100.0 * (
        float(testb_last["zeta_magnitude_mV"])
        / float(testb_first["zeta_magnitude_mV"])
        - 1.0
    )
    lines = [
        "# 27-08-26：不做hyd/time修正的apparent PB参数粗拟合",
        "",
        "## 直接结论",
        "",
        "这里按要求暂时忽略hydrodynamic force、acquisition-history relaxation和contact-state变化，把每条finite-speed F–D曲线直接当成equilibrium curve拟合。结果中的Debye length和surface potential因此是 **apparent/model-conditioned parameters**，不是已经测得的bulk Debye length或平衡surface potential。",
        "",
        "**直接回答：存在非常强的apparent时间依赖；当前没有可辨识的速度依赖。**",
        "",
        f"在无换液Test B中，line-corrected primary fit从acquisition #{int(testb_first['acquisition_order'])}到#{int(testb_last['acquisition_order'])}（{testb_duration_min:.1f} min），apparent λD由 {float(testb_first['lambda_D_nm']):.2f} 降到 {float(testb_last['lambda_D_nm']):.2f} nm（{lambda_drop_percent:.1f}%），|ψ|由 {float(testb_first['zeta_magnitude_mV']):.1f} 降到 {float(testb_last['zeta_magnitude_mV']):.1f} mV（{zeta_drop_percent:.1f}%）。这个趋势在两种far-field baseline和10–250/20–200 nm两个window下均保持。",
        "",
        "相反，Test B的五种time-adjusted map模型中没有一个速度系数的HC3 95% CI排除零；回文Fsym的0.2−0.05 µm/s参数差在map3/map4/map5之间变号，n=3 paired t-test也不显著。因此速度对apparent PB参数的独立影响至多是当前误差下的几个百分点。",
        "",
    ]

    for baseline in BASELINES:
        selected = sorted(
            [
                row
                for row in map_fits
                if row["baseline_method"] == baseline
                and row["fit_variant"] == PRIMARY_VARIANT
            ],
            key=lambda row: row["acquisition_order"],
        )
        first = selected[0]
        last = selected[-1]
        lines.append(
            f"- `{baseline}`：acquisition #1→#29的apparent λD为 "
            f"{first['lambda_D_nm']:.2f}→{last['lambda_D_nm']:.2f} nm，"
            f"|ψ|为 {first['zeta_magnitude_mV']:.1f}→{last['zeta_magnitude_mV']:.1f} mV。"
        )
    lines.extend(
        [
            "",
            "下面把时间趋势、time-adjusted speed coefficient以及回文Fsym的0.2−0.05 µm/s block contrast分开报告。若速度结果随baseline或time model变号，就不能称为稳定速度依赖。",
            "",
            "## 采用的旧公式",
            "",
            "同材料silica sphere–plane、equal constant-potential nonlinear PB Derjaguin模型：",
            "",
            "`F(D)=2πR ε (kBT/e)^2 κ G(κD,e|ψ|/kBT) − A_H R/(6D²) + b`。",
            "",
            f"固定 `R={base.PROBE_RADIUS_M*1e6:.6f} µm`, `A_H={base.HAMAKER_J:.3g} J`, `εr={EPSILON_R_WATER}`, `T={base.TEMPERATURE_C:.1f} °C`；每条曲线只自由拟合 `λD`, `|ψ|` 和常数baseline `b`。silica按负号约定输出 `ψ=-|ψ|`，但force本身不能从该同表面模型判定电势符号。primary window沿用之前的10–250 nm，同时保留20–200 nm sensitivity branch。",
            "",
            "## 时间和速度检验",
            "",
            "| baseline | parameter | all-29 time Spearman ρ / p | Test-B time ρ / p | time-model 0.2 vs 0.05 effect range | HC3 CI excludes 0 | palindrome blocks 3/4/5 high-vs-low | block t p |",
            "|:---|:---|:---|:---|:---|---:|:---|---:|",
        ]
    )
    for baseline in BASELINES:
        for parameter, label, _ in PARAMETERS:
            all_time = next(
                row
                for row in time_rows
                if row["baseline_method"] == baseline
                and row["fit_variant"] == PRIMARY_VARIANT
                and row["parameter"] == parameter
                and row["dataset_subset"] == "all_29_maps"
            )
            testb_time = next(
                row
                for row in time_rows
                if row["baseline_method"] == baseline
                and row["fit_variant"] == PRIMARY_VARIANT
                and row["parameter"] == parameter
                and row["dataset_subset"] == "no_refresh_testB_blocks3_5"
            )
            models = [
                row
                for row in speed_rows
                if row["baseline_method"] == baseline
                and row["fit_variant"] == PRIMARY_VARIANT
                and row["parameter"] == parameter
            ]
            effects = [float(row["high_0p2_vs_low_0p05_percent"]) for row in models]
            excludes = sum(
                float(row["high_vs_low_percent_CI_low"]) > 0.0
                or float(row["high_vs_low_percent_CI_high"]) < 0.0
                for row in models
            )
            contrasts = sorted(
                [
                    row
                    for row in contrast_rows
                    if row["baseline_method"] == baseline
                    and row["fit_variant"] == PRIMARY_VARIANT
                    and row["parameter"] == parameter
                ],
                key=lambda row: row["block"],
            )
            test = next(
                row
                for row in contrast_tests
                if row["baseline_method"] == baseline
                and row["fit_variant"] == PRIMARY_VARIANT
                and row["parameter"] == parameter
            )
            contrast_text = "/".join(
                f"{float(row['high_vs_low_percent']):+.1f}%" for row in contrasts
            )
            lines.append(
                f"| {baseline} | {label} | "
                f"{all_time['spearman_rho']:+.3f} / {fmt(float(all_time['spearman_two_sided_p']))} | "
                f"{testb_time['spearman_rho']:+.3f} / {fmt(float(testb_time['spearman_two_sided_p']))} | "
                f"{min(effects):+.1f}%…{max(effects):+.1f}% | {excludes}/5 | "
                f"{contrast_text} | {float(test['paired_t_two_sided_p']):.3f} |"
            )

    lines.extend(
        [
            "",
            "- time-model range来自Test B的17张map，分别使用global一至三次时间多项式以及block-fixed线性/二次within-block时间模型；它仍是顺序数据的敏感性分析。",
            "- palindrome contrast先拟合每个block/speed的Fsym曲线，再比较0.2与0.05 µm/s；三个百分数依次对应map3/map4/map5，实验重复数只有n=3。",
            "",
            "## 拟合质量与边界",
            "",
            "| dataset | baseline | window | fits | median R² | minimum R² | grid-boundary fits |",
            "|:---|:---|:---|---:|---:|---:|---:|",
        ]
    )
    for dataset in ("map", "palindrome_Fsym"):
        for baseline in BASELINES:
            for variant in FIT_VARIANTS:
                selected = [
                    row
                    for row in fit_rows
                    if row["dataset_type"] == dataset
                    and row["baseline_method"] == baseline
                    and row["fit_variant"] == variant
                ]
                r2 = np.asarray([float(row["r2"]) for row in selected])
                boundary = sum(row["fit_status"] != "interior_grid_solution" for row in selected)
                lines.append(
                    f"| {dataset} | {baseline} | {variant} | {len(selected)} | "
                    f"{np.median(r2):.4f} | {np.min(r2):.4f} | {boundary} |"
                )

    lines.extend(
        [
            "",
            "每个distance bin的权重尺度为 `max(pixel IQR/1.349, 2 pN)`，是空间离散度而不是mean的standard error；相邻D bins也高度相关。因此CSV中的profile ranges只能说明grid局部可辨识性，不能当正式95% CI。",
            "",
            "## 解释边界",
            "",
            "- 这个计算回答的是：如果强行把每条曲线解释成同一个PB+vdW equilibrium模型，拟合参数怎样随测量变化。它不回答真实λD或ψ是否随时间改变。",
            "- 同一真实表面在约数小时内出现大幅apparent参数变化，首先说明未建模force relaxation被PB参数吸收；不能优先解释成纯水离子强度或silica化学真的同步改变。",
            "- 速度依赖只有在不同block、baseline、fit window和time model下同号且量级稳定时才可信；当前表格直接显示这一稳定性。",
            "- `|ψ|`与`λD`在有限窗口内相关，且同材料force不识别ψ符号；负号仅来自silica convention。",
            "",
            "## 数值检查",
            "",
            f"- 共享PB网格与原 `total_equilibrium_force_pN` 的最大绝对差为 `{model_check['shared_grid_formula_max_abs_error_pN']:.3e} pN`，最大相对差为 `{model_check['shared_grid_formula_max_relative_error']:.3e}`。",
            f"- 网格为 {LAMBDA_GRID_NM.size}个λ点（{LAMBDA_GRID_NM[0]:g}–{LAMBDA_GRID_NM[-1]:g} nm）× {ZETA_GRID_MV.size}个|ψ|点（{ZETA_GRID_MV[0]:g}–{ZETA_GRID_MV[-1]:g} mV）；所有模型值、fit scales、RSS和输出参数均检查finite。",
            "- nonlinear PB evaluator沿用原脚本已经验证的small-potential linear limit和far-field asymptote；本脚本另做公式逐点identity检查。",
            "",
            "## 输出",
            "",
            "- `apparent_pb_fits.csv`: 29张map与14条palindrome Fsym曲线的两种baseline、两个window拟合。",
            "- `apparent_pb_fit_predictions.csv`: primary 10–250 nm曲线、EDL/vdW/baseline分量和residual。",
            "- `parameter_time_tests.csv`, `parameter_time_adjusted_speed_models.csv`: 时间相关与time-adjusted速度系数。",
            "- `palindrome_parameter_contrasts.csv`, `palindrome_parameter_tests.csv`: n=3 block-level 0.2−0.05 µm/s比较。",
            "- `parameter_summaries.csv`, `figures/`, `provenance.json`, `artifact_manifest.sha256`: 汇总、图、参数和身份记录。",
            "",
        ]
    )
    return "\n".join(lines)


def create_manifest(paths: list[Path], destination: Path) -> None:
    rows = [
        f"{sha256_file(path)}  {path.relative_to(ROOT)}"
        for path in sorted(paths, key=lambda item: str(item.relative_to(ROOT)))
    ]
    destination.write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    model_force, flat_lambda, flat_zeta, model_check = model_library()
    inventory, map_curve_rows, map_curves, pair_curves = load_curves()
    map_fits, map_predictions = fit_all_curves(
        map_curves, model_force, flat_lambda, flat_zeta
    )
    pair_fits, pair_predictions = fit_all_curves(
        pair_curves, model_force, flat_lambda, flat_zeta
    )
    fit_rows = map_fits + pair_fits
    prediction_rows = map_predictions + pair_predictions
    time_rows = time_dependence_rows(map_fits)
    speed_rows = time_adjusted_speed_rows(map_fits)
    contrast_rows, contrast_tests = palindrome_speed_contrasts(pair_fits)
    summary_rows = parameter_summaries(fit_rows)

    outputs = {
        "apparent_pb_fits.csv": fit_rows,
        "apparent_pb_fit_predictions.csv": prediction_rows,
        "parameter_time_tests.csv": time_rows,
        "parameter_time_adjusted_speed_models.csv": speed_rows,
        "palindrome_parameter_contrasts.csv": contrast_rows,
        "palindrome_parameter_tests.csv": contrast_tests,
        "parameter_summaries.csv": summary_rows,
    }
    for filename, rows in outputs.items():
        write_csv(RESULTS / filename, rows)

    make_map_chronology_figure(map_fits)
    make_palindrome_speed_figure(pair_fits)
    make_example_fit_figure(prediction_rows)
    make_time_model_speed_figure(speed_rows)

    report = build_report(
        fit_rows,
        time_rows,
        speed_rows,
        contrast_rows,
        contrast_tests,
        model_check,
    )
    (RESULTS / "REPORT.md").write_text(report, encoding="utf-8")

    upstream_paths = [
        UPSTREAM / "map_inventory_QC.csv",
        UPSTREAM / "map_force_curves.csv",
        UPSTREAM / "palindrome_pair_curves.csv",
        UPSTREAM / "provenance.json",
        UPSTREAM / "artifact_manifest.sha256",
        ROOT / "analysis" / "fit_glycerol_surface_forces.py",
        Path(__file__).resolve(),
    ]
    provenance = {
        "analysis": "brute per-map and palindrome-mean apparent PB parameter fits",
        "claim_status": "apparent_model_conditioned_finite_speed_parameters",
        "explicitly_ignored": [
            "hydrodynamic force",
            "acquisition-history relaxation",
            "contact-state dependence",
        ],
        "formula_source": "analysis/fit_glycerol_surface_forces.py",
        "model": MODEL,
        "temperature_C": base.TEMPERATURE_C,
        "epsilon_r": EPSILON_R_WATER,
        "probe_radius_m": base.PROBE_RADIUS_M,
        "hamaker_J": base.HAMAKER_J,
        "surface_geometry": "equal-potential silica sphere-plane",
        "surface_potential_sign": "force identifies magnitude only; negative assigned by silica convention",
        "fit_variants_nm": FIT_VARIANTS,
        "lambda_grid_nm": LAMBDA_GRID_NM.tolist(),
        "zeta_grid_mV": ZETA_GRID_MV.tolist(),
        "baseline_bounds_pN": BASELINE_BOUNDS_PN,
        "minimum_spatial_sigma_pN": MIN_SPATIAL_SIGMA_PN,
        "model_identity_check": model_check,
        "input_counts": {
            "maps": len(inventory),
            "map_force_rows": len(map_curve_rows),
            "map_fit_rows": len(map_fits),
            "palindrome_fit_rows": len(pair_fits),
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256_file(path) for path in upstream_paths
        },
    }
    (RESULTS / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    artifacts = [
        RESULTS / "REPORT.md",
        RESULTS / "provenance.json",
        *[RESULTS / filename for filename in outputs],
        RESULTS / "figures" / "map_apparent_parameters_chronology.png",
        RESULTS / "figures" / "palindrome_apparent_parameters_by_speed.png",
        RESULTS / "figures" / "representative_apparent_pb_fits.png",
        RESULTS / "figures" / "time_adjusted_speed_effect_sensitivity.png",
        Path(__file__).resolve(),
    ]
    create_manifest(artifacts, RESULTS / "artifact_manifest.sha256")
    print(f"Wrote {RESULTS}")
    print(f"map fits={len(map_fits)}, palindrome fits={len(pair_fits)}")
    print(
        "PB shared-grid max abs error = "
        f"{model_check['shared_grid_formula_max_abs_error_pN']:.3e} pN"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Subtract the empirical velocity-distance interaction and repeat PB fits.

The saved force table contains the full fitted K*v**alpha/(D+D0) subtraction.
For numerical PB fitting, each curve is shifted by the model value at 250 nm;
the exact constant-gauge transform leaves lambda and potential unchanged and
is reversed for the reported full-curve offset.  A hydrodynamic mechanism is
not uniquely identified by this operation.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import scipy
from scipy.optimize import least_squares
from scipy.stats import spearmanr, t as student_t

import analyze_12_09_26_997glycerol_hydrodynamic_scaling as hyd
import plot_12_09_26_997glycerol_force_pb_time_series as pb


ROOT = Path(__file__).resolve().parents[1]
FORCE_RESULTS = ROOT / "analysis" / "glycerol_99p7_D3_force_pb_time_results"
TIME_RESULTS = ROOT / "analysis" / "glycerol_99p7_D3_time_speed_results"
HYD_RESULTS = ROOT / "analysis" / "glycerol_99p7_D3_hydrodynamic_scaling_results"
OUT = ROOT / "analysis" / "glycerol_99p7_D3_empirical_hyd_pb_results"
FIG = OUT / "figures"

SPEEDS = np.asarray((0.1, 0.3, 0.9, 2.7), dtype=np.float64)
BLOCKS = (1, 2, 3, 4)
REFERENCE_DISTANCE_NM = 250.0
PRIMARY_VARIANT = "shared_empirical_20_200"
ORIGINAL_VARIANT = "original_no_subtraction"
CORRECTION_VARIANTS = (
    PRIMARY_VARIANT,
    "shared_empirical_20_250",
    "block_specific_empirical_20_200",
)
VARIANT_LABELS = {
    ORIGINAL_VARIANT: "Original",
    PRIMARY_VARIANT: "Shared full empirical subtraction",
    "shared_empirical_20_250": "Shared empirical 20-250 sensitivity",
    "block_specific_empirical_20_200": "Block-specific sensitivity",
}
PARAMETERS = (
    ("lambda_D_apparent_nm", "Apparent Debye length", "nm"),
    ("surface_potential_magnitude_apparent_mV", "Apparent surface potential", "mV"),
)
SPEED_COLORS = {0.1: "#198f7a", 0.3: "#3977a8", 0.9: "#d39b24", 2.7: "#c9533e"}
BLOCK_COLORS = {1: "#3d5a80", 2: "#2a9d8f", 3: "#d39b24", 4: "#c9533e"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def numeric(value) -> float:
    return float(value)


def meta_for_fit(row: dict[str, str]) -> dict:
    return {
        "acquisition_order": int(row["acquisition_order"]),
        "block": int(row["block"]),
        "position_in_block": int(row["position_in_block"]),
        "map_midpoint_time": row["map_midpoint_time"],
        "elapsed_midpoint_min": float(row["elapsed_midpoint_min"]),
        "speed_um_per_s": float(row["speed_um_per_s"]),
        "source": row["source"],
    }


def pair_force_matrices(
    pair_rows: list[dict[str, str]], minimum_nm: float, maximum_nm: float
) -> tuple[np.ndarray, list[np.ndarray]]:
    distances = np.arange(minimum_nm, maximum_nm + 0.1, 5.0, dtype=np.float64)
    lookup = {
        (int(row["block"]), float(row["distance_nm"]), float(row["speed_um_per_s"])): float(
            row["pair_mean_force_pN"]
        )
        for row in pair_rows
    }
    matrices = []
    for block in BLOCKS:
        matrix = np.asarray(
            [
                [lookup[(block, float(distance), float(speed))] for speed in SPEEDS]
                for distance in distances
            ],
            dtype=np.float64,
        )
        matrices.append(hyd.double_center(matrix))
    return distances, matrices


def shared_hyd_fit(
    pair_rows: list[dict[str, str]], minimum_nm: float, maximum_nm: float, label: str
) -> dict:
    distances, observed = pair_force_matrices(pair_rows, minimum_nm, maximum_nm)

    def predicted(parameters: np.ndarray) -> np.ndarray:
        return hyd.model_matrix(distances, float(parameters[0]), float(parameters[1]), float(parameters[2]))

    def residual(parameters: np.ndarray) -> np.ndarray:
        model = predicted(parameters)
        return np.concatenate([(model - matrix).ravel() for matrix in observed])

    starts = [
        np.asarray((coefficient, exponent, shift), dtype=np.float64)
        for coefficient in (5.0e4, 1.5e5, 4.0e5)
        for exponent in (0.4, 0.7, 1.0, 1.3)
        for shift in (0.0, 40.0, 120.0)
    ]
    lower = np.asarray((0.0, 0.05, -minimum_nm + 0.1))
    upper = np.asarray((1.0e7, 2.0, 500.0))
    solutions = [
        least_squares(
            residual,
            start,
            bounds=(lower, upper),
            method="trf",
            x_scale="jac",
            max_nfev=10000,
        )
        for start in starts
    ]
    solution = min(solutions, key=lambda item: float(np.sum(item.fun**2)))
    coefficient, exponent, shift = (float(value) for value in solution.x)
    model = predicted(solution.x)
    block_nrmse = []
    for matrix in observed:
        block_nrmse.append(
            float(np.sqrt(np.mean((matrix - model) ** 2)) / np.sqrt(np.mean(matrix**2)))
        )
    return {
        "hyd_parameter_set": label,
        "fit_scope": "shared_across_four_blocks",
        "fit_min_nm": minimum_nm,
        "fit_max_nm": maximum_nm,
        "coefficient_pN_nm_per_um_s_to_alpha": coefficient,
        "velocity_exponent_alpha": exponent,
        "distance_shift_D0_nm": shift,
        "pooled_residual_rmse_pN": float(np.sqrt(np.mean(solution.fun**2))),
        "block_normalized_rmse_min": float(np.min(block_nrmse)),
        "block_normalized_rmse_median": float(np.median(block_nrmse)),
        "block_normalized_rmse_max": float(np.max(block_nrmse)),
        "optimizer_success": bool(solution.success),
        "optimizer_nfev": int(solution.nfev),
        "optimizer_message": str(solution.message),
    }


def hyd_parameter_sets(pair_rows: list[dict[str, str]]) -> tuple[list[dict], dict[str, dict[int, dict]]]:
    shared_200 = shared_hyd_fit(pair_rows, 20.0, 200.0, PRIMARY_VARIANT)
    shared_250 = shared_hyd_fit(pair_rows, 20.0, 250.0, "shared_empirical_20_250")
    parameter_rows = [shared_200, shared_250]
    by_variant: dict[str, dict[int, dict]] = {
        PRIMARY_VARIANT: {block: shared_200 for block in BLOCKS},
        "shared_empirical_20_250": {block: shared_250 for block in BLOCKS},
    }

    prior_fits = read_csv(HYD_RESULTS / "interaction_model_fits.csv")
    block_parameters: dict[int, dict] = {}
    for block in BLOCKS:
        candidates = [
            row
            for row in prior_fits
            if row["fit_range"] == "primary_20_200"
            and int(row["block"]) == block
            and row["model"] == "velocity_power_and_shift"
        ]
        if len(candidates) != 1:
            raise RuntimeError(f"Block {block}: missing empirical hydrodynamic parameters")
        source = candidates[0]
        current = {
            "hyd_parameter_set": "block_specific_empirical_20_200",
            "fit_scope": f"block_{block}",
            "fit_min_nm": 20.0,
            "fit_max_nm": 200.0,
            "block": block,
            "coefficient_pN_nm_per_um_s_to_alpha": float(
                source["coefficient_pN_nm_per_um_s_to_alpha"]
            ),
            "velocity_exponent_alpha": float(source["velocity_exponent_alpha"]),
            "distance_shift_D0_nm": float(source["distance_shift_D0_nm"]),
            "pooled_residual_rmse_pN": float(source["residual_rmse_pN"]),
            "block_normalized_rmse_min": float(source["normalized_rmse"]),
            "block_normalized_rmse_median": float(source["normalized_rmse"]),
            "block_normalized_rmse_max": float(source["normalized_rmse"]),
            "optimizer_success": source["optimizer_success"] == "True",
            "optimizer_nfev": int(source["optimizer_nfev"]),
            "optimizer_message": source["optimizer_message"],
        }
        parameter_rows.append(current)
        block_parameters[block] = current
    by_variant["block_specific_empirical_20_200"] = block_parameters
    if not all(bool(row["optimizer_success"]) for row in parameter_rows):
        raise RuntimeError("At least one empirical hydrodynamic parameter fit failed")
    return parameter_rows, by_variant


def hyd_full_model_pN(distance_nm: float, speed_um_per_s: float, parameters: dict) -> float:
    coefficient = float(parameters["coefficient_pN_nm_per_um_s_to_alpha"])
    exponent = float(parameters["velocity_exponent_alpha"])
    shift = float(parameters["distance_shift_D0_nm"])
    distance = float(distance_nm)
    speed = float(speed_um_per_s)
    if distance + shift <= 0.0 or speed <= 0.0:
        raise ValueError("Empirical hydrodynamic correction received an invalid coordinate")
    return float(coefficient * speed**exponent / (distance + shift))


def hyd_shape_correction_pN(distance_nm: float, speed_um_per_s: float, parameters: dict) -> float:
    return hyd_full_model_pN(distance_nm, speed_um_per_s, parameters) - hyd_full_model_pN(
        REFERENCE_DISTANCE_NM, speed_um_per_s, parameters
    )


def corrected_force_rows(
    force_rows: list[dict[str, str]],
    chronology: list[dict[str, str]],
    parameters_by_variant: dict[str, dict[int, dict]],
) -> tuple[list[dict], dict[tuple[str, int], list[dict]]]:
    meta = {int(row["acquisition_order"]): row for row in chronology}
    by_task: dict[tuple[str, int], list[dict]] = {}
    output: list[dict] = []
    for variant in CORRECTION_VARIANTS:
        for row in force_rows:
            order = int(row["acquisition_order"])
            block = int(meta[order]["block"])
            speed = float(meta[order]["speed_um_per_s"])
            distance = float(row["distance_nm"])
            parameters = parameters_by_variant[variant][block]
            full_model = hyd_full_model_pN(distance, speed, parameters)
            reference_model = hyd_full_model_pN(
                REFERENCE_DISTANCE_NM, speed, parameters
            )
            shape_correction = full_model - reference_model
            full_corrected_median = float(row["map_median_force_pN"]) - full_model
            full_corrected_q25 = float(row["map_force_q25_pN"]) - full_model
            full_corrected_q75 = float(row["map_force_q75_pN"]) - full_model
            fit_input_median = full_corrected_median + reference_model
            fit_input_q25 = full_corrected_q25 + reference_model
            fit_input_q75 = full_corrected_q75 + reference_model
            fit_row = dict(row)
            fit_row["map_median_force_pN"] = fit_input_median
            fit_row["map_force_q25_pN"] = fit_input_q25
            fit_row["map_force_q75_pN"] = fit_input_q75
            by_task.setdefault((variant, order), []).append(fit_row)
            output.append(
                {
                    "analysis_variant": variant,
                    "acquisition_order": order,
                    "block": block,
                    "speed_um_per_s": speed,
                    "distance_nm": distance,
                    "available_pixels": int(row["available_pixels"]),
                    "original_map_median_force_pN": float(row["map_median_force_pN"]),
                    "empirical_hyd_full_model_pN": full_model,
                    "empirical_hyd_reference_value_pN": reference_model,
                    "empirical_hyd_shape_correction_pN": shape_correction,
                    "full_model_corrected_map_median_force_pN": full_corrected_median,
                    "full_model_corrected_map_force_q25_pN": full_corrected_q25,
                    "full_model_corrected_map_force_q75_pN": full_corrected_q75,
                    "pb_fit_gauge_centered_map_median_force_pN": fit_input_median,
                    "pb_fit_gauge_centered_map_force_q25_pN": fit_input_q25,
                    "pb_fit_gauge_centered_map_force_q75_pN": fit_input_q75,
                    "map_force_spatial_mad_pN": float(row["map_force_spatial_mad_pN"]),
                    "hyd_reference_distance_nm": REFERENCE_DISTANCE_NM,
                    "hyd_coefficient_pN_nm_per_um_s_to_alpha": float(
                        parameters["coefficient_pN_nm_per_um_s_to_alpha"]
                    ),
                    "hyd_velocity_exponent_alpha": float(parameters["velocity_exponent_alpha"]),
                    "hyd_distance_shift_D0_nm": float(parameters["distance_shift_D0_nm"]),
                }
            )
    expected = len(CORRECTION_VARIANTS) * 32 * 60
    if len(output) != expected or any(len(rows) != 60 for rows in by_task.values()):
        raise RuntimeError("Corrected force inventory is incomplete")
    return output, by_task


def fit_corrected_task(task: tuple[str, dict, list[dict], dict]) -> dict:
    variant, meta, rows, parameters = task
    pb.configure_model()
    fit = pb.fit_map(rows, meta)
    reference_model = hyd_full_model_pN(
        REFERENCE_DISTANCE_NM, float(meta["speed_um_per_s"]), parameters
    )
    fit.update(
        {
            "analysis_variant": variant,
            "hydrodynamic_subtraction_applied": True,
            "hydrodynamic_subtraction_interpretation": "full_empirical_model_subtracted_then_constant_gauge_recentered_for_pb_fit_not_unique_hydrodynamic_identification",
            "hyd_reference_distance_nm": REFERENCE_DISTANCE_NM,
            "pb_fit_constant_gauge_added_pN": reference_model,
            "fit_offset_full_model_corrected_curve_pN": float(fit["fit_offset_pN"])
            - reference_model,
            "hyd_coefficient_pN_nm_per_um_s_to_alpha": float(
                parameters["coefficient_pN_nm_per_um_s_to_alpha"]
            ),
            "hyd_velocity_exponent_alpha": float(parameters["velocity_exponent_alpha"]),
            "hyd_distance_shift_D0_nm": float(parameters["distance_shift_D0_nm"]),
        }
    )
    return fit


def original_fit_rows(path: Path) -> list[dict]:
    rows = read_csv(path)
    for row in rows:
        row.update(
            {
                "analysis_variant": ORIGINAL_VARIANT,
                "hydrodynamic_subtraction_applied": False,
                "hydrodynamic_subtraction_interpretation": "none",
                "hyd_reference_distance_nm": float("nan"),
                "pb_fit_constant_gauge_added_pN": 0.0,
                "fit_offset_full_model_corrected_curve_pN": float(row["fit_offset_pN"]),
                "hyd_coefficient_pN_nm_per_um_s_to_alpha": float("nan"),
                "hyd_velocity_exponent_alpha": float("nan"),
                "hyd_distance_shift_D0_nm": float("nan"),
            }
        )
    return rows


def fit_all_corrected(
    chronology: list[dict[str, str]],
    rows_by_task: dict[tuple[str, int], list[dict]],
    parameters_by_variant: dict[str, dict[int, dict]],
    workers: int,
) -> list[dict]:
    tasks = []
    meta_lookup = {int(row["acquisition_order"]): meta_for_fit(row) for row in chronology}
    for variant in CORRECTION_VARIANTS:
        for order in range(1, 33):
            block = int(meta_lookup[order]["block"])
            tasks.append(
                (
                    variant,
                    meta_lookup[order],
                    rows_by_task[(variant, order)],
                    parameters_by_variant[variant][block],
                )
            )
    if workers == 1:
        fits = [fit_corrected_task(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            fits = list(executor.map(fit_corrected_task, tasks, chunksize=1))
    fits.sort(key=lambda row: (CORRECTION_VARIANTS.index(row["analysis_variant"]), int(row["acquisition_order"])))
    if len(fits) != len(CORRECTION_VARIANTS) * 32:
        raise RuntimeError("Corrected PB fit inventory is incomplete")
    return fits


def finite_spearman(x: np.ndarray, y: np.ndarray) -> float:
    keep = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(keep) < 3 or np.unique(x[keep]).size < 2 or np.unique(y[keep]).size < 2:
        return float("nan")
    return float(spearmanr(x[keep], y[keep]).statistic)


def ols_hc3(design: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    x = np.asarray(design, dtype=np.float64)
    y = np.asarray(values, dtype=np.float64)
    if x.ndim != 2 or y.ndim != 1 or x.shape[0] != y.size or x.shape[0] <= x.shape[1]:
        raise ValueError("Invalid OLS design")
    inverse = np.linalg.inv(x.T @ x)
    coefficients = inverse @ x.T @ y
    residual = y - x @ coefficients
    leverage = np.sum(x * (x @ inverse), axis=1)
    adjusted = residual / np.maximum(1.0 - leverage, np.finfo(float).eps)
    covariance = inverse @ (x.T @ np.diag(adjusted**2) @ x) @ inverse
    return coefficients, covariance, x.shape[0] - x.shape[1]


def fit_quality_usable(row: dict) -> bool:
    return str(row["quality_state"]) != "unusable" and np.isfinite(float(row["lambda_D_apparent_nm"]))


def association_rows(all_fits: list[dict]) -> list[dict]:
    rows: list[dict] = []
    variants = (ORIGINAL_VARIANT, *CORRECTION_VARIANTS)
    for variant in variants:
        selected = sorted(
            (row for row in all_fits if row["analysis_variant"] == variant and fit_quality_usable(row)),
            key=lambda row: int(row["acquisition_order"]),
        )
        times = np.asarray([float(row["elapsed_midpoint_min"]) for row in selected])
        speeds = np.asarray([float(row["speed_um_per_s"]) for row in selected])
        for parameter, label, unit in PARAMETERS:
            values = np.asarray([float(row[parameter]) for row in selected])
            same_speed_rho = []
            for speed in SPEEDS:
                keep = np.isclose(speeds, speed)
                rho = finite_spearman(times[keep], values[keep])
                if np.isfinite(rho):
                    same_speed_rho.append(rho)
            time_hours = times / 60.0
            design = np.column_stack(
                (
                    np.ones(values.size),
                    time_hours,
                    np.isclose(speeds, 0.3),
                    np.isclose(speeds, 0.9),
                    np.isclose(speeds, 2.7),
                )
            ).astype(float)
            coefficients, covariance, dof = ols_hc3(design, np.log(values))
            critical = float(student_t.ppf(0.975, dof))
            slope = float(coefficients[1])
            slope_se = float(math.sqrt(max(0.0, covariance[1, 1])))
            rows.append(
                {
                    "analysis_variant": variant,
                    "parameter": parameter,
                    "parameter_label": label,
                    "unit": unit,
                    "usable_map_count": len(selected),
                    "global_spearman_rho_speed": finite_spearman(speeds, values),
                    "global_spearman_rho_elapsed_time": finite_spearman(times, values),
                    "same_speed_time_rho_median": float(np.median(same_speed_rho)),
                    "same_speed_time_rho_min": float(np.min(same_speed_rho)),
                    "same_speed_time_rho_max": float(np.max(same_speed_rho)),
                    "speed_fixed_effect_log_time_slope_per_hour": slope,
                    "speed_fixed_effect_percent_change_per_hour": 100.0 * math.expm1(slope),
                    "speed_fixed_effect_percent_change_ci95_low": 100.0
                    * math.expm1(slope - critical * slope_se),
                    "speed_fixed_effect_percent_change_ci95_high": 100.0
                    * math.expm1(slope + critical * slope_se),
                    "time_model_note": "OLS log(parameter) on elapsed hours plus speed categorical effects; HC3 interval is descriptive under serial acquisition",
                }
            )
    return rows


def same_speed_trend_rows(all_fits: list[dict]) -> list[dict]:
    rows: list[dict] = []
    for variant in (ORIGINAL_VARIANT, *CORRECTION_VARIANTS):
        for parameter, label, unit in PARAMETERS:
            for speed in SPEEDS:
                selected = sorted(
                    (
                        row
                        for row in all_fits
                        if row["analysis_variant"] == variant
                        and fit_quality_usable(row)
                        and np.isclose(float(row["speed_um_per_s"]), speed)
                    ),
                    key=lambda row: float(row["elapsed_midpoint_min"]),
                )
                times = np.asarray([float(row["elapsed_midpoint_min"]) for row in selected])
                values = np.asarray([float(row[parameter]) for row in selected])
                current = {
                    "analysis_variant": variant,
                    "parameter": parameter,
                    "parameter_label": label,
                    "unit": unit,
                    "speed_um_per_s": float(speed),
                    "usable_map_count": len(selected),
                    "spearman_rho_elapsed_time": finite_spearman(times, values),
                    "first_usable_map_order": int(selected[0]["acquisition_order"]) if selected else "",
                    "last_usable_map_order": int(selected[-1]["acquisition_order"]) if selected else "",
                    "first_usable_value": float(values[0]) if values.size else float("nan"),
                    "last_usable_value": float(values[-1]) if values.size else float("nan"),
                    "median_value": float(np.median(values)) if values.size else float("nan"),
                }
                if values.size >= 3:
                    design = np.column_stack((np.ones(values.size), times / 60.0))
                    coefficients, covariance, dof = ols_hc3(design, np.log(values))
                    critical = float(student_t.ppf(0.975, dof))
                    slope = float(coefficients[1])
                    slope_se = float(math.sqrt(max(0.0, covariance[1, 1])))
                    current.update(
                        {
                            "log_time_slope_per_hour": slope,
                            "percent_change_per_hour": 100.0 * math.expm1(slope),
                            "percent_change_ci95_low": 100.0 * math.expm1(slope - critical * slope_se),
                            "percent_change_ci95_high": 100.0 * math.expm1(slope + critical * slope_se),
                        }
                    )
                else:
                    current.update(
                        {
                            "log_time_slope_per_hour": float("nan"),
                            "percent_change_per_hour": float("nan"),
                            "percent_change_ci95_low": float("nan"),
                            "percent_change_ci95_high": float("nan"),
                        }
                    )
                rows.append(current)
    return rows


def palindrome_parameter_rows(all_fits: list[dict]) -> tuple[list[dict], list[dict]]:
    pair_rows: list[dict] = []
    summary_rows: list[dict] = []
    variants = (ORIGINAL_VARIANT, *CORRECTION_VARIANTS)
    for variant in variants:
        selected = [row for row in all_fits if row["analysis_variant"] == variant]
        for parameter, label, unit in PARAMETERS:
            for block in BLOCKS:
                block_rows = []
                for speed in SPEEDS:
                    maps = sorted(
                        (
                            row
                            for row in selected
                            if int(row["block"]) == block
                            and np.isclose(float(row["speed_um_per_s"]), speed)
                        ),
                        key=lambda row: int(row["position_in_block"]),
                    )
                    if len(maps) != 2 or int(maps[0]["position_in_block"]) + int(maps[1]["position_in_block"]) != 9:
                        raise RuntimeError("PB parameter maps do not form the expected palindrome pair")
                    usable = all(fit_quality_usable(row) for row in maps)
                    values = np.asarray([float(row[parameter]) for row in maps])
                    current = {
                        "analysis_variant": variant,
                        "parameter": parameter,
                        "parameter_label": label,
                        "unit": unit,
                        "block": block,
                        "speed_um_per_s": float(speed),
                        "early_map_order": int(maps[0]["acquisition_order"]),
                        "late_map_order": int(maps[1]["acquisition_order"]),
                        "both_maps_usable": usable,
                        "pair_mean": float(np.mean(values)),
                        "pair_half_difference": float(abs(values[1] - values[0]) / 2.0),
                    }
                    pair_rows.append(current)
                    if usable:
                        block_rows.append(current)
                if len(block_rows) >= 3:
                    speed_values = np.asarray([float(row["speed_um_per_s"]) for row in block_rows])
                    parameter_values = np.asarray([float(row["pair_mean"]) for row in block_rows])
                    low = next(
                        (row for row in block_rows if np.isclose(float(row["speed_um_per_s"]), 0.1)),
                        None,
                    )
                    high = next(
                        (row for row in block_rows if np.isclose(float(row["speed_um_per_s"]), 2.7)),
                        None,
                    )
                    summary_rows.append(
                        {
                            "analysis_variant": variant,
                            "parameter": parameter,
                            "block": block,
                            "usable_speed_pair_count": len(block_rows),
                            "complete_four_speed_pairs": len(block_rows) == 4,
                            "pair_mean_spearman_rho_speed": finite_spearman(speed_values, parameter_values),
                            "pair_mean_2p7_minus_0p1": (
                                float(high["pair_mean"] - low["pair_mean"])
                                if low is not None and high is not None
                                else float("nan")
                            ),
                            "unit": unit,
                        }
                    )
    return pair_rows, summary_rows


def paired_difference_statistics(values: list[float]) -> dict:
    differences = np.asarray(values, dtype=np.float64)
    if differences.size < 2 or not np.all(np.isfinite(differences)):
        raise ValueError("Paired comparison requires at least two finite block differences")
    standard_deviation = float(np.std(differences, ddof=1))
    standard_error = standard_deviation / math.sqrt(differences.size)
    critical = float(student_t.ppf(0.975, differences.size - 1))
    mean = float(np.mean(differences))
    return {
        "block_count": int(differences.size),
        "mean_delta_0p3_minus_0p1": mean,
        "block_delta_standard_deviation": standard_deviation,
        "paired_standard_error": standard_error,
        "descriptive_t95_low": mean - critical * standard_error,
        "descriptive_t95_high": mean + critical * standard_error,
        "absolute_mean_over_paired_standard_error": (
            abs(mean) / standard_error if standard_error > 0.0 else float("nan")
        ),
        "all_block_deltas_same_nonzero_sign": bool(
            np.all(differences > 0.0) or np.all(differences < 0.0)
        ),
        "block_deltas": ";".join(f"{value:.12g}" for value in differences),
    }


def low_speed_parameter_comparison(pair_rows: list[dict]) -> list[dict]:
    output: list[dict] = []
    for parameter, label, unit in PARAMETERS:
        lookup = {
            (int(row["block"]), float(row["speed_um_per_s"])): float(row["pair_mean"])
            for row in pair_rows
            if row["analysis_variant"] == PRIMARY_VARIANT
            and row["parameter"] == parameter
            and bool(row["both_maps_usable"])
            and float(row["speed_um_per_s"]) in (0.1, 0.3)
        }
        blocks = [
            block
            for block in BLOCKS
            if (block, 0.1) in lookup and (block, 0.3) in lookup
        ]
        differences = [lookup[(block, 0.3)] - lookup[(block, 0.1)] for block in blocks]
        output.append(
            {
                "analysis_variant": PRIMARY_VARIANT,
                "parameter": parameter,
                "parameter_label": label,
                "unit": unit,
                "speed_low_um_per_s": 0.1,
                "speed_high_um_per_s": 0.3,
                **paired_difference_statistics(differences),
                "inference_boundary": "four sequential palindrome blocks; descriptive paired SE is not an equivalence margin or pure random-error estimate",
            }
        )
    return output


def low_speed_force_shape_comparison(corrected_rows: list[dict]) -> list[dict]:
    selected = [
        row
        for row in corrected_rows
        if row["analysis_variant"] == PRIMARY_VARIANT
        and float(row["speed_um_per_s"]) in (0.1, 0.3)
    ]
    grouped: dict[tuple[int, float, float], list[float]] = {}
    for row in selected:
        key = (
            int(row["block"]),
            float(row["speed_um_per_s"]),
            float(row["distance_nm"]),
        )
        grouped.setdefault(key, []).append(
            float(row["pb_fit_gauge_centered_map_median_force_pN"])
        )

    speed_differences: dict[tuple[int, float], float] = {}
    for block in BLOCKS:
        for distance in np.arange(20.0, REFERENCE_DISTANCE_NM + 0.1, 5.0):
            low = grouped[(block, 0.1, float(distance))]
            high = grouped[(block, 0.3, float(distance))]
            if len(low) != 2 or len(high) != 2:
                raise RuntimeError("Low-speed force comparison lacks a palindrome map pair")
            speed_differences[(block, float(distance))] = float(
                np.mean(high) - np.mean(low)
            )

    output: list[dict] = []
    for distance in np.arange(20.0, REFERENCE_DISTANCE_NM + 0.1, 5.0):
        differences = [
            speed_differences[(block, float(distance))]
            - speed_differences[(block, REFERENCE_DISTANCE_NM)]
            for block in BLOCKS
        ]
        output.append(
            {
                "analysis_variant": PRIMARY_VARIANT,
                "observable": "constant_gauge_invariant_force_shape",
                "distance_nm": float(distance),
                "reference_distance_nm": REFERENCE_DISTANCE_NM,
                "unit": "pN",
                "speed_low_um_per_s": 0.1,
                "speed_high_um_per_s": 0.3,
                **paired_difference_statistics(differences),
                "delta_definition": "([F0.3(D)-F0.1(D)]-[F0.3(250nm)-F0.1(250nm)]) within each block",
                "inference_boundary": "distance bins are correlated; four sequential blocks are the descriptive repeated units",
            }
        )
    return output


def quality_summary(all_fits: list[dict]) -> list[dict]:
    rows: list[dict] = []
    for variant in (ORIGINAL_VARIANT, *CORRECTION_VARIANTS):
        selected = [row for row in all_fits if row["analysis_variant"] == variant]
        counts = {state: sum(str(row["quality_state"]) == state for row in selected) for state in ("pass", "weak", "unusable")}
        rows.append(
            {
                "analysis_variant": variant,
                "map_count": len(selected),
                "pass_count": counts["pass"],
                "weak_count": counts["weak"],
                "unusable_count": counts["unusable"],
                "median_r2_usable": float(
                    np.median([float(row["r2"]) for row in selected if fit_quality_usable(row)])
                ),
                "median_rmse_pN_usable": float(
                    np.median([float(row["residual_rmse_pN"]) for row in selected if fit_quality_usable(row)])
                ),
            }
        )
    return rows


def pair_force_curves(
    force_rows: list[dict[str, str]], chronology: list[dict[str, str]], primary_parameters: dict[int, dict]
) -> tuple[
    dict[tuple[float, float], np.ndarray],
    dict[tuple[float, float], np.ndarray],
    dict[tuple[float, float], np.ndarray],
]:
    meta = {int(row["acquisition_order"]): row for row in chronology}
    lookup = {
        (int(row["acquisition_order"]), float(row["distance_nm"])): float(row["map_median_force_pN"])
        for row in force_rows
    }
    raw: dict[tuple[float, float], list[float]] = {}
    full_corrected: dict[tuple[float, float], list[float]] = {}
    gauge_centered: dict[tuple[float, float], list[float]] = {}
    for block in BLOCKS:
        for speed in SPEEDS:
            maps = [
                order
                for order, row in meta.items()
                if int(row["block"]) == block and np.isclose(float(row["speed_um_per_s"]), speed)
            ]
            if len(maps) != 2:
                raise RuntimeError("Incomplete palindrome force pair")
            parameters = primary_parameters[block]
            for distance in np.arange(20.0, 250.0 + 0.1, 5.0):
                pair_mean = float(np.mean([lookup[(order, float(distance))] for order in maps]))
                full_model = hyd_full_model_pN(float(distance), float(speed), parameters)
                shape_correction = hyd_shape_correction_pN(
                    float(distance), float(speed), parameters
                )
                raw.setdefault((float(speed), float(distance)), []).append(pair_mean)
                full_corrected.setdefault((float(speed), float(distance)), []).append(
                    pair_mean - full_model
                )
                gauge_centered.setdefault((float(speed), float(distance)), []).append(
                    pair_mean - shape_correction
                )
    return (
        {key: np.asarray(value, dtype=np.float64) for key, value in raw.items()},
        {key: np.asarray(value, dtype=np.float64) for key, value in full_corrected.items()},
        {key: np.asarray(value, dtype=np.float64) for key, value in gauge_centered.items()},
    )


def force_comparison_figure(raw: dict, full_corrected: dict, gauge_centered: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.1), sharex=True)
    distances = np.arange(20.0, 250.0 + 0.1, 5.0)
    for axis, values, title in zip(
        axes,
        (raw, full_corrected, gauge_centered),
        (
            "Original baseline-corrected force",
            "Full empirical model subtracted",
            "PB-fit gauge: aligned at 250 nm",
        ),
        strict=True,
    ):
        for speed in SPEEDS:
            matrix = np.asarray([values[(float(speed), float(distance))] for distance in distances]).T
            median = np.median(matrix, axis=0)
            axis.plot(distances, median, color=SPEED_COLORS[float(speed)], lw=2.0, label=f"{speed:g} um/s")
            axis.fill_between(
                distances,
                np.min(matrix, axis=0),
                np.max(matrix, axis=0),
                color=SPEED_COLORS[float(speed)],
                alpha=0.12,
            )
        axis.axhline(0.0, color="0.35", lw=0.8)
        axis.set_title(title)
        axis.set_xlabel("Separation D (nm)")
        axis.set_ylabel("Palindrome pair-mean force (pN)")
        axis.grid(alpha=0.2)
    axes[0].legend(frameon=False)
    fig.suptitle("Empirical force subtraction and the constant-offset gauge", fontsize=14)
    fig.text(
        0.5,
        0.012,
        "Middle: F - K v^alpha/(D+D0). Right: the same curves plus the model value at 250 nm; this constant changes only the fitted offset, not lambda_D or |psi|.",
        ha="center",
        fontsize=9,
        color="0.3",
    )
    fig.tight_layout(rect=(0.02, 0.06, 0.99, 0.94))
    for suffix in ("png", "svg"):
        fig.savefig(FIG / f"force_curves_before_after_empirical_subtraction.{suffix}", dpi=220)
    plt.close(fig)


def parameter_time_figure(all_fits: list[dict]) -> None:
    variants = (ORIGINAL_VARIANT, PRIMARY_VARIANT)
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.4), sharex=True)
    for row_index, (parameter, label, unit) in enumerate(PARAMETERS):
        combined_usable = [
            float(row[parameter])
            for row in all_fits
            if row["analysis_variant"] in variants and fit_quality_usable(row)
        ]
        lower = max(0.05, 0.82 * min(combined_usable))
        upper = 1.22 * max(combined_usable)
        for column, variant in enumerate(variants):
            axis = axes[row_index, column]
            selected = sorted(
                (row for row in all_fits if row["analysis_variant"] == variant),
                key=lambda row: int(row["acquisition_order"]),
            )
            times = np.asarray([float(row["elapsed_midpoint_min"]) for row in selected])
            values = np.asarray([float(row[parameter]) for row in selected])
            states = np.asarray([str(row["quality_state"]) for row in selected])
            usable = states != "unusable"
            axis.plot(times, np.where(usable, values, np.nan), color="0.68", lw=0.9, zorder=1)
            for speed in SPEEDS:
                mask = np.asarray([np.isclose(float(row["speed_um_per_s"]), speed) for row in selected])
                passed = mask & (states == "pass")
                weak = mask & (states == "weak")
                axis.scatter(times[passed], values[passed], s=39, color=SPEED_COLORS[float(speed)], edgecolor="0.25", lw=0.45, zorder=3)
                axis.scatter(times[weak], values[weak], s=49, facecolors="none", edgecolors=SPEED_COLORS[float(speed)], lw=1.6, zorder=4)
            bad = ~usable
            axis.scatter(times[bad], np.full(np.count_nonzero(bad), lower), marker="x", color="#9d2f28", s=43, lw=1.4, zorder=4)
            axis.set_yscale("log")
            axis.set_ylim(lower, upper)
            axis.set_ylabel(f"{label} ({unit})")
            axis.set_title(VARIANT_LABELS[variant])
            axis.grid(alpha=0.2, which="both")
            if row_index == 1:
                axis.set_xlabel("Elapsed time from first map midpoint (min)")
    legend = [
        Line2D([], [], marker="o", color="none", markerfacecolor=color, markeredgecolor="0.25", label=f"{speed:g} um/s")
        for speed, color in SPEED_COLORS.items()
    ]
    legend += [
        Line2D([], [], marker="o", color="none", markerfacecolor="none", markeredgecolor="0.4", label="weak fit"),
        Line2D([], [], marker="x", color="none", markeredgecolor="#9d2f28", label="unusable"),
    ]
    axes[0, 0].legend(handles=legend, frameon=False, ncol=2, fontsize=8)
    fig.suptitle("Apparent PB parameters versus acquisition time", fontsize=15)
    fig.text(
        0.5,
        0.012,
        "Same PB+vdW+offset model and 20-250 nm QC on both sides. Log y axes share limits within each row.",
        ha="center",
        fontsize=9,
        color="0.3",
    )
    fig.tight_layout(rect=(0.02, 0.05, 0.99, 0.94))
    for suffix in ("png", "svg"):
        fig.savefig(FIG / f"pb_parameters_time_original_vs_empirical_corrected.{suffix}", dpi=220)
    plt.close(fig)


def corrected_parameter_time_detail_figure(all_fits: list[dict]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11.0, 8.0), sharex=True)
    for axis, (parameter, label, unit) in zip(axes, PARAMETERS, strict=True):
        selected = sorted(
            (row for row in all_fits if row["analysis_variant"] == PRIMARY_VARIANT),
            key=lambda row: int(row["acquisition_order"]),
        )
        usable_values = np.asarray(
            [float(row[parameter]) for row in selected if fit_quality_usable(row)]
        )
        span = float(np.ptp(usable_values))
        padding = max(0.08 * span, 0.5)
        lower = max(0.0, float(np.min(usable_values)) - padding)
        upper = float(np.max(usable_values)) + padding
        for speed in SPEEDS:
            speed_rows = [
                row
                for row in selected
                if np.isclose(float(row["speed_um_per_s"]), speed)
            ]
            times = np.asarray(
                [float(row["elapsed_midpoint_min"]) for row in speed_rows]
            )
            values = np.asarray([float(row[parameter]) for row in speed_rows])
            states = np.asarray([str(row["quality_state"]) for row in speed_rows])
            usable = states != "unusable"
            axis.plot(
                times,
                np.where(usable, values, np.nan),
                color=SPEED_COLORS[float(speed)],
                lw=1.2,
                alpha=0.72,
            )
            axis.scatter(
                times[states == "pass"],
                values[states == "pass"],
                s=46,
                color=SPEED_COLORS[float(speed)],
                edgecolor="0.25",
                lw=0.5,
                zorder=3,
                label=f"{speed:g} um/s",
            )
            axis.scatter(
                times[states == "weak"],
                values[states == "weak"],
                s=56,
                facecolors="none",
                edgecolors=SPEED_COLORS[float(speed)],
                lw=1.7,
                zorder=4,
            )
            axis.scatter(
                times[states == "unusable"],
                np.full(np.count_nonzero(states == "unusable"), lower),
                marker="x",
                color=SPEED_COLORS[float(speed)],
                s=48,
                lw=1.5,
                zorder=4,
            )
        axis.set_ylim(lower, upper)
        axis.set_ylabel(f"{label} ({unit})")
        axis.grid(alpha=0.2)
    legend = axes[0].get_legend_handles_labels()
    extra = [
        Line2D([], [], marker="o", color="none", markerfacecolor="none", markeredgecolor="0.4", label="weak fit"),
        Line2D([], [], marker="x", color="none", markeredgecolor="0.4", label="unusable"),
    ]
    axes[0].legend(
        [*legend[0], *extra],
        [*legend[1], "weak fit", "unusable"],
        frameon=False,
        ncol=3,
        fontsize=8,
    )
    axes[-1].set_xlabel("Elapsed time from first map midpoint (min)")
    fig.suptitle("Shared empirical correction: apparent PB parameters over time", fontsize=15)
    fig.text(
        0.5,
        0.012,
        "Lines connect usable observations at the same nominal speed; they are visual guides, not fitted kinetics.",
        ha="center",
        fontsize=9,
        color="0.3",
    )
    fig.tight_layout(rect=(0.03, 0.05, 0.99, 0.95))
    for suffix in ("png", "svg"):
        fig.savefig(FIG / f"pb_parameters_time_empirical_corrected_detail.{suffix}", dpi=220)
    plt.close(fig)


def parameter_speed_figure(pair_rows: list[dict]) -> None:
    variants = (ORIGINAL_VARIANT, PRIMARY_VARIANT)
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.2), sharex=True)
    for row_index, (parameter, label, unit) in enumerate(PARAMETERS):
        combined = [
            float(row["pair_mean"])
            for row in pair_rows
            if row["analysis_variant"] in variants
            and row["parameter"] == parameter
            and bool(row["both_maps_usable"])
        ]
        lower = max(0.05, 0.82 * min(combined))
        upper = 1.22 * max(combined)
        for column, variant in enumerate(variants):
            axis = axes[row_index, column]
            for block in BLOCKS:
                selected = sorted(
                    (
                        row
                        for row in pair_rows
                        if row["analysis_variant"] == variant
                        and row["parameter"] == parameter
                        and int(row["block"]) == block
                        and bool(row["both_maps_usable"])
                    ),
                    key=lambda row: float(row["speed_um_per_s"]),
                )
                if len(selected) >= 2:
                    speeds = np.asarray([float(row["speed_um_per_s"]) for row in selected])
                    values = np.asarray([float(row["pair_mean"]) for row in selected])
                    axis.plot(speeds, values, marker="o", color=BLOCK_COLORS[block], lw=1.2, alpha=0.78, label=f"block {block}")
                elif len(selected) == 1:
                    axis.scatter(
                        [float(selected[0]["speed_um_per_s"])],
                        [float(selected[0]["pair_mean"])],
                        color=BLOCK_COLORS[block],
                        marker="o",
                        alpha=0.78,
                        label=f"block {block}",
                    )
            median_speeds = []
            median_values = []
            for speed in SPEEDS:
                available = [
                    float(row["pair_mean"])
                    for row in pair_rows
                    if row["analysis_variant"] == variant
                    and row["parameter"] == parameter
                    and bool(row["both_maps_usable"])
                    and np.isclose(float(row["speed_um_per_s"]), speed)
                ]
                if available:
                    median_speeds.append(float(speed))
                    median_values.append(float(np.median(available)))
            if median_values:
                axis.plot(median_speeds, median_values, color="black", marker="D", lw=2.2, label="available-block median")
            axis.set_xscale("log")
            axis.set_yscale("log")
            axis.set_xticks(SPEEDS, [f"{speed:g}" for speed in SPEEDS])
            axis.set_ylim(lower, upper)
            axis.set_ylabel(f"Pair-mean {label.lower()} ({unit})")
            axis.set_title(VARIANT_LABELS[variant])
            axis.grid(alpha=0.2, which="both")
            if row_index == 1:
                axis.set_xlabel("Nominal approach speed (um/s)")
    axes[0, 0].legend(frameon=False, fontsize=8, ncol=2)
    fig.suptitle("Palindrome-controlled apparent PB parameter speed dependence", fontsize=15)
    fig.text(
        0.5,
        0.012,
        "Each point averages a symmetric same-speed map pair. Gaps mark pairs excluded by fit QC; medians use available blocks.",
        ha="center",
        fontsize=9,
        color="0.3",
    )
    fig.tight_layout(rect=(0.02, 0.05, 0.99, 0.94))
    for suffix in ("png", "svg"):
        fig.savefig(FIG / f"pb_parameters_speed_palindrome_original_vs_corrected.{suffix}", dpi=220)
    plt.close(fig)


def write_report(
    parameter_rows: list[dict],
    all_fits: list[dict],
    associations: list[dict],
    same_speed_trends: list[dict],
    palindrome_summary: list[dict],
    quality: list[dict],
) -> None:
    primary_hyd = next(row for row in parameter_rows if row["hyd_parameter_set"] == PRIMARY_VARIANT)
    quality_by_variant = {row["analysis_variant"]: row for row in quality}
    association_lookup = {
        (row["analysis_variant"], row["parameter"]): row for row in associations
    }

    comparison_lines = []
    for parameter, label, unit in PARAMETERS:
        raw = association_lookup[(ORIGINAL_VARIANT, parameter)]
        corrected = association_lookup[(PRIMARY_VARIANT, parameter)]
        speed_rows = [
            row
            for row in palindrome_summary
            if row["analysis_variant"] == PRIMARY_VARIANT and row["parameter"] == parameter
        ]
        if speed_rows:
            speed_rhos = np.asarray(
                [float(row["pair_mean_spearman_rho_speed"]) for row in speed_rows]
            )
            speed_rho_text = (
                f"{np.median(speed_rhos):.3f} "
                f"[{np.min(speed_rhos):.3f}, {np.max(speed_rhos):.3f}] "
                f"(n={len(speed_rhos)})"
            )
        else:
            speed_rho_text = "not available"
        comparison_lines.append(
            "| {label} | {raw_v:.3f} -> {corrected_v:.3f} | {raw_t:.3f} -> {corrected_t:.3f} | {raw_change:+.1f}% -> {corrected_change:+.1f}% | {pair_rho} |".format(
                label=label,
                raw_v=float(raw["global_spearman_rho_speed"]),
                corrected_v=float(corrected["global_spearman_rho_speed"]),
                raw_t=float(raw["same_speed_time_rho_median"]),
                corrected_t=float(corrected["same_speed_time_rho_median"]),
                raw_change=float(raw["speed_fixed_effect_percent_change_per_hour"]),
                corrected_change=float(corrected["speed_fixed_effect_percent_change_per_hour"]),
                pair_rho=speed_rho_text,
            )
        )

    sensitivity_lines = []
    for variant in CORRECTION_VARIANTS:
        lam = association_lookup[(variant, "lambda_D_apparent_nm")]
        potential = association_lookup[
            (variant, "surface_potential_magnitude_apparent_mV")
        ]
        sensitivity_lines.append(
            "| {label} | {lam_speed:+.3f} | {lam_time:+.1f}% | {psi_speed:+.3f} | {psi_time:+.1f}% |".format(
                label=VARIANT_LABELS[variant],
                lam_speed=float(lam["global_spearman_rho_speed"]),
                lam_time=float(lam["speed_fixed_effect_percent_change_per_hour"]),
                psi_speed=float(potential["global_spearman_rho_speed"]),
                psi_time=float(potential["speed_fixed_effect_percent_change_per_hour"]),
            )
        )

    primary_speed_lines = []
    for speed in SPEEDS:
        lambda_row = next(
            row
            for row in same_speed_trends
            if row["analysis_variant"] == PRIMARY_VARIANT
            and row["parameter"] == "lambda_D_apparent_nm"
            and np.isclose(float(row["speed_um_per_s"]), speed)
        )
        potential_row = next(
            row
            for row in same_speed_trends
            if row["analysis_variant"] == PRIMARY_VARIANT
            and row["parameter"] == "surface_potential_magnitude_apparent_mV"
            and np.isclose(float(row["speed_um_per_s"]), speed)
        )
        primary_speed_lines.append(
            "| {speed:g} | {n_lam} | {rho_lam:.3f} | {first_lam:.2f} -> {last_lam:.2f} | {n_psi} | {rho_psi:.3f} | {first_psi:.2f} -> {last_psi:.2f} |".format(
                speed=float(speed),
                n_lam=int(lambda_row["usable_map_count"]),
                rho_lam=float(lambda_row["spearman_rho_elapsed_time"]),
                first_lam=float(lambda_row["first_usable_value"]),
                last_lam=float(lambda_row["last_usable_value"]),
                n_psi=int(potential_row["usable_map_count"]),
                rho_psi=float(potential_row["spearman_rho_elapsed_time"]),
                first_psi=float(potential_row["first_usable_value"]),
                last_psi=float(potential_row["last_usable_value"]),
            )
        )

    quality_lines = []
    for variant in (ORIGINAL_VARIANT, *CORRECTION_VARIANTS):
        row = quality_by_variant[variant]
        quality_lines.append(
            f"| {VARIANT_LABELS[variant]} | {int(row['pass_count'])} | {int(row['weak_count'])} | {int(row['unusable_count'])} | {float(row['median_r2_usable']):.3f} | {float(row['median_rmse_pN_usable']):.1f} |"
        )

    lines = [
        "# 12-09-26 - empirical hydrodynamic subtraction and apparent-PB refit",
        "",
        "## 结论",
        "",
        "**经验式不能被称为完全、唯一地去除了 hydrodynamics。条件于该经验模型，可以逐点减去完整的 `K v^alpha/(D+D0)`；以下 Debye length 和 surface potential 仍只能称为 empirical-hyd-corrected、model-conditioned apparent parameters。**",
        "",
        f"主分支把四个 palindrome blocks 共享拟合为 `K v^alpha/(D+D0)`：`K={float(primary_hyd['coefficient_pN_nm_per_um_s_to_alpha']):.1f} pN nm/(um/s)^alpha`、`alpha={float(primary_hyd['velocity_exponent_alpha']):.4f}`、`D0={float(primary_hyd['distance_shift_D0_nm']):.2f} nm`。输出表实际保存 `F_full=F-K v^alpha/(D+D0)`。由于输入曲线此前已减过有限窗口 constant baseline，完整减法会留下一个 speed-dependent constant；PB 优化前精确加回该模型在 250 nm 的值，拟合后再从 offset 中减回。这个常数变换不改变 lambda_D、|psi|、残差或 R2。",
        "",
        "经验 hyd 项来自同一批 force curves，因此 speed dependence 变弱是部分 in-sample consequence，不能当作对 hydrodynamic mechanism 的独立验证。主分支使用所有 block 共用参数，以免把 block-to-block 时间变化直接吸收到 block-specific hyd amplitude；后者仅作敏感性检查。",
        "",
        "## 参数趋势比较",
        "",
        "global speed rho 使用 quality_state != unusable 的 map；same-speed time rho 先在四个速度内分别计算再取中位数。每小时变化来自 `log(parameter) ~ elapsed_hours + speed categorical effects`，HC3 interval 只作描述。",
        "",
        "| parameter | global rho(speed), original -> corrected | median same-speed rho(time), original -> corrected | speed-adjusted change/h, original -> corrected | corrected palindrome rho(speed), block median |",
        "|---|---:|---:|---:|---:|",
        *comparison_lines,
        "",
        "完整的整体关联和 HC3 描述区间见 `parameter_time_speed_associations.csv`。速度可视化优先使用同 block 的 symmetric pair means，从而压低一阶 acquisition-time drift；校正后并非每个 block 的四速度 pairs 都通过 QC，所以汇总 rho 使用至少三个可用速度的 block。",
        "",
        "Hyd 参数作用域与拟合范围的敏感性不是小修正，尤其会改变 lambda_D 的时间方向：",
        "",
        "| corrected branch | rho_speed lambda | adjusted lambda change/h | rho_speed potential | adjusted potential change/h |",
        "|---|---:|---:|---:|---:|",
        *sensitivity_lines,
        "",
        "主分支逐速度时间趋势如下；first/last 指该速度首末两个 usable fits，不保证对应整个实验的共同端点：",
        "",
        "| speed (um/s) | n lambda | rho_time lambda | lambda first -> last (nm) | n potential | rho_time potential | potential first -> last (mV) |",
        "|---:|---:|---:|---:|---:|---:|---:|",
        *primary_speed_lines,
        "",
        "## Fit quality",
        "",
        "| variant | pass | weak | unusable | median usable R2 | median usable RMSE (pN) |",
        "|---|---:|---:|---:|---:|---:|",
        *quality_lines,
        "",
        "所有 corrected branches 都复用原始 nonlinear PB Derjaguin + vdW + constant offset、20-250 nm fit window、spatial-MAD weights、四起点 optimizer 和相同 quality flags。fit quality 只说明该投影的数值状态，不证明 corrected curve 是 equilibrium PB force。",
        "",
        "## 可解释范围",
        "",
        "1. corrected force curves 可以用来判断：移除经验 interaction 后，原来很强的 apparent-parameter speed ordering 还剩多少，以及同速度随时间的趋势是否保留。",
        "2. 不能把 corrected lambda_D 当作 bulk Debye length、corrected |psi| 当作已识别 equilibrium surface potential；empirical subtraction、contact zero、epsilon_r=42.5、R、Hamaker constant 和 force scale 都是条件。",
        "3. block-specific subtraction 若比 shared subtraction 更大幅消除时间趋势，说明它同时吸收了 history；这不是更完全的 hydrodynamic correction。",
        "4. 真正的 physical subtraction 仍需 raw `D(t)` 和 `U_gap(t)=-dD/dt`，并用独立 viscosity、temperature、blank/control 和 approach-retract sign reversal 检查。当前 raw JPK 不在本地 checkout。",
        "",
        "## 文件",
        "",
        "- `map_force_empirical_hyd_corrected.csv`: 三个 branches 的完整模型减法曲线、250 nm gauge-centered PB 输入曲线及两种 correction 分量。",
        "- `map_pb_fits_comparison.csv`: 原始与三个 corrected branches 的全部 PB fits、full-curve offset 和 QC。",
        "- `empirical_hyd_parameter_sets.csv`: shared 与 block-specific 经验参数。",
        "- `parameter_time_speed_associations.csv`: time/speed rank association 与 speed-adjusted log-time trend。",
        "- `same_speed_parameter_time_trends.csv`: 每个速度的 usable-count、首末值、rho 和 log-time slope。",
        "- `palindrome_parameter_pair_means.csv`, `palindrome_parameter_speed_summary.csv`: block 内速度对照。",
        "- `low_speed_parameter_pair_comparison.csv`: 0.3-0.1 um/s 的四-block PB 参数 paired SE 与描述区间。",
        "- `low_speed_force_shape_comparison.csv`: 消去 250 nm constant gauge 后的逐距离低速 force-shape 对照。",
        "- `fit_quality_summary.csv`: pass/weak/unusable 与 residual summary。",
        "- `figures/force_curves_before_after_empirical_subtraction.*`: 原始、完整模型减法、PB gauge-centered 三联图。",
        "- `figures/pb_parameters_time_original_vs_empirical_corrected.*`: Debye/potential 时间趋势。",
        "- `figures/pb_parameters_time_empirical_corrected_detail.*`: corrected-only 同速度时间趋势放大图。",
        "- `figures/pb_parameters_speed_palindrome_original_vs_corrected.*`: palindrome-controlled 速度比较。",
        "- `provenance.json`, `artifact_manifest.sha256`: 输入 hashes、模型定义、软件和产物校验。",
    ]
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_manifest() -> None:
    targets = sorted(
        path
        for path in OUT.rglob("*")
        if path.is_file() and path.name != "artifact_manifest.sha256"
    )
    lines = [f"{sha256_file(path)}  {path.relative_to(ROOT).as_posix()}" for path in targets]
    (OUT / "artifact_manifest.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument(
        "--postprocess-only",
        action="store_true",
        help="Reuse saved corrected force/PB fits and rebuild summaries, report, figures, and manifest.",
    )
    arguments = parser.parse_args()
    if arguments.workers < 1:
        raise ValueError("--workers must be positive")

    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    force_path = FORCE_RESULTS / "map_force_by_separation.csv"
    chronology_path = TIME_RESULTS / "map_time_speed_observables.csv"
    original_fit_path = FORCE_RESULTS / "map_pb_apparent_fits.csv"
    pair_path = HYD_RESULTS / "palindrome_pair_force_matrix.csv"
    force_rows = read_csv(force_path)
    chronology = sorted(read_csv(chronology_path), key=lambda row: int(row["acquisition_order"]))
    pair_rows = read_csv(pair_path)
    if len(force_rows) != 32 * 60 or len(chronology) != 32:
        raise RuntimeError("Input force or chronology inventory is incomplete")

    parameter_rows, parameters_by_variant = hyd_parameter_sets(pair_rows)
    if arguments.postprocess_only:
        corrected_force_path = OUT / "map_force_empirical_hyd_corrected.csv"
        comparison_fit_path = OUT / "map_pb_fits_comparison.csv"
        if not corrected_force_path.exists() or not comparison_fit_path.exists():
            raise RuntimeError("--postprocess-only requires existing corrected force and PB fit CSVs")
        corrected_rows = read_csv(corrected_force_path)
        all_fits = read_csv(comparison_fit_path)
        corrected_fits = [row for row in all_fits if row["analysis_variant"] in CORRECTION_VARIANTS]
    else:
        corrected_rows, rows_by_task = corrected_force_rows(
            force_rows, chronology, parameters_by_variant
        )
        corrected_fits = fit_all_corrected(
            chronology, rows_by_task, parameters_by_variant, arguments.workers
        )
        original_fits = original_fit_rows(original_fit_path)
        all_fits = [*original_fits, *corrected_fits]
    associations = association_rows(all_fits)
    same_speed_trends = same_speed_trend_rows(all_fits)
    palindrome_pairs, palindrome_summary = palindrome_parameter_rows(all_fits)
    low_speed_parameters = low_speed_parameter_comparison(palindrome_pairs)
    low_speed_force = low_speed_force_shape_comparison(corrected_rows)
    quality = quality_summary(all_fits)

    write_csv(OUT / "empirical_hyd_parameter_sets.csv", parameter_rows)
    if not arguments.postprocess_only:
        write_csv(OUT / "map_force_empirical_hyd_corrected.csv", corrected_rows)
        write_csv(OUT / "map_pb_fits_comparison.csv", all_fits)
    write_csv(OUT / "parameter_time_speed_associations.csv", associations)
    write_csv(OUT / "same_speed_parameter_time_trends.csv", same_speed_trends)
    write_csv(OUT / "palindrome_parameter_pair_means.csv", palindrome_pairs)
    write_csv(OUT / "palindrome_parameter_speed_summary.csv", palindrome_summary)
    write_csv(OUT / "low_speed_parameter_pair_comparison.csv", low_speed_parameters)
    write_csv(OUT / "low_speed_force_shape_comparison.csv", low_speed_force)
    write_csv(OUT / "fit_quality_summary.csv", quality)

    raw_curves, full_corrected_curves, gauge_centered_curves = pair_force_curves(
        force_rows, chronology, parameters_by_variant[PRIMARY_VARIANT]
    )
    force_comparison_figure(raw_curves, full_corrected_curves, gauge_centered_curves)
    parameter_time_figure(all_fits)
    corrected_parameter_time_detail_figure(all_fits)
    parameter_speed_figure(palindrome_pairs)
    write_report(
        parameter_rows,
        all_fits,
        associations,
        same_speed_trends,
        palindrome_summary,
        quality,
    )

    provenance = {
        "script": str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
        "inputs": {
            "force_csv": str(force_path.relative_to(ROOT)).replace("\\", "/"),
            "force_csv_sha256": sha256_file(force_path),
            "chronology_csv": str(chronology_path.relative_to(ROOT)).replace("\\", "/"),
            "chronology_csv_sha256": sha256_file(chronology_path),
            "original_pb_fit_csv": str(original_fit_path.relative_to(ROOT)).replace("\\", "/"),
            "original_pb_fit_csv_sha256": sha256_file(original_fit_path),
            "hyd_pair_matrix_csv": str(pair_path.relative_to(ROOT)).replace("\\", "/"),
            "hyd_pair_matrix_csv_sha256": sha256_file(pair_path),
            "hyd_interaction_fit_csv": str((HYD_RESULTS / "interaction_model_fits.csv").relative_to(ROOT)).replace("\\", "/"),
            "hyd_interaction_fit_csv_sha256": sha256_file(HYD_RESULTS / "interaction_model_fits.csv"),
        },
        "primary_subtraction": {
            "variant": PRIMARY_VARIANT,
            "full_corrected_curve_formula": "F_full_corrected=F_original-K*v**alpha/(D+D0)",
            "pb_fit_input_formula": "F_pb_input=F_full_corrected+K*v**alpha/(250nm+D0)",
            "reference_distance_nm": REFERENCE_DISTANCE_NM,
            "parameter_scope": "shared across four palindrome blocks; fitted on double-centred 20-200 nm interaction",
            "constant_gauge_transform": "the added 250 nm value is removed from the reported full-curve fit offset and cannot change lambda, potential, residuals, or R2",
            "identification_boundary": "complete only conditional on the empirical formula; the data do not uniquely identify all hydrodynamic force",
        },
        "pb_refit": {
            "implementation": "plot_12_09_26_997glycerol_force_pb_time_series.fit_map",
            "model": pb.MODEL,
            "fit_range_nm": [pb.FIT_MIN_NM, pb.FIT_MAX_NM],
            "epsilon_r": pb.EPSILON_R_997_GLYCEROL,
            "probe_radius_m": pb.PROBE_RADIUS_M,
            "hamaker_J": pb.HAMAKER_J,
            "temperature_C": pb.TEMPERATURE_C,
            "optimizer_starts": "unchanged four starts",
            "quality_rules": "unchanged from original 12-09-26 force/PB analysis",
        },
        "parallel_workers": arguments.workers,
        "postprocess_only": arguments.postprocess_only,
        "claim_boundary": "empirical-hyd-corrected apparent PB trends; not complete hydrodynamic removal or equilibrium parameter identification",
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
    }
    (OUT / "provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    write_manifest()
    print(
        json.dumps(
            {
                "output": str(OUT),
                "workers": arguments.workers,
                "corrected_fit_count": len(corrected_fits),
                "quality": quality,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

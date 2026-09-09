#!/usr/bin/env python3
"""Advanced same-pixel analysis of the four 31-08-26 D5 force maps.

The analysis deliberately separates model-light evidence from model-conditioned
surface-force fits.  It uses all 64 original approach curves from each of
map2--map5 over 20--200 nm, keeps the exact physical-pixel pairing, and writes
machine-readable tables, figures, provenance, and a numerical audit report.

Primary physical fit: equal constant-potential silica sphere--silica plane,
nonlinear 1:1 Poisson--Boltzmann Derjaguin force plus sphere--plane van der
Waals, fixed no-slip lubrication force, and a free per-curve constant baseline.
The fitted lambda_D and |psi| remain apparent/model-conditioned parameters.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
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
import pandas as pd
import scipy
from scipy import stats
from scipy.constants import Boltzmann, elementary_charge, epsilon_0
from scipy.optimize import curve_fit


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import analyze_31_08_26_maps as raw  # noqa: E402
import fit_31_08_26_colloid_colloid as colloid  # noqa: E402
import fit_glycerol_surface_forces as base  # noqa: E402


UPSTREAM = ROOT / "analysis" / "maps_31_08_26_results"
COLLOID_UPSTREAM = ROOT / "analysis" / "maps_31_08_26_colloid_colloid_results"
CALIBRATION_SUMMARY = (
    ROOT / "analysis" / "palindrome_27_08_26_full_results" / "calibration_summary.csv"
)
RESULTS = ROOT / "analysis" / "maps_31_08_26_advanced_results"
FIGURES = RESULTS / "figures"

TEMPERATURE_K = base.TEMPERATURE_K
EPSILON_R = 78.5
RADIUS_M = raw.CONTEXT_PROBE_RADIUS_M
EFFECTIVE_EQUAL_SPHERE_RADIUS_M = RADIUS_M / 2.0
HAMAKER_J = base.HAMAKER_J
VISCOSITY_PA_S = base.cheng_viscosity_mPa_s(0.0, raw.CONTEXT_TEMPERATURE_C) * 1e-3
FIT_MIN_NM = 20.0
FIT_MAX_NM = 200.0
LAMBDA_GRID_NM = np.geomspace(1.0, 300.0, 181, dtype=np.float64)
ZETA_GRID_MV = np.linspace(1.0, 250.0, 167, dtype=np.float64)
CONTACT_OFFSETS_NM = np.asarray([-10.0, -5.0, -2.0, 0.0, 2.0, 5.0, 10.0])
NORMALIZATION_MIN_DENOMINATOR_PN = 50.0
PROFILE_DELTA_TWO_PARAMETERS = 2.30
RANDOM_SEED = 20260901
MORAN_PERMUTATIONS = 4999
BOOTSTRAP_SAMPLES = 20000
PRIMARY_MODEL = "sphere_plane_nonlinear_PB_fixed_no_slip_hydrodynamics"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def finite_quantile(values: np.ndarray, q: float) -> float:
    selected = np.asarray(values, dtype=np.float64)
    selected = selected[np.isfinite(selected)]
    return float(np.quantile(selected, q)) if selected.size else float("nan")


def stable_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    valid = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(valid) < 4 or np.ptp(x[valid]) == 0 or np.ptp(y[valid]) == 0:
        return float("nan"), float("nan"), int(np.count_nonzero(valid))
    result = stats.spearmanr(x[valid], y[valid])
    return float(result.statistic), float(result.pvalue), int(np.count_nonzero(valid))


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    values = np.asarray(p_values, dtype=np.float64)
    adjusted = np.full(values.shape, np.nan, dtype=np.float64)
    valid_indices = np.flatnonzero(np.isfinite(values))
    if not valid_indices.size:
        return adjusted.tolist()
    ordered_local = np.argsort(values[valid_indices])
    ordered_indices = valid_indices[ordered_local]
    ordered_p = values[ordered_indices]
    count = ordered_p.size
    raw_adjusted = ordered_p * count / np.arange(1, count + 1)
    monotone = np.minimum.accumulate(raw_adjusted[::-1])[::-1]
    adjusted[ordered_indices] = np.minimum(monotone, 1.0)
    return adjusted.tolist()


def exact_row_signflip_p(difference: np.ndarray, physical_rows: np.ndarray) -> float:
    """Two-sided exact sign-flip p-value over eight physical-row effects."""
    difference = np.asarray(difference, dtype=np.float64)
    physical_rows = np.asarray(physical_rows, dtype=int)
    effects: list[float] = []
    for row in range(8):
        selected = difference[(physical_rows == row) & np.isfinite(difference)]
        if selected.size:
            effects.append(float(np.median(selected)))
    row_effects = np.asarray(effects, dtype=np.float64)
    if row_effects.size < 4 or not np.all(np.isfinite(row_effects)):
        return float("nan")
    observed = abs(float(np.mean(row_effects)))
    statistics = []
    for signs in itertools.product((-1.0, 1.0), repeat=row_effects.size):
        statistics.append(abs(float(np.mean(row_effects * np.asarray(signs)))))
    return float(np.mean(np.asarray(statistics) >= observed - 1e-15))


def row_cluster_bootstrap_ci(
    difference: np.ndarray,
    physical_rows: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float]:
    difference = np.asarray(difference, dtype=np.float64)
    physical_rows = np.asarray(physical_rows, dtype=int)
    row_values = [
        difference[(physical_rows == row) & np.isfinite(difference)] for row in range(8)
    ]
    row_values = [values for values in row_values if values.size]
    if len(row_values) < 2:
        return float("nan"), float("nan")
    estimates = np.empty(BOOTSTRAP_SAMPLES, dtype=np.float64)
    for sample in range(BOOTSTRAP_SAMPLES):
        selected_rows = rng.integers(0, len(row_values), size=len(row_values))
        values = np.concatenate([row_values[index] for index in selected_rows])
        estimates[sample] = np.nanmedian(values)
    return float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))


def physical_matrix(values: np.ndarray, rows: np.ndarray, columns: np.ndarray) -> np.ndarray:
    output = np.full((8, 8), np.nan, dtype=np.float64)
    for value, row, column in zip(values, rows, columns, strict=True):
        output[int(row), int(column)] = float(value)
    return output


def neighbor_pairs() -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    horizontal: list[tuple[int, int]] = []
    vertical: list[tuple[int, int]] = []
    for row in range(8):
        for column in range(8):
            index = row * 8 + column
            if column < 7:
                horizontal.append((index, index + 1))
            if row < 7:
                vertical.append((index, index + 8))
    return horizontal, vertical


HORIZONTAL_PAIRS, VERTICAL_PAIRS = neighbor_pairs()
DIRECTED_NEIGHBORS = [
    pair for edge in HORIZONTAL_PAIRS + VERTICAL_PAIRS for pair in (edge, edge[::-1])
]


def moran_i(values: np.ndarray) -> float:
    vector = np.asarray(values, dtype=np.float64).reshape(-1)
    valid = np.isfinite(vector)
    if vector.size != 64 or np.count_nonzero(valid) < 4:
        return float("nan")
    centered = np.full_like(vector, np.nan)
    centered[valid] = vector[valid] - np.mean(vector[valid])
    denominator = float(np.sum(centered[valid] ** 2))
    if denominator <= 0.0:
        return float("nan")
    valid_neighbors = [
        (i, j) for i, j in DIRECTED_NEIGHBORS if valid[i] and valid[j]
    ]
    if not valid_neighbors:
        return float("nan")
    numerator = sum(centered[i] * centered[j] for i, j in valid_neighbors)
    return float(
        np.count_nonzero(valid) / len(valid_neighbors) * numerator / denominator
    )


def moran_permutation_test(
    values: np.ndarray, rng: np.random.Generator
) -> tuple[float, float]:
    observed = moran_i(values)
    if not np.isfinite(observed):
        return observed, float("nan")
    expected = -1.0 / 63.0
    deviations = np.empty(MORAN_PERMUTATIONS, dtype=np.float64)
    vector = np.asarray(values, dtype=np.float64).reshape(-1)
    valid = np.isfinite(vector)
    expected = -1.0 / (np.count_nonzero(valid) - 1.0)
    for index in range(MORAN_PERMUTATIONS):
        candidate = vector.copy()
        candidate[valid] = rng.permutation(vector[valid])
        deviations[index] = abs(moran_i(candidate) - expected)
    p_value = (1.0 + np.count_nonzero(deviations >= abs(observed - expected))) / (
        MORAN_PERMUTATIONS + 1.0
    )
    return observed, float(p_value)


def pb_sphere_plane_library(distance_nm: np.ndarray) -> np.ndarray:
    """Nonlinear equal-potential 1:1 PB + vdW sphere-plane library, in pN."""
    distance_nm = np.asarray(distance_nm, dtype=np.float64)
    if distance_nm.ndim != 1 or np.any(distance_nm <= 0.0):
        raise ValueError("model separations must be a positive 1-D array")
    distance_m = distance_nm * 1e-9
    lambda_m = LAMBDA_GRID_NM * 1e-9
    kappa = 1.0 / lambda_m
    dimensionless_gap = distance_nm[None, :] / LAMBDA_GRID_NM[:, None]
    absolute_permittivity = EPSILON_R * epsilon_0
    prefactor = (
        2.0
        * np.pi
        * RADIUS_M
        * absolute_permittivity
        * (Boltzmann * TEMPERATURE_K / elementary_charge) ** 2
    )
    vdw_pN = -HAMAKER_J * RADIUS_M / (6.0 * distance_m**2) * 1e12
    library = np.empty(
        (ZETA_GRID_MV.size, LAMBDA_GRID_NM.size, distance_nm.size),
        dtype=np.float64,
    )
    for zeta_index, zeta_mV in enumerate(ZETA_GRID_MV):
        dimensionless_surface = (
            elementary_charge * zeta_mV * 1e-3 / (Boltzmann * TEMPERATURE_K)
        )
        dimensionless_integral = base.pb_dimensionless_G(
            dimensionless_gap, dimensionless_surface
        )
        edl_pN = prefactor * kappa[:, None] * dimensionless_integral * 1e12
        library[zeta_index] = edl_pN + vdw_pN[None, :]
    flattened = library.reshape(-1, distance_nm.size)
    if not np.all(np.isfinite(flattened)):
        raise FloatingPointError("nonlinear PB library contains non-finite values")
    return flattened


def linear_hhf_library(distance_nm: np.ndarray, radius_m: float) -> np.ndarray:
    """Equal-potential linear HHF + vdW library for a Derjaguin radius."""
    distance_nm = np.asarray(distance_nm, dtype=np.float64)
    distance_m = distance_nm * 1e-9
    lambda_m = LAMBDA_GRID_NM * 1e-9
    kappa = 1.0 / lambda_m
    H = distance_nm[None, :] / LAMBDA_GRID_NM[:, None]
    reciprocal = np.where(H > 40.0, np.exp(-H), 1.0 / (np.exp(H) + 1.0))
    zeta_V2 = (ZETA_GRID_MV * 1e-3) ** 2
    amplitude = (
        4.0
        * np.pi
        * radius_m
        * EPSILON_R
        * epsilon_0
        * kappa[:, None]
        * reciprocal
        * 1e12
    )
    edl = zeta_V2[:, None, None] * amplitude[None, :, :]
    vdw = -HAMAKER_J * radius_m / (6.0 * distance_m**2) * 1e12
    library = edl + vdw[None, None, :]
    flattened = library.reshape(-1, distance_nm.size)
    if not np.all(np.isfinite(flattened)):
        raise FloatingPointError("linear HHF library contains non-finite values")
    return flattened


FLAT_LAMBDA_NM = np.tile(LAMBDA_GRID_NM, ZETA_GRID_MV.size)
FLAT_ZETA_MV = np.repeat(ZETA_GRID_MV, LAMBDA_GRID_NM.size)


def hydrodynamic_force_pN(
    distance_nm: np.ndarray, speed_um_per_s: np.ndarray, radius_m: float
) -> np.ndarray:
    distance_m = np.asarray(distance_nm, dtype=np.float64)[None, :] * 1e-9
    speed_m_per_s = np.asarray(speed_um_per_s, dtype=np.float64)[:, None] * 1e-6
    force = (
        6.0
        * np.pi
        * VISCOSITY_PA_S
        * radius_m**2
        * speed_m_per_s
        / distance_m
        * 1e12
    )
    if not np.all(np.isfinite(force)) or np.any(force <= 0.0):
        raise FloatingPointError("invalid fixed hydrodynamic force")
    return force


def batch_grid_fit(
    observed_pN: np.ndarray,
    equilibrium_library_pN: np.ndarray,
    hydrodynamic_pN: np.ndarray,
    keep_scores: bool = False,
) -> dict[str, np.ndarray]:
    """OLS grid fit with an analytical free baseline for every curve."""
    observed = np.asarray(observed_pN, dtype=np.float64)
    library = np.asarray(equilibrium_library_pN, dtype=np.float64)
    hydrodynamic = np.asarray(hydrodynamic_pN, dtype=np.float64)
    if observed.ndim != 2 or library.ndim != 2 or hydrodynamic.shape != observed.shape:
        raise ValueError("batch fit shape contract failed")
    if library.shape[1] != observed.shape[1]:
        raise ValueError("model and observation distance axes differ")
    if not all(np.all(np.isfinite(value)) for value in (observed, library, hydrodynamic)):
        raise FloatingPointError("batch fit input is non-finite")

    adjusted = observed - hydrodynamic
    model_mean = np.mean(library, axis=1)
    data_mean = np.mean(adjusted, axis=1)
    centered_model = library - model_mean[:, None]
    centered_data = adjusted - data_mean[:, None]
    model_norm = np.einsum("ij,ij->i", centered_model, centered_model)
    data_norm = np.einsum("ij,ij->i", centered_data, centered_data)
    scores = (
        model_norm[:, None]
        + data_norm[None, :]
        - 2.0 * centered_model @ centered_data.T
    )
    scores = np.maximum(scores, 0.0)
    best_index = np.argmin(scores, axis=0)
    curve_index = np.arange(observed.shape[0])
    baseline = data_mean - model_mean[best_index]
    prediction = library[best_index] + hydrodynamic + baseline[:, None]
    residual = observed - prediction
    direct_rss = np.einsum("ij,ij->i", residual, residual)
    selected_rss = scores[best_index, curve_index]
    if not np.allclose(direct_rss, selected_rss, rtol=2e-10, atol=1e-7):
        raise AssertionError("batch score and direct residual closure failed")
    if np.any(np.abs(baseline) > 500.0):
        raise RuntimeError("selected baseline exceeds the documented +/-500 pN range")
    total = np.sum((observed - np.mean(observed, axis=1)[:, None]) ** 2, axis=1)
    r2 = np.where(total > 0.0, 1.0 - direct_rss / total, np.nan)
    n = observed.shape[1]
    safe_rss = np.maximum(direct_rss, np.finfo(np.float64).tiny)
    result = {
        "best_index": best_index,
        "lambda_D_nm": FLAT_LAMBDA_NM[best_index],
        "zeta_magnitude_mV": FLAT_ZETA_MV[best_index],
        "baseline_pN": baseline,
        "prediction_pN": prediction,
        "residual_pN": residual,
        "rss_pN2": direct_rss,
        "rmse_pN": np.sqrt(direct_rss / n),
        "r2": r2,
        "aic": n * np.log(safe_rss / n) + 2.0 * 3.0,
        "bic": n * np.log(safe_rss / n) + math.log(n) * 3.0,
    }
    if keep_scores:
        result["scores_pN2"] = scores
    return result


def load_inputs() -> dict:
    archive = np.load(UPSTREAM / "pixel_force_curves.npz", allow_pickle=False)
    distance_all = np.asarray(archive["distance_nm"], dtype=np.float64)
    mask = (distance_all >= FIT_MIN_NM) & (distance_all <= FIT_MAX_NM)
    distance_nm = distance_all[mask]
    force = np.asarray(archive["force_linear_drift_corrected_pN"], dtype=np.float64)[
        1:5, :, mask
    ]
    sources = np.asarray(archive["source"])[1:5]
    if distance_nm.size != 37 or force.shape != (4, 64, 37):
        raise RuntimeError(f"unexpected force cube {force.shape}, distance {distance_nm.shape}")
    if not np.all(np.isfinite(force)):
        raise FloatingPointError("force cube contains non-finite values")

    inventory = pd.read_csv(UPSTREAM / "map_inventory_QC.csv").sort_values("map_order")
    inventory = inventory[inventory["map_order"].between(2, 5)].copy()
    if inventory.shape[0] != 4 or list(inventory["source"]) != list(sources):
        raise RuntimeError("NPZ source order and map inventory do not match")
    times = (
        pd.to_datetime(inventory["timestamp"])
        - pd.to_datetime(inventory["timestamp"]).iloc[0]
    ).dt.total_seconds().to_numpy(dtype=np.float64) / 60.0

    qc = pd.read_csv(UPSTREAM / "pixel_QC.csv")
    qc = qc[qc["map_order"].between(2, 5)].sort_values(
        ["map_order", "point_index"]
    )
    if qc.shape[0] != 256:
        raise RuntimeError("expected 256 map2--5 QC rows")
    rows = qc[qc["map_order"] == 2]["row"].to_numpy(dtype=int)
    columns = qc[qc["map_order"] == 2]["column"].to_numpy(dtype=int)
    for map_order in range(3, 6):
        selected = qc[qc["map_order"] == map_order]
        if not np.array_equal(rows, selected["row"].to_numpy(dtype=int)) or not np.array_equal(
            columns, selected["column"].to_numpy(dtype=int)
        ):
            raise RuntimeError("physical pixel mapping differs across map2--5")
    speeds = qc["gap_speed_20_200nm_um_per_s"].to_numpy(dtype=np.float64).reshape(4, 64)
    noise = qc["far_noise_pN"].to_numpy(dtype=np.float64).reshape(4, 64)
    if np.any(speeds <= 0.0) or np.any(noise <= 0.0):
        raise FloatingPointError("speed/noise QC must be positive")

    calibration = pd.read_csv(CALIBRATION_SUMMARY)
    d5 = calibration[calibration["cantilever"] == "D5"]
    if d5.shape[0] != 1:
        raise RuntimeError("D5 calibration row is missing or duplicated")

    return {
        "distance_nm": distance_nm,
        "force_pN": force,
        "sources": sources,
        "times_min": times,
        "inventory": inventory,
        "qc": qc,
        "rows": rows,
        "columns": columns,
        "speeds_um_per_s": speeds,
        "noise_pN": noise,
        "spring_N_per_m": float(d5.iloc[0]["spring_constant_N_per_m"]),
        "spring_sd_N_per_m": float(
            d5.iloc[0]["spring_constant_repeatability_sd_N_per_m"]
        ),
    }


def reconstruct_contact_heights() -> np.ndarray:
    paths = sorted((ROOT / "31-08-26").glob("*.jpk-force-map"))
    if len(paths) != 5:
        raise RuntimeError("contact-height reconstruction requires five raw maps")
    sources = [base.load_source(path.resolve(), 0) for path in paths]
    sources.sort(key=lambda source: source.timestamp)
    raw.contact_analysis(sources)
    heights = np.empty((4, 64), dtype=np.float64)
    for map_index, source in enumerate(sources[1:]):
        curves = sorted(source.curves, key=lambda curve: int(curve.point_index))
        if len(curves) != 64:
            raise RuntimeError("contact-height reconstruction lost curves")
        for curve in curves:
            if curve.contact_fit is None or curve.point_index is None:
                raise RuntimeError("contact fit was not attached to a raw curve")
            heights[map_index, int(curve.point_index)] = (
                -curve.contact_fit.intercept_V / curve.contact_fit.slope_V_per_m * 1e9
            )
    if not np.all(np.isfinite(heights)):
        raise FloatingPointError("contact heights contain non-finite values")
    return heights


def model_free_analysis(data: dict, rng: np.random.Generator) -> dict:
    distance = data["distance_nm"]
    force = data["force_pN"]
    noise = data["noise_pN"]
    rows = data["rows"]
    times = data["times_min"]
    near_index = int(np.flatnonzero(np.isclose(distance, 20.0))[0])
    far_index = int(np.flatnonzero(np.isclose(distance, 200.0))[0])
    denominator = force[:, :, near_index] - force[:, :, far_index]
    threshold = np.maximum(3.0 * noise, NORMALIZATION_MIN_DENOMINATOR_PN)
    valid_normalization = denominator > threshold
    normalized = (
        force - force[:, :, far_index][:, :, None]
    ) / denominator[:, :, None]
    normalized[~valid_normalization] = np.nan

    curve_rows: list[dict] = []
    summary_rows: list[dict] = []
    for map_index, map_order in enumerate(range(2, 6)):
        valid = valid_normalization[map_index]
        for point in range(64):
            for distance_index, distance_nm in enumerate(distance):
                curve_rows.append(
                    {
                        "map_order": map_order,
                        "relative_time_min": float(times[map_index]),
                        "point_index": point,
                        "row": int(data["rows"][point]),
                        "column": int(data["columns"][point]),
                        "distance_nm": float(distance_nm),
                        "normalization_valid": bool(valid[point]),
                        "normalization_denominator_pN": float(denominator[map_index, point]),
                        "normalized_force": float(normalized[map_index, point, distance_index]),
                    }
                )
        for distance_index, distance_nm in enumerate(distance):
            values = normalized[map_index, :, distance_index]
            selected = values[np.isfinite(values)]
            summary_rows.append(
                {
                    "map_order": map_order,
                    "relative_time_min": float(times[map_index]),
                    "distance_nm": float(distance_nm),
                    "valid_pixels": int(selected.size),
                    "normalized_force_median": float(np.median(selected)) if selected.size else float("nan"),
                    "normalized_force_q25": finite_quantile(selected, 0.25),
                    "normalized_force_q75": finite_quantile(selected, 0.75),
                }
            )

    local_rows: list[dict] = []
    local_values: dict[str, np.ndarray] = {}
    far_mask = (distance >= 150.0) & (distance <= 200.0)
    windows = {"20_60nm": (20.0, 60.0), "30_80nm": (30.0, 80.0), "40_100nm": (40.0, 100.0)}
    for label, (low, high) in windows.items():
        output = np.full((4, 64), np.nan, dtype=np.float64)
        for map_index, map_order in enumerate(range(2, 6)):
            for point in range(64):
                baseline = float(np.median(force[map_index, point, far_mask]))
                signal = force[map_index, point] - baseline
                selected = (
                    (distance >= low)
                    & (distance <= high)
                    & (signal > max(2.0 * noise[map_index, point], 5.0))
                )
                fit_status = "insufficient_positive_signal"
                slope = intercept = r2 = apparent_lambda = float("nan")
                if np.count_nonzero(selected) >= 5:
                    fit = stats.linregress(distance[selected], np.log(signal[selected]))
                    slope = float(fit.slope)
                    intercept = float(fit.intercept)
                    r2 = float(fit.rvalue**2)
                    if slope < 0.0:
                        apparent_lambda = -1.0 / slope
                        fit_status = "negative_log_slope"
                        output[map_index, point] = apparent_lambda
                    else:
                        fit_status = "nondecaying_log_slope"
                local_rows.append(
                    {
                        "map_order": map_order,
                        "relative_time_min": float(times[map_index]),
                        "point_index": point,
                        "row": int(data["rows"][point]),
                        "column": int(data["columns"][point]),
                        "window": label,
                        "far_baseline_150_200nm_pN": baseline,
                        "positive_fit_bins": int(np.count_nonzero(selected)),
                        "log_force_slope_per_nm": slope,
                        "log_force_intercept": intercept,
                        "log_linear_r2": r2,
                        "apparent_local_decay_length_nm": apparent_lambda,
                        "fit_status": fit_status,
                    }
                )
        local_values[label] = output

    paired_rows: list[dict] = []
    for target_nm in (30.0, 50.0, 80.0, 100.0):
        index = int(np.flatnonzero(np.isclose(distance, target_nm))[0])
        before = normalized[0, :, index]
        after = normalized[3, :, index]
        valid = np.isfinite(before) & np.isfinite(after)
        difference = after[valid] - before[valid]
        valid_rows = rows[valid]
        low, high = row_cluster_bootstrap_ci(difference, valid_rows, rng)
        wilcoxon = stats.wilcoxon(difference) if difference.size and np.any(difference != 0.0) else None
        paired_rows.append(
            {
                "metric": "normalized_force",
                "distance_nm": target_nm,
                "paired_pixels": int(difference.size),
                "map2_median": float(np.median(before[valid])),
                "map5_median": float(np.median(after[valid])),
                "map5_minus_map2_median": float(np.median(difference)),
                "difference_q25": finite_quantile(difference, 0.25),
                "difference_q75": finite_quantile(difference, 0.75),
                "wilcoxon_p_naive_pixels": float(wilcoxon.pvalue) if wilcoxon else float("nan"),
                "exact_physical_row_signflip_p": exact_row_signflip_p(difference, valid_rows),
                "row_cluster_bootstrap_median_95CI_low": low,
                "row_cluster_bootstrap_median_95CI_high": high,
            }
        )
    for label, values in local_values.items():
        valid = np.isfinite(values[0]) & np.isfinite(values[3])
        difference = values[3, valid] - values[0, valid]
        valid_rows = rows[valid]
        low, high = row_cluster_bootstrap_ci(difference, valid_rows, rng)
        test = stats.wilcoxon(difference) if difference.size and np.any(difference != 0.0) else None
        paired_rows.append(
            {
                "metric": f"local_decay_{label}",
                "distance_nm": float("nan"),
                "paired_pixels": int(difference.size),
                "map2_median": finite_quantile(values[0], 0.5),
                "map5_median": finite_quantile(values[3], 0.5),
                "map5_minus_map2_median": finite_quantile(difference, 0.5),
                "difference_q25": finite_quantile(difference, 0.25),
                "difference_q75": finite_quantile(difference, 0.75),
                "wilcoxon_p_naive_pixels": float(test.pvalue) if test else float("nan"),
                "exact_physical_row_signflip_p": exact_row_signflip_p(difference, valid_rows),
                "row_cluster_bootstrap_median_95CI_low": low,
                "row_cluster_bootstrap_median_95CI_high": high,
            }
        )

    paired_q = benjamini_hochberg(
        [float(row["exact_physical_row_signflip_p"]) for row in paired_rows]
    )
    for row, q_value in zip(paired_rows, paired_q, strict=True):
        row["exact_row_signflip_BH_q_across_shape_tests"] = q_value

    write_csv(RESULTS / "model_free_normalized_curves.csv", curve_rows)
    write_csv(RESULTS / "model_free_normalized_summary.csv", summary_rows)
    write_csv(RESULTS / "model_free_local_decay_fits.csv", local_rows)
    write_csv(RESULTS / "model_free_paired_shape_tests.csv", paired_rows)

    figure, axes = plt.subplots(1, 2, figsize=(14.0, 5.4))
    colors = ("#264653", "#2a9d8f", "#e9c46a", "#e76f51")
    for map_index, (map_order, color) in enumerate(zip(range(2, 6), colors, strict=True)):
        values = normalized[map_index]
        median = np.nanmedian(values, axis=0)
        q25 = np.nanquantile(values, 0.25, axis=0)
        q75 = np.nanquantile(values, 0.75, axis=0)
        axes[0].fill_between(distance, q25, q75, color=color, alpha=0.12)
        axes[0].plot(distance, median, color=color, lw=2.1, label=f"map{map_order}, {times[map_index]:.1f} min")
    axes[0].axhline(0.0, color="0.5", lw=0.8)
    axes[0].set(xlabel="Nominal separation D (nm)", ylabel="Endpoint-normalized force", title="Model-free curve shape: median and IQR")
    axes[0].grid(alpha=0.2)
    axes[0].legend(frameon=False)
    label = "30_80nm"
    for point in range(64):
        axes[1].plot(times, local_values[label][:, point], color="0.65", alpha=0.18, lw=0.8)
    medians = np.asarray([finite_quantile(values, 0.5) for values in local_values[label]])
    q25 = np.asarray([finite_quantile(values, 0.25) for values in local_values[label]])
    q75 = np.asarray([finite_quantile(values, 0.75) for values in local_values[label]])
    axes[1].errorbar(times, medians, yerr=np.vstack([medians - q25, q75 - medians]), fmt="o-", color="#d1495b", lw=2.2, capsize=4)
    axes[1].set(xlabel="Time from map2 start (min)", ylabel="Apparent local decay length (nm)", title="30–80 nm log-slope diagnostic")
    axes[1].grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(FIGURES / "model_free_curve_shape_and_local_decay.png", dpi=230)
    plt.close(figure)

    return {
        "normalized": normalized,
        "normalization_valid": valid_normalization,
        "local_decay": local_values,
        "paired_rows": paired_rows,
    }


def fit_model_variants(data: dict) -> dict:
    distance = data["distance_nm"]
    observed = data["force_pN"].reshape(256, 37)
    speed = data["speeds_um_per_s"].reshape(256)
    nonlinear_library = pb_sphere_plane_library(distance)
    sphere_plane_hyd = hydrodynamic_force_pN(distance, speed, RADIUS_M)
    equal_sphere_hyd = hydrodynamic_force_pN(
        distance, speed, EFFECTIVE_EQUAL_SPHERE_RADIUS_M
    )
    libraries = {
        PRIMARY_MODEL: (nonlinear_library, sphere_plane_hyd),
        "sphere_plane_nonlinear_PB_no_hydrodynamics": (
            nonlinear_library,
            np.zeros_like(sphere_plane_hyd),
        ),
        "sphere_plane_linear_HHF_fixed_no_slip_hydrodynamics": (
            linear_hhf_library(distance, RADIUS_M),
            sphere_plane_hyd,
        ),
        "equal_sphere_linear_HHF_fixed_no_slip_hydrodynamics": (
            linear_hhf_library(distance, EFFECTIVE_EQUAL_SPHERE_RADIUS_M),
            equal_sphere_hyd,
        ),
    }
    fits: dict[str, dict[str, np.ndarray]] = {}
    for name, (library, hydrodynamic) in libraries.items():
        fits[name] = batch_grid_fit(
            observed,
            library,
            hydrodynamic,
            keep_scores=name == PRIMARY_MODEL,
        )

    all_aic = np.column_stack([fits[name]["aic"] for name in libraries])
    best_aic = np.min(all_aic, axis=1)
    rows: list[dict] = []
    summary_rows: list[dict] = []
    for model_index, name in enumerate(libraries):
        fit = fits[name]
        for curve_index in range(256):
            map_index = curve_index // 64
            point = curve_index % 64
            rows.append(
                {
                    "model": name,
                    "map_order": map_index + 2,
                    "relative_time_min": float(data["times_min"][map_index]),
                    "point_index": point,
                    "row": int(data["rows"][point]),
                    "column": int(data["columns"][point]),
                    "lambda_D_nm": float(fit["lambda_D_nm"][curve_index]),
                    "zeta_magnitude_mV": float(fit["zeta_magnitude_mV"][curve_index]),
                    "zeta_signed_silica_mV": -float(fit["zeta_magnitude_mV"][curve_index]),
                    "baseline_pN": float(fit["baseline_pN"][curve_index]),
                    "rss_pN2": float(fit["rss_pN2"][curve_index]),
                    "rmse_pN": float(fit["rmse_pN"][curve_index]),
                    "r2": float(fit["r2"][curve_index]),
                    "aic_correlated_bins_diagnostic": float(fit["aic"][curve_index]),
                    "bic_correlated_bins_diagnostic": float(fit["bic"][curve_index]),
                    "delta_aic_to_best_same_curve": float(
                        fit["aic"][curve_index] - best_aic[curve_index]
                    ),
                    "aic_winner_tie_inclusive": bool(
                        np.isclose(fit["aic"][curve_index], best_aic[curve_index], atol=1e-9)
                    ),
                }
            )
        for map_index, map_order in enumerate(range(2, 6)):
            selected = slice(map_index * 64, (map_index + 1) * 64)
            summary_rows.append(
                {
                    "model": name,
                    "map_order": map_order,
                    "relative_time_min": float(data["times_min"][map_index]),
                    "curves": 64,
                    "lambda_D_nm_median": finite_quantile(fit["lambda_D_nm"][selected], 0.5),
                    "lambda_D_nm_q25": finite_quantile(fit["lambda_D_nm"][selected], 0.25),
                    "lambda_D_nm_q75": finite_quantile(fit["lambda_D_nm"][selected], 0.75),
                    "zeta_magnitude_mV_median": finite_quantile(fit["zeta_magnitude_mV"][selected], 0.5),
                    "zeta_magnitude_mV_q25": finite_quantile(fit["zeta_magnitude_mV"][selected], 0.25),
                    "zeta_magnitude_mV_q75": finite_quantile(fit["zeta_magnitude_mV"][selected], 0.75),
                    "rmse_pN_median": finite_quantile(fit["rmse_pN"][selected], 0.5),
                    "r2_median": finite_quantile(fit["r2"][selected], 0.5),
                    "aic_winner_fraction": float(
                        np.mean(np.isclose(fit["aic"][selected], best_aic[selected], atol=1e-9))
                    ),
                }
            )
    write_csv(RESULTS / "surface_model_comparison_per_curve.csv", rows)
    write_csv(RESULTS / "surface_model_comparison_summary.csv", summary_rows)

    primary = fits[PRIMARY_MODEL]
    prediction_rows: list[dict] = []
    for curve_index in range(256):
        map_index = curve_index // 64
        point = curve_index % 64
        for distance_index, distance_nm in enumerate(distance):
            prediction_rows.append(
                {
                    "map_order": map_index + 2,
                    "relative_time_min": float(data["times_min"][map_index]),
                    "point_index": point,
                    "row": int(data["rows"][point]),
                    "column": int(data["columns"][point]),
                    "distance_nm": float(distance_nm),
                    "observed_force_pN": float(observed[curve_index, distance_index]),
                    "fitted_total_force_pN": float(primary["prediction_pN"][curve_index, distance_index]),
                    "residual_observed_minus_fitted_pN": float(primary["residual_pN"][curve_index, distance_index]),
                    "fixed_sphere_plane_hydrodynamic_force_pN": float(
                        sphere_plane_hyd[curve_index, distance_index]
                    ),
                    "fitted_baseline_pN": float(primary["baseline_pN"][curve_index]),
                    "lambda_D_nm": float(primary["lambda_D_nm"][curve_index]),
                    "zeta_magnitude_mV": float(primary["zeta_magnitude_mV"][curve_index]),
                }
            )
    write_csv(RESULTS / "sphere_plane_nonlinear_predictions.csv", prediction_rows)

    figure, axes = plt.subplots(1, 3, figsize=(16.0, 5.0))
    model_colors = {
        PRIMARY_MODEL: "#2a9d8f",
        "sphere_plane_nonlinear_PB_no_hydrodynamics": "#457b9d",
        "sphere_plane_linear_HHF_fixed_no_slip_hydrodynamics": "#e9c46a",
        "equal_sphere_linear_HHF_fixed_no_slip_hydrodynamics": "#e76f51",
    }
    for name in libraries:
        summary = [row for row in summary_rows if row["model"] == name]
        label = name.replace("_fixed_no_slip_hydrodynamics", "+hyd").replace("_no_hydrodynamics", "−hyd")
        axes[0].plot(data["times_min"], [row["lambda_D_nm_median"] for row in summary], "o-", color=model_colors[name], label=label)
        axes[1].plot(data["times_min"], [row["zeta_magnitude_mV_median"] for row in summary], "o-", color=model_colors[name])
        axes[2].plot(data["times_min"], [row["rmse_pN_median"] for row in summary], "o-", color=model_colors[name])
    axes[0].set_ylabel("Median apparent Debye length (nm)")
    axes[1].set_ylabel("Median apparent |psi| (mV)")
    axes[2].set_ylabel("Median per-curve RMSE (pN)")
    for axis in axes:
        axis.set_xlabel("Time from map2 start (min)")
        axis.grid(alpha=0.2)
    axes[0].legend(frameon=False, fontsize=8)
    figure.suptitle("Uniform-OLS comparison of geometry, PB linearization, and hydrodynamics")
    figure.tight_layout()
    figure.savefig(FIGURES / "surface_model_comparison.png", dpi=230)
    plt.close(figure)

    return {
        "fits": fits,
        "libraries": libraries,
        "nonlinear_library": nonlinear_library,
        "sphere_plane_hyd_pN": sphere_plane_hyd,
        "summary_rows": summary_rows,
    }


def contact_zero_analysis(data: dict, nonlinear_library_at_zero: np.ndarray) -> dict:
    observed = data["force_pN"].reshape(256, 37)
    speed = data["speeds_um_per_s"].reshape(256)
    fit_rows: list[dict] = []
    summary_rows: list[dict] = []
    offset_results: dict[float, dict[str, np.ndarray]] = {}
    for offset_nm in CONTACT_OFFSETS_NM:
        true_distance = data["distance_nm"] + offset_nm
        if np.any(true_distance <= 0.0):
            raise RuntimeError("contact offset generated non-positive separation")
        library = (
            nonlinear_library_at_zero
            if np.isclose(offset_nm, 0.0)
            else pb_sphere_plane_library(true_distance)
        )
        hydrodynamic = hydrodynamic_force_pN(true_distance, speed, RADIUS_M)
        fit = batch_grid_fit(observed, library, hydrodynamic, keep_scores=False)
        offset_results[float(offset_nm)] = fit
        for curve_index in range(256):
            map_index = curve_index // 64
            point = curve_index % 64
            fit_rows.append(
                {
                    "fixed_contact_offset_nm": float(offset_nm),
                    "distance_interpretation": "D_true=D_nominal+fixed_contact_offset",
                    "map_order": map_index + 2,
                    "relative_time_min": float(data["times_min"][map_index]),
                    "point_index": point,
                    "row": int(data["rows"][point]),
                    "column": int(data["columns"][point]),
                    "lambda_D_nm": float(fit["lambda_D_nm"][curve_index]),
                    "zeta_magnitude_mV": float(fit["zeta_magnitude_mV"][curve_index]),
                    "baseline_pN": float(fit["baseline_pN"][curve_index]),
                    "rmse_pN": float(fit["rmse_pN"][curve_index]),
                    "r2": float(fit["r2"][curve_index]),
                    "rss_pN2": float(fit["rss_pN2"][curve_index]),
                }
            )
        for map_index, map_order in enumerate(range(2, 6)):
            selected = slice(map_index * 64, (map_index + 1) * 64)
            summary_rows.append(
                {
                    "fixed_contact_offset_nm": float(offset_nm),
                    "map_order": map_order,
                    "relative_time_min": float(data["times_min"][map_index]),
                    "lambda_D_nm_median": finite_quantile(fit["lambda_D_nm"][selected], 0.5),
                    "lambda_D_nm_q25": finite_quantile(fit["lambda_D_nm"][selected], 0.25),
                    "lambda_D_nm_q75": finite_quantile(fit["lambda_D_nm"][selected], 0.75),
                    "zeta_magnitude_mV_median": finite_quantile(fit["zeta_magnitude_mV"][selected], 0.5),
                    "zeta_magnitude_mV_q25": finite_quantile(fit["zeta_magnitude_mV"][selected], 0.25),
                    "zeta_magnitude_mV_q75": finite_quantile(fit["zeta_magnitude_mV"][selected], 0.75),
                    "rmse_pN_median": finite_quantile(fit["rmse_pN"][selected], 0.5),
                }
            )

    rss_by_offset = np.vstack(
        [offset_results[float(offset)]["rss_pN2"] for offset in CONTACT_OFFSETS_NM]
    )
    best_offset_index = np.argmin(rss_by_offset, axis=0)
    nuisance_rows: list[dict] = []
    for curve_index in range(256):
        map_index = curve_index // 64
        point = curve_index % 64
        offset_index = int(best_offset_index[curve_index])
        offset = float(CONTACT_OFFSETS_NM[offset_index])
        fit = offset_results[offset]
        zero_fit = offset_results[0.0]
        nuisance_rows.append(
            {
                "map_order": map_index + 2,
                "relative_time_min": float(data["times_min"][map_index]),
                "point_index": point,
                "row": int(data["rows"][point]),
                "column": int(data["columns"][point]),
                "best_discrete_contact_offset_nm": offset,
                "best_offset_at_search_boundary": bool(
                    offset_index in (0, CONTACT_OFFSETS_NM.size - 1)
                ),
                "lambda_D_nm_with_profiled_contact_offset": float(
                    fit["lambda_D_nm"][curve_index]
                ),
                "zeta_magnitude_mV_with_profiled_contact_offset": float(
                    fit["zeta_magnitude_mV"][curve_index]
                ),
                "baseline_pN_with_profiled_contact_offset": float(
                    fit["baseline_pN"][curve_index]
                ),
                "rmse_pN_with_profiled_contact_offset": float(fit["rmse_pN"][curve_index]),
                "rss_improvement_vs_zero_offset_pN2": float(
                    zero_fit["rss_pN2"][curve_index] - fit["rss_pN2"][curve_index]
                ),
                "zero_offset_lambda_D_nm": float(zero_fit["lambda_D_nm"][curve_index]),
                "zero_offset_zeta_magnitude_mV": float(
                    zero_fit["zeta_magnitude_mV"][curve_index]
                ),
            }
        )
    write_csv(RESULTS / "contact_zero_fixed_offset_fits.csv", fit_rows)
    write_csv(RESULTS / "contact_zero_fixed_offset_summary.csv", summary_rows)
    write_csv(RESULTS / "contact_zero_profiled_nuisance.csv", nuisance_rows)

    figure, axes = plt.subplots(1, 3, figsize=(16.0, 5.0))
    colors = ("#264653", "#2a9d8f", "#e9c46a", "#e76f51")
    for map_index, (map_order, color) in enumerate(zip(range(2, 6), colors, strict=True)):
        selected = [row for row in summary_rows if row["map_order"] == map_order]
        axes[0].plot(
            CONTACT_OFFSETS_NM,
            [row["lambda_D_nm_median"] for row in selected],
            "o-",
            color=color,
            label=f"map{map_order}",
        )
        axes[1].plot(
            CONTACT_OFFSETS_NM,
            [row["zeta_magnitude_mV_median"] for row in selected],
            "o-",
            color=color,
        )
    axes[0].set_ylabel("Median apparent Debye length (nm)")
    axes[1].set_ylabel("Median apparent |psi| (mV)")
    axes[0].legend(frameon=False)
    profiled = pd.DataFrame(nuisance_rows)
    axes[2].boxplot(
        [
            profiled[profiled["map_order"] == map_order][
                "best_discrete_contact_offset_nm"
            ].to_numpy()
            for map_order in range(2, 6)
        ],
        tick_labels=[f"map{order}" for order in range(2, 6)],
        showfliers=False,
    )
    axes[2].set_ylabel("Best discrete contact offset (nm)")
    for axis in axes[:2]:
        axis.axvline(0.0, color="0.5", lw=0.8)
        axis.set_xlabel("Fixed D offset (nm); Dtrue = Dnominal + offset")
        axis.grid(alpha=0.2)
    axes[2].grid(axis="y", alpha=0.2)
    figure.suptitle("Contact/separation-zero sensitivity of sphere–plane nonlinear-PB fits")
    figure.tight_layout()
    figure.savefig(FIGURES / "contact_zero_sensitivity.png", dpi=230)
    plt.close(figure)
    return {
        "offset_results": offset_results,
        "summary_rows": summary_rows,
        "nuisance_rows": nuisance_rows,
    }


def identifiability_analysis(data: dict, model_results: dict) -> dict:
    primary = model_results["fits"][PRIMARY_MODEL]
    scores = primary["scores_pN2"]
    best_rss = primary["rss_pN2"]
    noise = data["noise_pN"].reshape(256)
    dof = data["distance_nm"].size - 3
    interval_rows: list[dict] = []
    for curve_index in range(256):
        scale2 = max(noise[curve_index] ** 2, best_rss[curve_index] / dof)
        threshold = best_rss[curve_index] + PROFILE_DELTA_TWO_PARAMETERS * scale2
        selected = scores[:, curve_index] <= threshold
        profile_lambda = FLAT_LAMBDA_NM[selected]
        profile_zeta = FLAT_ZETA_MV[selected]
        correlation = (
            float(np.corrcoef(np.log(profile_lambda), profile_zeta)[0, 1])
            if profile_lambda.size >= 3
            and np.ptp(profile_lambda) > 0.0
            and np.ptp(profile_zeta) > 0.0
            else float("nan")
        )
        map_index = curve_index // 64
        point = curve_index % 64
        interval_rows.append(
            {
                "map_order": map_index + 2,
                "relative_time_min": float(data["times_min"][map_index]),
                "point_index": point,
                "row": int(data["rows"][point]),
                "column": int(data["columns"][point]),
                "best_lambda_D_nm": float(primary["lambda_D_nm"][curve_index]),
                "lambda_profile_low_nm": float(np.min(profile_lambda)),
                "lambda_profile_high_nm": float(np.max(profile_lambda)),
                "best_zeta_magnitude_mV": float(primary["zeta_magnitude_mV"][curve_index]),
                "zeta_profile_low_mV": float(np.min(profile_zeta)),
                "zeta_profile_high_mV": float(np.max(profile_zeta)),
                "profile_grid_points": int(np.count_nonzero(selected)),
                "profile_log_lambda_zeta_correlation": correlation,
                "profile_lambda_hits_grid_boundary": bool(
                    np.isclose(np.min(profile_lambda), LAMBDA_GRID_NM[0])
                    or np.isclose(np.max(profile_lambda), LAMBDA_GRID_NM[-1])
                ),
                "profile_zeta_hits_grid_boundary": bool(
                    np.isclose(np.min(profile_zeta), ZETA_GRID_MV[0])
                    or np.isclose(np.max(profile_zeta), ZETA_GRID_MV[-1])
                ),
                "profile_status": "diagnostic_correlated_distance_bins_not_formal_CI",
            }
        )

    library = model_results["nonlinear_library"]
    hyd = model_results["sphere_plane_hyd_pN"].reshape(4, 64, 37)
    force = data["force_pN"]
    profile_rows: list[dict] = []
    map_profile_summary: list[dict] = []
    contour_payload: list[tuple[np.ndarray, float, float]] = []
    for map_index, map_order in enumerate(range(2, 6)):
        observed = np.median(force[map_index], axis=0)
        fixed_hyd = np.median(hyd[map_index], axis=0)
        adjusted = observed - fixed_hyd
        q25 = np.quantile(force[map_index], 0.25, axis=0)
        q75 = np.quantile(force[map_index], 0.75, axis=0)
        sigma = np.maximum((q75 - q25) / 1.349, 2.0)
        weights = 1.0 / sigma**2
        baseline = np.sum((adjusted[None, :] - library) * weights[None, :], axis=1)
        baseline /= np.sum(weights)
        residual = library + baseline[:, None] - adjusted[None, :]
        objective = np.einsum("ij,j,ij->i", residual, weights, residual)
        best_index = int(np.argmin(objective))
        reduced = float(objective[best_index] / (distance_size := (data["distance_nm"].size - 3)))
        scale = max(reduced, 1.0)
        delta = objective - objective[best_index]
        selected = delta <= PROFILE_DELTA_TWO_PARAMETERS * scale
        matrix = delta.reshape(ZETA_GRID_MV.size, LAMBDA_GRID_NM.size)
        contour_payload.append((matrix, scale, float(objective[best_index])))
        for lambda_index, lambda_nm in enumerate(LAMBDA_GRID_NM):
            profile_rows.append(
                {
                    "map_order": map_order,
                    "parameter": "lambda_D_nm",
                    "parameter_value": float(lambda_nm),
                    "profile_delta_objective": float(np.min(matrix[:, lambda_index])),
                    "diagnostic_scale": scale,
                }
            )
        for zeta_index, zeta_mV in enumerate(ZETA_GRID_MV):
            profile_rows.append(
                {
                    "map_order": map_order,
                    "parameter": "zeta_magnitude_mV",
                    "parameter_value": float(zeta_mV),
                    "profile_delta_objective": float(np.min(matrix[zeta_index, :])),
                    "diagnostic_scale": scale,
                }
            )
        map_profile_summary.append(
            {
                "map_order": map_order,
                "relative_time_min": float(data["times_min"][map_index]),
                "best_lambda_D_nm": float(FLAT_LAMBDA_NM[best_index]),
                "lambda_profile_low_nm": float(np.min(FLAT_LAMBDA_NM[selected])),
                "lambda_profile_high_nm": float(np.max(FLAT_LAMBDA_NM[selected])),
                "best_zeta_magnitude_mV": float(FLAT_ZETA_MV[best_index]),
                "zeta_profile_low_mV": float(np.min(FLAT_ZETA_MV[selected])),
                "zeta_profile_high_mV": float(np.max(FLAT_ZETA_MV[selected])),
                "best_baseline_pN": float(baseline[best_index]),
                "weighted_objective": float(objective[best_index]),
                "reduced_objective_using_spatial_IQR_scale": reduced,
                "profile_status": "map_median_spatial_IQR_weighted_diagnostic_not_formal_CI",
            }
        )
    write_csv(RESULTS / "per_curve_parameter_profile_intervals.csv", interval_rows)
    write_csv(RESULTS / "map_median_parameter_profiles.csv", profile_rows)
    write_csv(RESULTS / "map_median_parameter_profile_summary.csv", map_profile_summary)

    figure, axes = plt.subplots(1, 4, figsize=(18.0, 4.2), sharex=True, sharey=True)
    for axis, map_order, payload, summary in zip(
        axes, range(2, 6), contour_payload, map_profile_summary, strict=True
    ):
        matrix, scale, _ = payload
        levels = np.asarray([2.30, 6.18, 11.83]) * scale
        axis.contour(
            LAMBDA_GRID_NM,
            ZETA_GRID_MV,
            matrix,
            levels=levels,
            colors=("#2a9d8f", "#e9c46a", "#e76f51"),
        )
        axis.scatter(
            summary["best_lambda_D_nm"],
            summary["best_zeta_magnitude_mV"],
            marker="x",
            color="black",
            s=45,
        )
        axis.set_xscale("log")
        axis.set_title(f"map{map_order}")
        axis.set_xlabel("Apparent Debye length (nm)")
        axis.grid(alpha=0.15)
    axes[0].set_ylabel("Apparent |psi| (mV)")
    figure.suptitle("Map-median lambda–potential profile contours (diagnostic spatial-IQR scale)")
    figure.tight_layout()
    figure.savefig(FIGURES / "parameter_identifiability_profiles.png", dpi=230)
    plt.close(figure)
    return {
        "interval_rows": interval_rows,
        "map_profile_summary": map_profile_summary,
    }


def paired_parameter_and_spatial_analysis(
    data: dict,
    model_free: dict,
    model_results: dict,
    contact_heights_nm: np.ndarray,
    contact_results: dict,
    rng: np.random.Generator,
) -> dict:
    primary = model_results["fits"][PRIMARY_MODEL]
    parameter_arrays = {
        "apparent_lambda_D_nm": primary["lambda_D_nm"].reshape(4, 64),
        "apparent_zeta_magnitude_mV": primary["zeta_magnitude_mV"].reshape(4, 64),
        "fitted_baseline_pN": primary["baseline_pN"].reshape(4, 64),
        "fit_RMSE_pN": primary["rmse_pN"].reshape(4, 64),
        "contact_height_nm": contact_heights_nm,
        "local_contact_InvOLS_nm_per_V": data["qc"][
            "local_contact_InvOLS_nm_per_V"
        ].to_numpy(dtype=np.float64).reshape(4, 64),
        "far_slope_pN_per_100nm": data["qc"]["far_slope_pN_per_100nm"].to_numpy(dtype=np.float64).reshape(4, 64),
    }
    for target in (20.0, 50.0, 100.0, 200.0):
        index = int(np.flatnonzero(np.isclose(data["distance_nm"], target))[0])
        parameter_arrays[f"force_{int(target)}nm_pN"] = data["force_pN"][:, :, index]
    parameter_arrays["model_free_local_decay_30_80nm"] = model_free["local_decay"]["30_80nm"]
    normalized_50_index = int(
        np.flatnonzero(np.isclose(data["distance_nm"], 50.0))[0]
    )
    parameter_arrays["model_free_normalized_force_50nm"] = model_free[
        "normalized"
    ][:, :, normalized_50_index]
    profiled = pd.DataFrame(contact_results["nuisance_rows"])
    parameter_arrays["best_discrete_contact_offset_nm"] = profiled[
        "best_discrete_contact_offset_nm"
    ].to_numpy(dtype=np.float64).reshape(4, 64)

    times = data["times_min"]
    rows = data["rows"]
    columns = data["columns"]
    slope_rows: list[dict] = []
    paired_rows: list[dict] = []
    for metric, values in parameter_arrays.items():
        for point in range(64):
            y = values[:, point]
            valid = np.isfinite(y)
            slope = intercept = r2 = fitted_change = float("nan")
            if np.count_nonzero(valid) >= 3 and np.ptp(times[valid]) > 0.0:
                fit = stats.linregress(times[valid], y[valid])
                slope = float(fit.slope)
                intercept = float(fit.intercept)
                r2 = float(fit.rvalue**2)
                fitted_change = slope * float(times[-1] - times[0])
            slope_rows.append(
                {
                    "metric": metric,
                    "point_index": point,
                    "row": int(rows[point]),
                    "column": int(columns[point]),
                    "valid_timepoints": int(np.count_nonzero(valid)),
                    "linear_slope_per_min": slope,
                    "linear_intercept": intercept,
                    "linear_r2": r2,
                    "fitted_map5_minus_map2_change": fitted_change,
                    "observed_map5_minus_map2": float(y[3] - y[0])
                    if np.isfinite(y[3]) and np.isfinite(y[0])
                    else float("nan"),
                    "strictly_decreasing": bool(np.all(np.diff(y) < 0.0))
                    if np.all(np.isfinite(y))
                    else False,
                    "strictly_increasing": bool(np.all(np.diff(y) > 0.0))
                    if np.all(np.isfinite(y))
                    else False,
                }
            )
        valid = np.isfinite(values[0]) & np.isfinite(values[3])
        difference = values[3, valid] - values[0, valid]
        valid_rows = rows[valid]
        low, high = row_cluster_bootstrap_ci(difference, valid_rows, rng)
        test = stats.wilcoxon(difference) if difference.size and np.any(difference != 0.0) else None
        rho, rho_p, count = stable_spearman(values[0], values[3])
        paired_rows.append(
            {
                "metric": metric,
                "paired_pixels": int(difference.size),
                "map2_median": finite_quantile(values[0], 0.5),
                "map5_median": finite_quantile(values[3], 0.5),
                "map5_minus_map2_median": finite_quantile(difference, 0.5),
                "difference_q25": finite_quantile(difference, 0.25),
                "difference_q75": finite_quantile(difference, 0.75),
                "negative_difference_fraction": float(np.mean(difference < 0.0)),
                "strict_decrease_fraction": float(np.mean(np.all(np.diff(values[:, valid], axis=0) < 0.0, axis=0))),
                "wilcoxon_p_naive_pixels": float(test.pvalue) if test else float("nan"),
                "exact_physical_row_signflip_p": exact_row_signflip_p(difference, valid_rows),
                "row_cluster_bootstrap_median_95CI_low": low,
                "row_cluster_bootstrap_median_95CI_high": high,
                "map2_map5_spatial_spearman_rho": rho,
                "map2_map5_spatial_spearman_p_naive": rho_p,
                "map2_map5_spatial_spearman_pixels": count,
            }
        )
    paired_q = benjamini_hochberg(
        [float(row["exact_physical_row_signflip_p"]) for row in paired_rows]
    )
    for row, q_value in zip(paired_rows, paired_q, strict=True):
        row["exact_row_signflip_BH_q_across_parameter_tests"] = q_value
    write_csv(RESULTS / "same_pixel_parameter_time_slopes.csv", slope_rows)
    write_csv(RESULTS / "same_pixel_paired_parameter_tests.csv", paired_rows)

    spatial_rows: list[dict] = []
    variogram_rows: list[dict] = []
    for metric, values in parameter_arrays.items():
        datasets = [(f"map{order}", values[index]) for index, order in enumerate(range(2, 6))]
        datasets.append(("map5_minus_map2", values[3] - values[0]))
        for dataset_label, vector in datasets:
            if np.count_nonzero(np.isfinite(vector)) < 16:
                continue
            matrix = physical_matrix(vector, rows, columns)
            flattened = matrix.reshape(-1)
            moran, moran_p = moran_permutation_test(flattened, rng)
            h_x = np.asarray([flattened[i] for i, _ in HORIZONTAL_PAIRS])
            h_y = np.asarray([flattened[j] for _, j in HORIZONTAL_PAIRS])
            v_x = np.asarray([flattened[i] for i, _ in VERTICAL_PAIRS])
            v_y = np.asarray([flattened[j] for _, j in VERTICAL_PAIRS])
            hrho, hp, _ = stable_spearman(h_x, h_y)
            vrho, vp, _ = stable_spearman(v_x, v_y)
            spatial_rows.append(
                {
                    "metric": metric,
                    "dataset": dataset_label,
                    "Moran_I_rook_neighbors": moran,
                    "Moran_permutation_p_two_sided": moran_p,
                    "horizontal_neighbor_Spearman_rho": hrho,
                    "horizontal_neighbor_p_naive": hp,
                    "vertical_neighbor_Spearman_rho": vrho,
                    "vertical_neighbor_p_naive": vp,
                    "finite_pixels": int(np.count_nonzero(np.isfinite(flattened))),
                    "variance": float(np.nanvar(flattened, ddof=1)),
                }
            )
            pair_groups: dict[float, list[float]] = defaultdict(list)
            coords = np.asarray([(r, c) for r in range(8) for c in range(8)], dtype=float)
            for first in range(64):
                for second in range(first + 1, 64):
                    lag = float(np.linalg.norm(coords[first] - coords[second]))
                    if np.isfinite(flattened[first]) and np.isfinite(flattened[second]):
                        pair_groups[round(lag, 6)].append(
                            0.5 * (flattened[first] - flattened[second]) ** 2
                        )
            for lag, semivariances in sorted(pair_groups.items()):
                variogram_rows.append(
                    {
                        "metric": metric,
                        "dataset": dataset_label,
                        "lag_grid_units": lag,
                        "pairs": len(semivariances),
                        "semivariance": float(np.mean(semivariances)),
                        "semivariance_over_sample_variance": float(
                            np.mean(semivariances) / np.nanvar(flattened, ddof=1)
                        )
                        if np.nanvar(flattened, ddof=1) > 0.0
                        else float("nan"),
                    }
                )
    spatial_q = benjamini_hochberg(
        [float(row["Moran_permutation_p_two_sided"]) for row in spatial_rows]
    )
    for row, q_value in zip(spatial_rows, spatial_q, strict=True):
        row["Moran_BH_q_across_spatial_tests"] = q_value
    write_csv(RESULTS / "spatial_statistics.csv", spatial_rows)
    write_csv(RESULTS / "spatial_variograms.csv", variogram_rows)

    global_rows: list[dict] = []
    for metric in (
        "force_20nm_pN",
        "force_50nm_pN",
        "apparent_lambda_D_nm",
        "apparent_zeta_magnitude_mV",
        "model_free_local_decay_30_80nm",
        "model_free_normalized_force_50nm",
    ):
        values = parameter_arrays[metric]
        x = values[0]
        y = values[3]
        valid = np.isfinite(x) & np.isfinite(y)
        x = x[valid]
        y = y[valid]
        ordinary = stats.linregress(x, y) if x.size >= 4 and np.ptp(x) > 0.0 else None
        robust = stats.theilslopes(y, x, alpha=0.95) if x.size >= 4 and np.ptp(x) > 0.0 else None
        scale = float(np.dot(x, y) / np.dot(x, x)) if np.dot(x, x) > 0.0 else float("nan")
        zero_residual = y - scale * x
        affine_residual = y - (ordinary.intercept + ordinary.slope * x) if ordinary else np.full_like(y, np.nan)
        rho, rho_p, _ = stable_spearman(x, y)
        global_rows.append(
            {
                "metric": metric,
                "paired_pixels": int(x.size),
                "spatial_spearman_rho": rho,
                "spatial_spearman_p_naive": rho_p,
                "OLS_affine_slope": float(ordinary.slope) if ordinary else float("nan"),
                "OLS_affine_intercept": float(ordinary.intercept) if ordinary else float("nan"),
                "OLS_affine_r2": float(ordinary.rvalue**2) if ordinary else float("nan"),
                "Theil_Sen_slope": float(robust.slope) if robust else float("nan"),
                "Theil_Sen_slope_95CI_low": float(robust.low_slope) if robust else float("nan"),
                "Theil_Sen_slope_95CI_high": float(robust.high_slope) if robust else float("nan"),
                "zero_intercept_scale_map5_over_map2": scale,
                "zero_intercept_scaling_RMSE": float(np.sqrt(np.mean(zero_residual**2))),
                "affine_residual_RMSE": float(np.sqrt(np.nanmean(affine_residual**2))),
                "interpretation": "high correlation plus small residual supports global scaling; residual spatial structure supports patch change",
            }
        )
    write_csv(RESULTS / "global_scaling_vs_spatial_change.csv", global_rows)

    qc_arrays = {
        column: data["qc"][column].to_numpy(dtype=np.float64).reshape(4, 64)
        for column in (
            "local_contact_InvOLS_nm_per_V",
            "far_slope_pN_per_100nm",
            "far_noise_pN",
            "terminal_load_nN",
            "gap_speed_20_200nm_um_per_s",
            "retract_pull_off_force_nN",
            "retract_detachment_piezo_travel_nm",
        )
    }
    qc_arrays["contact_height_nm"] = contact_heights_nm
    association_rows: list[dict] = []
    for response_name in (
        "force_20nm_pN",
        "force_50nm_pN",
        "apparent_lambda_D_nm",
        "apparent_zeta_magnitude_mV",
        "best_discrete_contact_offset_nm",
        "model_free_normalized_force_50nm",
    ):
        response_change = parameter_arrays[response_name][3] - parameter_arrays[response_name][0]
        for qc_name, qc_values in qc_arrays.items():
            qc_change = qc_values[3] - qc_values[0]
            rho, p_value, count = stable_spearman(response_change, qc_change)
            association_rows.append(
                {
                    "response_map5_minus_map2": response_name,
                    "QC_map5_minus_map2": qc_name,
                    "paired_pixels": count,
                    "Spearman_rho": rho,
                    "Spearman_p_naive_pixels": p_value,
                    "claim_status": "association_diagnostic_not_causal",
                }
            )
    association_q = benjamini_hochberg(
        [float(row["Spearman_p_naive_pixels"]) for row in association_rows]
    )
    for row, q_value in zip(association_rows, association_q, strict=True):
        row["Spearman_BH_q_across_QC_associations"] = q_value
    write_csv(RESULTS / "contact_and_QC_change_associations.csv", association_rows)

    figure, axes = plt.subplots(3, 5, figsize=(18.0, 10.0), layout="constrained")
    for row_index, (metric, label) in enumerate(
        (
            ("apparent_lambda_D_nm", "Apparent lambda_D (nm)"),
            ("apparent_zeta_magnitude_mV", "Apparent |psi| (mV)"),
            ("force_50nm_pN", "Force at 50 nm (pN)"),
        )
    ):
        values = parameter_arrays[metric]
        combined = values[np.isfinite(values)]
        vmin, vmax = np.quantile(combined, [0.02, 0.98])
        for map_index, map_order in enumerate(range(2, 6)):
            matrix = physical_matrix(values[map_index], rows, columns)
            image = axes[row_index, map_index].imshow(
                matrix, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax
            )
            axes[row_index, map_index].set_title(f"map{map_order}")
            figure.colorbar(image, ax=axes[row_index, map_index], shrink=0.76)
        slopes = np.asarray(
            [
                row["fitted_map5_minus_map2_change"]
                for row in slope_rows
                if row["metric"] == metric
            ],
            dtype=np.float64,
        )
        limit = np.nanquantile(np.abs(slopes), 0.98)
        image = axes[row_index, 4].imshow(
            physical_matrix(slopes, rows, columns),
            origin="lower",
            cmap="coolwarm",
            vmin=-limit,
            vmax=limit,
        )
        axes[row_index, 4].set_title("four-timepoint fitted change")
        figure.colorbar(image, ax=axes[row_index, 4], shrink=0.76)
        axes[row_index, 0].set_ylabel(label)
    for axis in axes.ravel():
        axis.set_xlabel("physical column")
        axis.set_ylabel(axis.get_ylabel() or "physical row")
    figure.suptitle("Same-pixel parameter maps and temporal-change maps")
    figure.savefig(FIGURES / "same_pixel_parameter_maps_and_slopes.png", dpi=220)
    plt.close(figure)
    return {
        "parameter_arrays": parameter_arrays,
        "paired_rows": paired_rows,
        "spatial_rows": spatial_rows,
        "global_rows": global_rows,
        "association_rows": association_rows,
    }


def select_hierarchical_parameters(map_grid_scores: np.ndarray) -> dict[str, list[int]]:
    """Select parameter-grid indices for four shared-parameter structures."""
    if map_grid_scores.shape != (ZETA_GRID_MV.size, LAMBDA_GRID_NM.size, 4):
        raise ValueError("hierarchical score cube has an unexpected shape")
    flattened = map_grid_scores.reshape(-1, 4)
    all_index = int(np.argmin(np.sum(flattened, axis=1)))
    map_indices = [int(np.argmin(flattened[:, map_index])) for map_index in range(4)]

    lambda_objective = np.sum(np.min(map_grid_scores, axis=0), axis=1)
    shared_lambda_index = int(np.argmin(lambda_objective))
    shared_lambda_map_zeta = [
        int(np.argmin(map_grid_scores[:, shared_lambda_index, map_index]))
        * LAMBDA_GRID_NM.size
        + shared_lambda_index
        for map_index in range(4)
    ]

    zeta_objective = np.sum(np.min(map_grid_scores, axis=1), axis=1)
    shared_zeta_index = int(np.argmin(zeta_objective))
    map_lambda_shared_zeta = [
        shared_zeta_index * LAMBDA_GRID_NM.size
        + int(np.argmin(map_grid_scores[shared_zeta_index, :, map_index]))
        for map_index in range(4)
    ]
    return {
        "all_shared_lambda_and_zeta": [all_index] * 4,
        "map_shared_lambda_and_zeta": map_indices,
        "shared_lambda_map_specific_zeta": shared_lambda_map_zeta,
        "map_specific_lambda_shared_zeta": map_lambda_shared_zeta,
    }


def hierarchical_analysis(data: dict, model_results: dict) -> dict:
    primary = model_results["fits"][PRIMARY_MODEL]
    scores = primary["scores_pN2"]
    candidate_count = scores.shape[0]
    score_by_map_pixel = scores.reshape(candidate_count, 4, 64)
    map_scores = np.sum(score_by_map_pixel, axis=2)
    map_grid_scores = map_scores.reshape(ZETA_GRID_MV.size, LAMBDA_GRID_NM.size, 4)
    selected = select_hierarchical_parameters(map_grid_scores)

    total_observations = 4 * 64 * data["distance_nm"].size
    model_rows: list[dict] = []
    parameter_rows: list[dict] = []
    structures = {
        "independent_curve_lambda_zeta": {
            "rss": float(np.sum(np.min(scores, axis=0))),
            "parameters": 3 * 256,
            "indices": None,
        },
        "all_shared_lambda_and_zeta": {
            "parameters": 256 + 2,
            "indices": selected["all_shared_lambda_and_zeta"],
        },
        "map_shared_lambda_and_zeta": {
            "parameters": 256 + 8,
            "indices": selected["map_shared_lambda_and_zeta"],
        },
        "shared_lambda_map_specific_zeta": {
            "parameters": 256 + 5,
            "indices": selected["shared_lambda_map_specific_zeta"],
        },
        "map_specific_lambda_shared_zeta": {
            "parameters": 256 + 5,
            "indices": selected["map_specific_lambda_shared_zeta"],
        },
    }
    for name, structure in structures.items():
        indices = structure["indices"]
        if indices is not None:
            rss = float(
                sum(map_scores[int(indices[map_index]), map_index] for map_index in range(4))
            )
            structure["rss"] = rss
            for map_index, index in enumerate(indices):
                parameter_rows.append(
                    {
                        "joint_model": name,
                        "map_order": map_index + 2,
                        "relative_time_min": float(data["times_min"][map_index]),
                        "lambda_D_nm": float(FLAT_LAMBDA_NM[int(index)]),
                        "zeta_magnitude_mV": float(FLAT_ZETA_MV[int(index)]),
                    }
                )
        rss = float(structure["rss"])
        parameters = int(structure["parameters"])
        safe_rss = max(rss, np.finfo(np.float64).tiny)
        model_rows.append(
            {
                "joint_model": name,
                "distance_bins": data["distance_nm"].size,
                "curves": 256,
                "nominal_observations_treating_bins_as_independent": total_observations,
                "free_parameters_including_256_curve_baselines": parameters,
                "total_RSS_pN2": rss,
                "global_RMSE_pN": math.sqrt(rss / total_observations),
                "pseudo_AIC_correlated_bins": total_observations
                * math.log(safe_rss / total_observations)
                + 2.0 * parameters,
                "pseudo_BIC_correlated_bins": total_observations
                * math.log(safe_rss / total_observations)
                + math.log(total_observations) * parameters,
                "information_criterion_status": "diagnostic_only_distance_bins_correlated",
            }
        )

    observed = data["force_pN"]
    hyd = model_results["sphere_plane_hyd_pN"].reshape(4, 64, 37)
    library = model_results["nonlinear_library"]
    far_mask = data["distance_nm"] >= 150.0
    test_mask = data["distance_nm"] < 150.0
    cv_rows: list[dict] = []
    cv_structures = tuple(selected)
    for point in range(64):
        training_map_scores = map_scores - score_by_map_pixel[:, :, point]
        training_grid = training_map_scores.reshape(
            ZETA_GRID_MV.size, LAMBDA_GRID_NM.size, 4
        )
        fold_selection = select_hierarchical_parameters(training_grid)
        for model_name in cv_structures:
            indices = fold_selection[model_name]
            for map_index, index in enumerate(indices):
                adjusted = observed[map_index, point] - hyd[map_index, point]
                baseline = float(np.mean(adjusted[far_mask] - library[index, far_mask]))
                prediction = library[index] + hyd[map_index, point] + baseline
                residual = observed[map_index, point, test_mask] - prediction[test_mask]
                cv_rows.append(
                    {
                        "joint_model": model_name,
                        "held_out_point_index": point,
                        "held_out_row": int(data["rows"][point]),
                        "held_out_column": int(data["columns"][point]),
                        "map_order": map_index + 2,
                        "training_pixels_per_map": 63,
                        "held_out_baseline_fit_window_nm": "150-200",
                        "held_out_prediction_window_nm": "20-145",
                        "selected_lambda_D_nm": float(FLAT_LAMBDA_NM[index]),
                        "selected_zeta_magnitude_mV": float(FLAT_ZETA_MV[index]),
                        "held_out_baseline_pN": baseline,
                        "held_out_shape_RMSE_pN": float(np.sqrt(np.mean(residual**2))),
                    }
                )
    cv_frame = pd.DataFrame(cv_rows)
    for row in model_rows:
        selected_cv = cv_frame[cv_frame["joint_model"] == row["joint_model"]]
        row["leave_one_physical_pixel_out_shape_RMSE_mean_pN"] = (
            float(selected_cv["held_out_shape_RMSE_pN"].mean())
            if not selected_cv.empty
            else float("nan")
        )
        row["leave_one_physical_pixel_out_shape_RMSE_median_pN"] = (
            float(selected_cv["held_out_shape_RMSE_pN"].median())
            if not selected_cv.empty
            else float("nan")
        )
    write_csv(RESULTS / "hierarchical_joint_model_comparison.csv", model_rows)
    write_csv(RESULTS / "hierarchical_joint_selected_parameters.csv", parameter_rows)
    write_csv(RESULTS / "hierarchical_leave_one_pixel_out.csv", cv_rows)

    comparison = pd.DataFrame(model_rows)
    plot_frame = comparison[comparison["joint_model"] != "independent_curve_lambda_zeta"].copy()
    figure, axes = plt.subplots(1, 2, figsize=(13.0, 4.8))
    labels = [
        value.replace("_lambda_and_zeta", "").replace("_specific", "-specific").replace("_", "\n")
        for value in plot_frame["joint_model"]
    ]
    axes[0].bar(labels, plot_frame["pseudo_BIC_correlated_bins"], color="#457b9d")
    axes[0].set_ylabel("Pseudo-BIC (correlated-bin diagnostic)")
    axes[1].bar(
        labels,
        plot_frame["leave_one_physical_pixel_out_shape_RMSE_mean_pN"],
        color="#e76f51",
    )
    axes[1].set_ylabel("Leave-one-pixel-out shape RMSE (pN)")
    for axis in axes:
        axis.tick_params(axis="x", labelsize=8)
        axis.grid(axis="y", alpha=0.2)
    figure.suptitle("Joint parameter-sharing models")
    figure.tight_layout()
    figure.savefig(FIGURES / "hierarchical_joint_model_comparison.png", dpi=230)
    plt.close(figure)
    return {"model_rows": model_rows, "parameter_rows": parameter_rows, "cv_rows": cv_rows}


def runs_test(residual: np.ndarray) -> tuple[int, float, float]:
    values = np.asarray(residual, dtype=np.float64)
    signs = values >= 0.0
    n1 = int(np.count_nonzero(signs))
    n2 = int(signs.size - n1)
    runs = int(1 + np.count_nonzero(signs[1:] != signs[:-1]))
    if n1 == 0 or n2 == 0 or signs.size < 3:
        return runs, float("nan"), float("nan")
    mean = 1.0 + 2.0 * n1 * n2 / signs.size
    variance = (
        2.0
        * n1
        * n2
        * (2.0 * n1 * n2 - signs.size)
        / (signs.size**2 * (signs.size - 1.0))
    )
    if variance <= 0.0:
        return runs, float("nan"), float("nan")
    z = (runs - mean) / math.sqrt(variance)
    return runs, float(z), float(2.0 * stats.norm.sf(abs(z)))


def residual_analysis(data: dict, model_results: dict, rng: np.random.Generator) -> dict:
    primary = model_results["fits"][PRIMARY_MODEL]
    residual = primary["residual_pN"].reshape(4, 64, 37)
    diagnostic_rows: list[dict] = []
    for map_index, map_order in enumerate(range(2, 6)):
        for point in range(64):
            values = residual[map_index, point]
            lag_rho = (
                float(stats.pearsonr(values[:-1], values[1:]).statistic)
                if np.ptp(values[:-1]) > 0.0 and np.ptp(values[1:]) > 0.0
                else float("nan")
            )
            denominator = float(np.sum(values**2))
            durbin_watson = (
                float(np.sum(np.diff(values) ** 2) / denominator)
                if denominator > 0.0
                else float("nan")
            )
            runs, z, p_value = runs_test(values)
            diagnostic_rows.append(
                {
                    "map_order": map_order,
                    "relative_time_min": float(data["times_min"][map_index]),
                    "point_index": point,
                    "row": int(data["rows"][point]),
                    "column": int(data["columns"][point]),
                    "rmse_pN": float(primary["rmse_pN"][map_index * 64 + point]),
                    "mean_residual_pN": float(np.mean(values)),
                    "median_residual_pN": float(np.median(values)),
                    "lag1_residual_Pearson_rho": lag_rho,
                    "Durbin_Watson": durbin_watson,
                    "residual_sign_runs": runs,
                    "runs_test_z": z,
                    "runs_test_p_naive_correlated_bins": p_value,
                }
            )
    distance_rows: list[dict] = []
    for map_index, map_order in enumerate(range(2, 6)):
        for distance_index, distance_nm in enumerate(data["distance_nm"]):
            values = residual[map_index, :, distance_index]
            test = stats.ttest_1samp(values, 0.0)
            distance_rows.append(
                {
                    "map_order": map_order,
                    "relative_time_min": float(data["times_min"][map_index]),
                    "distance_nm": float(distance_nm),
                    "residual_mean_pN": float(np.mean(values)),
                    "residual_median_pN": float(np.median(values)),
                    "residual_q25_pN": float(np.quantile(values, 0.25)),
                    "residual_q75_pN": float(np.quantile(values, 0.75)),
                    "residual_std_pN": float(np.std(values, ddof=1)),
                    "one_sample_t_p_naive_pixels": float(test.pvalue),
                    "exact_physical_row_signflip_p": exact_row_signflip_p(
                        values, data["rows"]
                    ),
                }
            )
    distance_q = benjamini_hochberg(
        [float(row["exact_physical_row_signflip_p"]) for row in distance_rows]
    )
    for row, q_value in zip(distance_rows, distance_q, strict=True):
        row["exact_row_signflip_BH_q_across_distance_tests"] = q_value
    write_csv(RESULTS / "residual_diagnostics_per_curve.csv", diagnostic_rows)
    write_csv(RESULTS / "residual_by_distance.csv", distance_rows)

    figure, axes = plt.subplots(2, 4, figsize=(17.0, 8.0), layout="constrained")
    for map_index, map_order in enumerate(range(2, 6)):
        values = residual[map_index]
        median = np.median(values, axis=0)
        q25 = np.quantile(values, 0.25, axis=0)
        q75 = np.quantile(values, 0.75, axis=0)
        axes[0, map_index].fill_between(data["distance_nm"], q25, q75, color="#8ecae6", alpha=0.4)
        axes[0, map_index].plot(data["distance_nm"], median, color="#023047", lw=2.0)
        axes[0, map_index].axhline(0.0, color="0.4", lw=0.8)
        axes[0, map_index].set_title(f"map{map_order}")
        axes[0, map_index].set_xlabel("D (nm)")
        axes[0, map_index].set_ylabel("Residual (pN)")
        axes[0, map_index].grid(alpha=0.2)
        rmse = primary["rmse_pN"].reshape(4, 64)[map_index]
        image = axes[1, map_index].imshow(
            physical_matrix(rmse, data["rows"], data["columns"]),
            origin="lower",
            cmap="magma",
        )
        axes[1, map_index].set_xlabel("physical column")
        axes[1, map_index].set_ylabel("physical row")
        figure.colorbar(image, ax=axes[1, map_index], label="RMSE (pN)", shrink=0.77)
    figure.suptitle("Sphere–plane nonlinear-PB residual structure")
    figure.savefig(FIGURES / "sphere_plane_residual_diagnostics.png", dpi=230)
    plt.close(figure)
    return {"diagnostic_rows": diagnostic_rows, "distance_rows": distance_rows}


def fit_time_laws(times: np.ndarray, values: np.ndarray, metric: str) -> list[dict]:
    times = np.asarray(times, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    if times.shape != (4,) or values.shape != (4,) or not np.all(np.isfinite(values)):
        return []

    def exponential_zero(t: np.ndarray, amplitude: float, tau: float) -> np.ndarray:
        return amplitude * np.exp(-t / tau)

    def exponential_offset(
        t: np.ndarray, asymptote: float, amplitude: float, tau: float
    ) -> np.ndarray:
        return asymptote + amplitude * np.exp(-t / tau)

    models: list[tuple[str, int, np.ndarray, str]] = []
    ordinary = stats.linregress(times, values)
    models.append(
        (
            "linear",
            2,
            ordinary.intercept + ordinary.slope * times,
            json.dumps({"intercept": ordinary.intercept, "slope_per_min": ordinary.slope}),
        )
    )
    try:
        parameters, _ = curve_fit(
            exponential_zero,
            times,
            values,
            p0=(max(values[0], 1e-6), max(np.ptp(times) / 2.0, 1.0)),
            bounds=([0.0, 0.1], [np.inf, 1e5]),
            maxfev=20000,
        )
        models.append(
            (
                "exponential_to_zero",
                2,
                exponential_zero(times, *parameters),
                json.dumps({"amplitude": parameters[0], "tau_min": parameters[1]}),
            )
        )
    except (RuntimeError, ValueError, FloatingPointError):
        pass
    try:
        parameters, _ = curve_fit(
            exponential_offset,
            times,
            values,
            p0=(max(values[-1] * 0.5, 0.0), max(values[0] - values[-1] * 0.5, 1e-6), max(np.ptp(times) / 2.0, 1.0)),
            bounds=([-2.0 * np.max(np.abs(values)), 0.0, 0.1], [2.0 * np.max(np.abs(values)), np.inf, 1e5]),
            maxfev=40000,
        )
        models.append(
            (
                "exponential_with_asymptote",
                3,
                exponential_offset(times, *parameters),
                json.dumps(
                    {
                        "asymptote": parameters[0],
                        "amplitude": parameters[1],
                        "tau_min": parameters[2],
                    }
                ),
            )
        )
    except (RuntimeError, ValueError, FloatingPointError):
        pass

    rows: list[dict] = []
    for name, parameters, prediction, parameter_json in models:
        rss = float(np.sum((values - prediction) ** 2))
        safe = max(rss, np.finfo(np.float64).tiny)
        aic = 4.0 * math.log(safe / 4.0) + 2.0 * parameters
        aicc = (
            aic + 2.0 * parameters * (parameters + 1.0) / (4.0 - parameters - 1.0)
            if 4 > parameters + 1
            else float("nan")
        )
        rows.append(
            {
                "metric": metric,
                "time_model": name,
                "timepoints": 4,
                "free_parameters": parameters,
                "parameters_json": parameter_json,
                "rss": rss,
                "rmse": math.sqrt(rss / 4.0),
                "AIC_four_points_diagnostic": aic,
                "AICc_four_points_diagnostic": aicc,
                "claim_status": "exploratory_four_timepoints_cannot_identify_kinetic_law",
            }
        )
    return rows


def time_law_analysis(data: dict, paired_results: dict) -> list[dict]:
    arrays = paired_results["parameter_arrays"]
    metric_arrays = {
        "apparent_lambda_D_nm_map_median": arrays["apparent_lambda_D_nm"],
        "apparent_zeta_magnitude_mV_map_median": arrays[
            "apparent_zeta_magnitude_mV"
        ],
        "force_20nm_pN_map_median": arrays["force_20nm_pN"],
        "force_50nm_pN_map_median": arrays["force_50nm_pN"],
        "model_free_local_decay_30_80nm_map_median": arrays[
            "model_free_local_decay_30_80nm"
        ],
    }
    rows: list[dict] = []
    for metric, values in metric_arrays.items():
        medians = np.asarray([finite_quantile(item, 0.5) for item in values])
        rows.extend(fit_time_laws(data["times_min"], medians, metric))
    write_csv(RESULTS / "four_timepoint_time_law_diagnostics.csv", rows)
    return rows


def calibration_scale_sensitivity(data: dict, model_results: dict) -> list[dict]:
    relative_sd = data["spring_sd_N_per_m"] / data["spring_N_per_m"]
    factors = (1.0 - relative_sd, 1.0, 1.0 + relative_sd)
    library = model_results["nonlinear_library"]
    speed = data["speeds_um_per_s"].reshape(256)
    hyd = hydrodynamic_force_pN(data["distance_nm"], speed, RADIUS_M)
    rows: list[dict] = []
    for factor in factors:
        fit = batch_grid_fit(
            data["force_pN"].reshape(256, 37) * factor,
            library,
            hyd,
        )
        for map_index, map_order in enumerate(range(2, 6)):
            selected = slice(map_index * 64, (map_index + 1) * 64)
            rows.append(
                {
                    "observed_force_common_scale_factor": factor,
                    "spring_constant_interpretation": "k_minus_1SD"
                    if factor < 1.0
                    else "k_plus_1SD"
                    if factor > 1.0
                    else "nominal_k",
                    "map_order": map_order,
                    "lambda_D_nm_median": finite_quantile(fit["lambda_D_nm"][selected], 0.5),
                    "zeta_magnitude_mV_median": finite_quantile(fit["zeta_magnitude_mV"][selected], 0.5),
                    "baseline_pN_median": finite_quantile(fit["baseline_pN"][selected], 0.5),
                    "rmse_pN_median": finite_quantile(fit["rmse_pN"][selected], 0.5),
                    "uncertainty_scope": "common_force_scale_only; InvOLS distance-axis and contact-zero handled separately",
                }
            )
    write_csv(RESULTS / "calibration_common_force_scale_sensitivity.csv", rows)
    return rows


def numerical_audit(
    data: dict,
    model_free: dict,
    model_results: dict,
    contact_results: dict,
    hierarchy: dict,
) -> dict:
    distance = data["distance_nm"]
    library = model_results["nonlinear_library"]
    identity_errors: list[float] = []
    for zeta_index, lambda_index in ((0, 0), (50, 80), (100, 130), (166, 180)):
        flat_index = zeta_index * LAMBDA_GRID_NM.size + lambda_index
        reference = base.total_equilibrium_force_pN(
            distance,
            float(LAMBDA_GRID_NM[lambda_index]),
            float(ZETA_GRID_MV[zeta_index]),
            0.0,
            EPSILON_R,
            "nonlinear_pb_derjaguin",
        )
        identity_errors.append(float(np.max(np.abs(library[flat_index] - reference))))
    formula_identity_error = max(identity_errors)

    limiting_distance_nm = np.asarray([20.0, 50.0, 100.0, 200.0])
    nonlinear_small_potential = base.edl_force_pN(
        limiting_distance_nm,
        20.0,
        0.1,
        EPSILON_R,
        "nonlinear_pb_derjaguin",
    )
    linear_small_potential = base.edl_force_pN(
        limiting_distance_nm,
        20.0,
        0.1,
        EPSILON_R,
        "linear_hhf_equal_potential",
    )
    low_potential_limit_max_relative_error = float(
        np.max(
            np.abs(nonlinear_small_potential - linear_small_potential)
            / np.maximum(np.abs(linear_small_potential), 1e-30)
        )
    )

    primary = model_results["fits"][PRIMARY_MODEL]
    observed = data["force_pN"].reshape(256, 37)
    prediction_closure = float(
        np.max(
            np.abs(
                primary["residual_pN"]
                - (observed - primary["prediction_pN"])
            )
        )
    )
    normalized = model_free["normalized"]
    valid = model_free["normalization_valid"]
    index20 = int(np.flatnonzero(np.isclose(distance, 20.0))[0])
    index200 = int(np.flatnonzero(np.isclose(distance, 200.0))[0])
    normalization_endpoint_error = max(
        float(np.nanmax(np.abs(normalized[:, :, index20][valid] - 1.0))),
        float(np.nanmax(np.abs(normalized[:, :, index200][valid]))),
    )

    hierarchy_frame = pd.DataFrame(hierarchy["model_rows"]).set_index("joint_model")
    independent_rss = float(
        hierarchy_frame.loc["independent_curve_lambda_zeta", "total_RSS_pN2"]
    )
    shared_rss = hierarchy_frame.drop(index="independent_curve_lambda_zeta")[
        "total_RSS_pN2"
    ].to_numpy(dtype=np.float64)
    hierarchy_order_pass = bool(np.all(independent_rss <= shared_rss + 1e-6))

    median_speed = float(np.median(data["speeds_um_per_s"]))
    hyd_probe = hydrodynamic_force_pN(
        np.asarray([20.0, 50.0, 100.0, 200.0]),
        np.asarray([median_speed]),
        RADIUS_M,
    )[0]
    hyd_at_one = hydrodynamic_force_pN(
        limiting_distance_nm, np.asarray([1.0]), RADIUS_M
    )[0]
    hyd_at_two = hydrodynamic_force_pN(
        limiting_distance_nm, np.asarray([2.0]), RADIUS_M
    )[0]
    hyd_inverse_distance_relative_range = float(
        np.ptp(hyd_at_one * limiting_distance_nm)
        / np.mean(hyd_at_one * limiting_distance_nm)
    )
    hyd_linear_speed_max_relative_error = float(
        np.max(np.abs(hyd_at_two - 2.0 * hyd_at_one) / hyd_at_two)
    )
    dimensionless_potential = (
        elementary_charge
        * primary["zeta_magnitude_mV"]
        * 1e-3
        / (Boltzmann * TEMPERATURE_K)
    )
    checks = {
        "claimed_output": "model-light same-pixel shape/spatial evidence plus apparent sphere-plane PB parameters",
        "input_domain": "four same-location 8x8 maps, 256 curves, 20-200 nm in 5 nm bins, 25.6 C pure-water environment",
        "units": "distance nm, force pN, potential mV, time min, k N/m, InvOLS nm/V",
        "acceptable_numerical_error": "formula/prediction closure <=1e-7 pN; qualitative inference uses spatial and model sensitivity",
        "reproducibility": f"deterministic except seeded spatial permutation/bootstrap; seed={RANDOM_SEED}",
        "nonlinear_library_formula_identity_max_abs_error_pN": formula_identity_error,
        "nonlinear_to_linear_HHF_0p1mV_limit_max_relative_error": low_potential_limit_max_relative_error,
        "hydrodynamic_inverse_distance_relative_range": hyd_inverse_distance_relative_range,
        "hydrodynamic_linear_speed_max_relative_error": hyd_linear_speed_max_relative_error,
        "prediction_residual_closure_max_abs_error_pN": prediction_closure,
        "model_free_normalization_endpoint_max_abs_error": normalization_endpoint_error,
        "all_primary_parameters_finite": bool(
            all(
                np.all(np.isfinite(primary[field]))
                for field in (
                    "lambda_D_nm",
                    "zeta_magnitude_mV",
                    "baseline_pN",
                    "rmse_pN",
                    "r2",
                )
            )
        ),
        "hierarchical_independent_RSS_not_above_shared_models": hierarchy_order_pass,
        "contact_offset_search_nm": CONTACT_OFFSETS_NM.tolist(),
        "contact_profile_boundary_fraction": float(
            np.mean(
                [
                    bool(row["best_offset_at_search_boundary"])
                    for row in contact_results["nuisance_rows"]
                ]
            )
        ),
        "fixed_sphere_plane_no_slip_hydrodynamic_force_at_median_speed_pN": {
            str(distance_nm): float(force_pN)
            for distance_nm, force_pN in zip(
                (20.0, 50.0, 100.0, 200.0), hyd_probe, strict=True
            )
        },
        "dimensionless_potential_epsi_over_kT_median": float(
            np.median(dimensionless_potential)
        ),
        "dimensionless_potential_epsi_over_kT_q25": float(
            np.quantile(dimensionless_potential, 0.25)
        ),
        "dimensionless_potential_epsi_over_kT_q75": float(
            np.quantile(dimensionless_potential, 0.75)
        ),
        "grid": {
            "lambda_points": int(LAMBDA_GRID_NM.size),
            "lambda_range_nm": [float(LAMBDA_GRID_NM[0]), float(LAMBDA_GRID_NM[-1])],
            "zeta_points": int(ZETA_GRID_MV.size),
            "zeta_range_mV": [float(ZETA_GRID_MV[0]), float(ZETA_GRID_MV[-1])],
            "distance_bins": int(distance.size),
        },
    }
    core_pass = bool(
        formula_identity_error <= 1e-7
        and low_potential_limit_max_relative_error <= 5e-3
        and hyd_inverse_distance_relative_range <= 1e-12
        and hyd_linear_speed_max_relative_error <= 1e-12
        and prediction_closure <= 1e-7
        and normalization_endpoint_error <= 1e-12
        and checks["all_primary_parameters_finite"]
        and hierarchy_order_pass
    )
    checks["core_numerical_checks_pass"] = core_pass
    if not core_pass:
        raise AssertionError(f"core numerical audit failed: {checks}")

    lines = [
        "# Numerical audit: 31-08-26 advanced analysis",
        "",
        "## Computational claim",
        "",
        f"- Claimed output: {checks['claimed_output']}.",
        f"- Domain: {checks['input_domain']}.",
        f"- Units: {checks['units']}.",
        f"- Numerical contract: {checks['acceptable_numerical_error']}.",
        f"- Reproducibility: {checks['reproducibility']}.",
        "",
        "## Direct probes",
        "",
        f"- Nonlinear-PB library identity against the existing sphere-plane evaluator: max `{formula_identity_error:.3e} pN` — PASS.",
        f"- Nonlinear PB to linear-HHF limit at `|psi|=0.1 mV`, `lambda=20 nm`: max relative difference `{low_potential_limit_max_relative_error:.3e}` — PASS.",
        f"- Fixed lubrication limiting laws: relative range of `F*D` `{hyd_inverse_distance_relative_range:.3e}`; max relative error in `F(2v)=2F(v)` `{hyd_linear_speed_max_relative_error:.3e}` — PASS.",
        f"- Saved prediction/residual algebraic closure: max `{prediction_closure:.3e} pN` — PASS.",
        f"- Endpoint normalization closure, valid curves: max `{normalization_endpoint_error:.3e}` — PASS.",
        f"- Independent-fit RSS is no larger than every parameter-sharing model: `{hierarchy_order_pass}` — PASS.",
        f"- All 256 primary fitted parameter and diagnostic arrays finite: `{checks['all_primary_parameters_finite']}` — PASS.",
        "",
        "## Physical-scale checks",
        "",
        f"At median measured gap speed `{median_speed:.6f} µm/s`, the fixed ideal no-slip sphere-plane lubrication term is "
        + ", ".join(
            f"{distance_nm:g} nm: {force_pN:.2f} pN"
            for distance_nm, force_pN in zip(
                (20.0, 50.0, 100.0, 200.0), hyd_probe, strict=True
            )
        )
        + ".",
        f"The fitted dimensionless potential `e|psi|/kBT` has median `{np.median(dimensionless_potential):.3f}` "
        f"[IQR `{np.quantile(dimensionless_potential,0.25):.3f}`, `{np.quantile(dimensionless_potential,0.75):.3f}`]. "
        "The primary evaluator is nonlinear PB, so this is a regime descriptor rather than a linearization gate.",
        "",
        "## Claim limits",
        "",
        "- Adjacent 5 nm bins from one curve are correlated; profile regions and AIC/BIC are diagnostic, not formal independent-observation confidence statements.",
        "- The nonlinear model assumes equal constant potential, symmetric 1:1 PB electrolyte, ideal Derjaguin geometry, fixed Hamaker constant, and ideal no-slip hydrodynamics.",
        "- Spring-constant common scaling, contact-zero shifts, baseline choice, spatial clustering, and model form are audited separately; no single branch is promoted to a bulk Debye-length or zeta-potential measurement.",
        "- InvOLS uncertainty changes both force and separation and is not represented by the common-force-scale branch alone.",
        "",
    ]
    (RESULTS / "NUMERICS_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")
    return checks


def render_report(
    data: dict,
    model_free: dict,
    model_results: dict,
    contact_results: dict,
    identifiability: dict,
    paired: dict,
    hierarchy: dict,
    residual: dict,
    time_rows: list[dict],
    calibration_rows: list[dict],
    audit: dict,
) -> str:
    primary_summary = [
        row for row in model_results["summary_rows"] if row["model"] == PRIMARY_MODEL
    ]
    paired_by_metric = {row["metric"]: row for row in paired["paired_rows"]}
    shape50 = next(
        row
        for row in model_free["paired_rows"]
        if row["metric"] == "normalized_force"
        and np.isclose(float(row["distance_nm"]), 50.0)
    )
    local_decay_30_80 = next(
        row
        for row in model_free["paired_rows"]
        if row["metric"] == "local_decay_30_80nm"
    )
    contact_summary = pd.DataFrame(contact_results["summary_rows"])
    nuisance = pd.DataFrame(contact_results["nuisance_rows"])
    hierarchy_frame = pd.DataFrame(hierarchy["model_rows"])
    shared_only = hierarchy_frame[
        hierarchy_frame["joint_model"] != "independent_curve_lambda_zeta"
    ]
    best_bic = shared_only.loc[shared_only["pseudo_BIC_correlated_bins"].idxmin()]
    best_cv = shared_only.loc[
        shared_only["leave_one_physical_pixel_out_shape_RMSE_mean_pN"].idxmin()
    ]
    residual_frame = pd.DataFrame(residual["diagnostic_rows"])
    spatial_frame = pd.DataFrame(paired["spatial_rows"])
    association_frame = pd.DataFrame(paired["association_rows"])
    model_summary = pd.DataFrame(model_results["summary_rows"])

    lambda_fixed_shift_changes = []
    for offset in CONTACT_OFFSETS_NM:
        selected = contact_summary[
            np.isclose(contact_summary["fixed_contact_offset_nm"], offset)
        ].sort_values("map_order")
        lambda_fixed_shift_changes.append(
            float(selected.iloc[-1]["lambda_D_nm_median"] - selected.iloc[0]["lambda_D_nm_median"])
        )
    zeta_fixed_shift_changes = []
    for offset in CONTACT_OFFSETS_NM:
        selected = contact_summary[
            np.isclose(contact_summary["fixed_contact_offset_nm"], offset)
        ].sort_values("map_order")
        zeta_fixed_shift_changes.append(
            float(
                selected.iloc[-1]["zeta_magnitude_mV_median"]
                - selected.iloc[0]["zeta_magnitude_mV_median"]
            )
        )
    calibration = pd.DataFrame(calibration_rows)
    calibration_nominal = calibration[calibration["spring_constant_interpretation"] == "nominal_k"]
    calibration_low = calibration[calibration["spring_constant_interpretation"] == "k_minus_1SD"]
    calibration_high = calibration[calibration["spring_constant_interpretation"] == "k_plus_1SD"]

    lines = [
        "# 31-08-26 map2–5 advanced same-pixel analysis",
        "",
        "## 结论先行",
        "",
        "现有数据确认20–50 nm曲线不只是整体幅值下降，归一化后的曲线形状也发生变化；但contact/separation零点、空间斑块、baseline和模型参数耦合仍足以显著改变拟合参数。因此这里报告的是 **apparent/model-conditioned** `lambda_D` 与 `|psi|`，不是已经验证的bulk Debye length或zeta potential。",
        "",
        f"Primary sphere–plane nonlinear-PB fit给出的map2→map5中位apparent `lambda_D`为 `{primary_summary[0]['lambda_D_nm_median']:.2f} → {primary_summary[-1]['lambda_D_nm_median']:.2f} nm`，`|psi|`为 `{primary_summary[0]['zeta_magnitude_mV_median']:.1f} → {primary_summary[-1]['zeta_magnitude_mV_median']:.1f} mV`。",
        "",
        f"Model-free endpoint-normalized force在50 nm的map5−map2 paired median为 `{shape50['map5_minus_map2_median']:+.3f}`，row-signflip `p={shape50['exact_physical_row_signflip_p']:.4g}`，row-cluster 95% CI `[{shape50['row_cluster_bootstrap_median_95CI_low']:+.3f}, {shape50['row_cluster_bootstrap_median_95CI_high']:+.3f}]`；30–80 nm local log-slope得到的decay length变化为 `{local_decay_30_80['map5_minus_map2_median']:+.2f} nm`，但其row-signflip `p={local_decay_30_80['exact_physical_row_signflip_p']:.4g}`且cluster CI跨零，因此该window本身不构成稳健证据。",
        "",
        f"把所有曲线统一平移 `D0=-10…+10 nm` 后，map5−map2中位lambda变化仍处在 `{min(lambda_fixed_shift_changes):+.2f}…{max(lambda_fixed_shift_changes):+.2f} nm`；对应potential变化 `{min(zeta_fixed_shift_changes):+.1f}…{max(zeta_fixed_shift_changes):+.1f} mV`。但允许每条曲线从七档D0中自由选择时，`{np.mean(nuisance['best_offset_at_search_boundary']):.1%}` 落在搜索边界，说明D0与PB参数不能由当前20–200 nm窗口稳定共同识别。",
        "",
        "## Primary sphere–plane nonlinear-PB结果",
        "",
        "模型为equal constant-potential silica sphere–plane nonlinear 1:1 PB Derjaguin + sphere–plane vdW + fixed ideal no-slip lubrication + per-curve constant baseline。固定 `R=4.546849 µm`, `A_H=2.4e-21 J`, `epsilon_r=78.5`, `T=25.6 °C`。",
        "",
        "| map | time (min) | apparent lambda_D (nm), median [IQR] | apparent |psi| (mV), median [IQR] | RMSE median (pN) | R2 median |",
        "|---:|---:|:---|:---|---:|---:|",
    ]
    for row in primary_summary:
        lines.append(
            f"| {row['map_order']} | {row['relative_time_min']:.2f} | "
            f"{row['lambda_D_nm_median']:.2f} [{row['lambda_D_nm_q25']:.2f}, {row['lambda_D_nm_q75']:.2f}] | "
            f"{row['zeta_magnitude_mV_median']:.1f} [{row['zeta_magnitude_mV_q25']:.1f}, {row['zeta_magnitude_mV_q75']:.1f}] | "
            f"{row['rmse_pN_median']:.2f} | {row['r2_median']:.3f} |"
        )
    lines += [
        "",
        "### Same-pixel map2→map5统计",
        "",
        "| metric | map2 median | map5 median | paired delta median [IQR] | row-signflip p | row-cluster 95% CI |",
        "|:---|---:|---:|:---|---:|:---|",
    ]
    for metric in (
        "apparent_lambda_D_nm",
        "apparent_zeta_magnitude_mV",
        "force_20nm_pN",
        "force_50nm_pN",
        "contact_height_nm",
        "local_contact_InvOLS_nm_per_V",
    ):
        row = paired_by_metric[metric]
        lines.append(
            f"| {metric} | {row['map2_median']:.3g} | {row['map5_median']:.3g} | "
            f"{row['map5_minus_map2_median']:+.3g} [{row['difference_q25']:+.3g}, {row['difference_q75']:+.3g}] | "
            f"{row['exact_physical_row_signflip_p']:.4g} | "
            f"[{row['row_cluster_bootstrap_median_95CI_low']:+.3g}, {row['row_cluster_bootstrap_median_95CI_high']:+.3g}] |"
        )
    lines += [
        "",
        "## Contact-zero与参数可识别性",
        "",
        "固定D0敏感性与per-curve离散D0 nuisance结果分别保存在 `contact_zero_fixed_offset_summary.csv` 与 `contact_zero_profiled_nuisance.csv`。固定共同D0不能消除map2→map5趋势；但per-curve D0 profile大量触边，不能把其最优D0直接解释为真实contact drift。",
        "",
        "Map-median及逐曲线lambda–potential profile均使用far-noise或spatial-IQR scale；相邻distance bins相关，所以profile范围不是formal 95% CI。",
        "",
        f"尤其map5的map-median spatial-IQR profile在lambda上达到 `{identifiability['map_profile_summary'][-1]['lambda_profile_high_nm']:.0f} nm` grid上界、potential达到 `{identifiability['map_profile_summary'][-1]['zeta_profile_high_mV']:.0f} mV`上界；这说明用一条map median曲线做formal参数区间是不可靠的。逐曲线profiles通常更窄，但仍受相邻D bins相关影响。",
        "",
        "## Joint/hierarchical parameter-sharing",
        "",
        f"在共享参数模型中，pseudo-BIC最低的是 `{best_bic['joint_model']}`；leave-one-physical-pixel-out shape RMSE最低的也是 `{best_cv['joint_model']}`（mean `{best_cv['leave_one_physical_pixel_out_shape_RMSE_mean_pN']:.2f} pN`）。两项均要求map-specific lambda和potential；一个跨四图完全共享的参数对不足以描述数据。",
        "",
        "每条曲线仍保留自己的constant baseline。CV每次留出同一physical pixel的四条曲线，用其150–200 nm只估baseline，再预测20–145 nm shape；因此没有把held-out lambda/zeta重新拟合回来。",
        "",
        "## 空间结构与QC关联",
        "",
        f"共 `{len(spatial_frame)}` 个map/difference空间检验中，Moran permutation在BH校正后 `q<0.05` 的有 `{int(np.count_nonzero(spatial_frame['Moran_BH_q_across_spatial_tests'] < 0.05))}` 个。64个pixel因此不能普遍视为64个独立重复；报告以row-signflip和row-cluster bootstrap为主。",
        "",
        f"map2→map5的spatial correlation在20 nm force为 `{next(row for row in paired['global_rows'] if row['metric']=='force_20nm_pN')['spatial_spearman_rho']:+.3f}`，到50 nm仅 `{next(row for row in paired['global_rows'] if row['metric']=='force_50nm_pN')['spatial_spearman_rho']:+.3f}`；apparent lambda甚至为 `{next(row for row in paired['global_rows'] if row['metric']=='apparent_lambda_D_nm')['spatial_spearman_rho']:+.3f}`，而potential仍为 `{next(row for row in paired['global_rows'] if row['metric']=='apparent_zeta_magnitude_mV')['spatial_spearman_rho']:+.3f}`。因此变化不是一个统一global scale factor：50 nm和lambda的空间pattern发生了重排。",
        "",
        f"对force/fit-parameter变化与InvOLS、far slope、noise、load、speed、pull-off、contact height变化的 `{len(association_frame)}` 个Spearman诊断中，BH `q<0.05` 的有 `{int(np.count_nonzero(association_frame['Spearman_BH_q_across_QC_associations'] < 0.05))}` 个；这些仍是association，不是因果校正。",
        "",
        f"其中50 nm force变化与far-field slope变化的关联最强：Spearman `rho={association_frame[(association_frame['response_map5_minus_map2']=='force_50nm_pN') & (association_frame['QC_map5_minus_map2']=='far_slope_pN_per_100nm')]['Spearman_rho'].iloc[0]:+.3f}`，BH `q={association_frame[(association_frame['response_map5_minus_map2']=='force_50nm_pN') & (association_frame['QC_map5_minus_map2']=='far_slope_pN_per_100nm')]['Spearman_BH_q_across_QC_associations'].iloc[0]:.3g}`。这说明50 nm绝对force下降中有明显baseline/drift耦合，不能全部解释为EDL变化。",
        "",
        f"相反，endpoint-normalized 50 nm shape变化与far-field slope变化没有显示同样的单调关联：`rho={association_frame[(association_frame['response_map5_minus_map2']=='model_free_normalized_force_50nm') & (association_frame['QC_map5_minus_map2']=='far_slope_pN_per_100nm')]['Spearman_rho'].iloc[0]:+.3f}`，BH `q={association_frame[(association_frame['response_map5_minus_map2']=='model_free_normalized_force_50nm') & (association_frame['QC_map5_minus_map2']=='far_slope_pN_per_100nm')]['Spearman_BH_q_across_QC_associations'].iloc[0]:.3g}`。因此已观察到的normalized shape separation不能仅由“far-slope越变、50 nm force越变”这一单一关系解释，但这并不排除其他baseline curvature或history systematics。",
        "",
        "Raw fitted contact-height scanner coordinate从map2到map5整体增加约427 nm，并保持很强空间相关；该坐标在每条曲线构造separation时已被逐条减去，因此它反映scanner/sample coordinate creep或offset，而不是未校正的427 nm物理gap变化，也不能直接等同于PB fit中的D0。",
        "",
        "## Model comparison与residual",
        "",
    ]
    for model_name in model_summary["model"].unique():
        selected_model = model_summary[model_summary["model"] == model_name]
        lines.append(
            f"- `{model_name}`: four-map median RMSE range "
            f"`{selected_model['rmse_pN_median'].min():.2f}–{selected_model['rmse_pN_median'].max():.2f} pN`; "
            f"per-curve AIC winner fraction range `{selected_model['aic_winner_fraction'].min():.1%}–{selected_model['aic_winner_fraction'].max():.1%}`."
        )
    lines += [
        "",
        f"Primary residual的lag-1 correlation中位数为 `{residual_frame['lag1_residual_Pearson_rho'].median():+.3f}`；naive runs test `p<0.05` 的曲线占 `{np.mean(residual_frame['runs_test_p_naive_correlated_bins'] < 0.05):.1%}`。这表明高R2不能代替residual结构检查。",
        "",
        "## Calibration与time-law敏感性",
        "",
        f"D5 `k` repeatability SD对应common force scale约 `±{data['spring_sd_N_per_m']/data['spring_N_per_m']:.1%}`。在这两个端点下，四图中位lambda相对nominal的最大偏移为 `"
        f"{max(np.max(np.abs(calibration_low['lambda_D_nm_median'].to_numpy()-calibration_nominal['lambda_D_nm_median'].to_numpy())), np.max(np.abs(calibration_high['lambda_D_nm_median'].to_numpy()-calibration_nominal['lambda_D_nm_median'].to_numpy()))):.2f} nm`；"
        f"potential最大偏移 `{max(np.max(np.abs(calibration_low['zeta_magnitude_mV_median'].to_numpy()-calibration_nominal['zeta_magnitude_mV_median'].to_numpy())), np.max(np.abs(calibration_high['zeta_magnitude_mV_median'].to_numpy()-calibration_nominal['zeta_magnitude_mV_median'].to_numpy()))):.1f} mV`。InvOLS还会改变distance axis，不能由此common-scale branch完全代表。",
        "",
        f"四时点对linear、exponential-to-zero及exponential-with-asymptote共生成 `{len(time_rows)}` 个诊断拟合。四点不足以识别kinetic law；tau或asymptote只能作为探索性描述，不能用来证明液体浓度随时间变化。",
        "",
        f"在median gap speed下，固定ideal no-slip sphere–plane hydrodynamic force为20 nm `{audit['fixed_sphere_plane_no_slip_hydrodynamic_force_at_median_speed_pN']['20.0']:.2f} pN`、50 nm `{audit['fixed_sphere_plane_no_slip_hydrodynamic_force_at_median_speed_pN']['50.0']:.2f} pN`。它分别约占map2中位measured force的 `{audit['fixed_sphere_plane_no_slip_hydrodynamic_force_at_median_speed_pN']['20.0']/paired_by_metric['force_20nm_pN']['map2_median']:.1%}` / `{audit['fixed_sphere_plane_no_slip_hydrodynamic_force_at_median_speed_pN']['50.0']/paired_by_metric['force_50nm_pN']['map2_median']:.1%}`，到map5则约为 `{audit['fixed_sphere_plane_no_slip_hydrodynamic_force_at_median_speed_pN']['20.0']/paired_by_metric['force_20nm_pN']['map5_median']:.1%}` / `{audit['fixed_sphere_plane_no_slip_hydrodynamic_force_at_median_speed_pN']['50.0']/paired_by_metric['force_50nm_pN']['map5_median']:.1%}`。这些只是固定no-slip理论项相对line-corrected measured force的量级，不是数据已独立识别出的hyd比例。",
        "",
        "## 数值与物理边界",
        "",
        f"- PB library formula identity max error `{audit['nonlinear_library_formula_identity_max_abs_error_pN']:.3e} pN`; prediction closure `{audit['prediction_residual_closure_max_abs_error_pN']:.3e} pN`; core checks PASS.",
        "- Primary nonlinear PB避免了linear-HHF的small-potential限制，但仍假定symmetric 1:1 electrolyte、equal constant potential、Derjaguin、固定Hamaker和ideal no-slip。纯水中实际离子种类/CO2并未由force curve确定。",
        "- Hydrodynamic term是理论固定项；当前same-speed四图不能从数据本身独立验证其幅度或slip boundary condition。",
        "- AIC/BIC、profile和distance-wise p值均受同一曲线内distance correlation影响，已明确标为diagnostic。",
        "- 接触零点、baseline、共同force scale、空间相关和模型选择分别做了sensitivity branch；没有把任一branch提升为bulk Debye length或zeta potential测量。",
        "",
        "## 主要输出",
        "",
        "- `model_free_*`: endpoint-normalized shape、local log-slope及paired tests。",
        "- `contact_zero_*`: 七档fixed D0及per-curve discrete nuisance profile。",
        "- `surface_model_comparison_*`, `sphere_plane_nonlinear_predictions.csv`: 四模型统一OLS比较和primary逐点预测。",
        "- `same_pixel_*`, `spatial_*`, `global_scaling_vs_spatial_change.csv`: paired trajectories、空间统计和global-vs-patch变化。",
        "- `per_curve_parameter_profile_intervals.csv`, `map_median_parameter_profiles.csv`: identifiability诊断。",
        "- `hierarchical_*`: parameter-sharing、pseudo-information criteria和leave-one-pixel-out shape prediction。",
        "- `residual_*`, `four_timepoint_time_law_diagnostics.csv`, `calibration_common_force_scale_sensitivity.csv`。",
        "- `NUMERICS_AUDIT.md`, `provenance.json`, `artifact_manifest.sha256`, `figures/`。",
        "",
    ]
    return "\n".join(lines)


def create_manifest(paths: list[Path], destination: Path) -> None:
    unique = sorted(set(path.resolve() for path in paths), key=lambda path: str(path))
    destination.write_text(
        "\n".join(
            f"{sha256_file(path)}  {path.relative_to(ROOT)}" for path in unique
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)
    data = load_inputs()
    contact_heights_nm = reconstruct_contact_heights()
    model_free = model_free_analysis(data, rng)
    model_results = fit_model_variants(data)
    contact_results = contact_zero_analysis(data, model_results["nonlinear_library"])
    identifiability = identifiability_analysis(data, model_results)
    paired = paired_parameter_and_spatial_analysis(
        data,
        model_free,
        model_results,
        contact_heights_nm,
        contact_results,
        rng,
    )
    hierarchy = hierarchical_analysis(data, model_results)
    residual = residual_analysis(data, model_results, rng)
    time_rows = time_law_analysis(data, paired)
    calibration_rows = calibration_scale_sensitivity(data, model_results)
    audit = numerical_audit(
        data,
        model_free,
        model_results,
        contact_results,
        hierarchy,
    )
    report = render_report(
        data,
        model_free,
        model_results,
        contact_results,
        identifiability,
        paired,
        hierarchy,
        residual,
        time_rows,
        calibration_rows,
        audit,
    )
    (RESULTS / "REPORT.md").write_text(report, encoding="utf-8")

    input_files = [
        UPSTREAM / "pixel_force_curves.npz",
        UPSTREAM / "pixel_QC.csv",
        UPSTREAM / "map_inventory_QC.csv",
        UPSTREAM / "provenance.json",
        COLLOID_UPSTREAM / "per_curve_colloid_colloid_fits.csv",
        CALIBRATION_SUMMARY,
    ]
    provenance = {
        "analysis": "advanced model-light, spatial, contact-zero, identifiability, sphere-plane PB, joint-model, residual and time-law analysis of 31-08-26 map2--5",
        "claim_status": "apparent_model_conditioned_same_pixel_diagnostic",
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256_file(path) for path in input_files
        },
        "raw_map_hashes": {
            str(row.source): str(row.sha256) for row in data["inventory"].itertuples()
        },
        "temperature_C": raw.CONTEXT_TEMPERATURE_C,
        "cantilever": {
            "inferred_upstream": "D5",
            "spring_constant_N_per_m": data["spring_N_per_m"],
            "spring_repeatability_sd_N_per_m": data["spring_sd_N_per_m"],
            "global_liquid_InvOLS_nm_per_V": float(
                data["inventory"]["global_InvOLS_used_nm_per_V"].iloc[0]
            ),
        },
        "experimental_pairing": {
            "maps": [2, 3, 4, 5],
            "relative_times_min": data["times_min"].tolist(),
            "curves_per_map": 64,
            "physical_coordinates": "exactly paired upstream row/column and source-coordinate audit",
        },
        "primary_model": {
            "name": PRIMARY_MODEL,
            "geometry": "silica sphere--silica plane",
            "electrostatics": "equal constant potential nonlinear symmetric 1:1 Poisson--Boltzmann Derjaguin",
            "radius_m": RADIUS_M,
            "epsilon_r": EPSILON_R,
            "temperature_K": TEMPERATURE_K,
            "Hamaker_J": HAMAKER_J,
            "water_viscosity_Pa_s": VISCOSITY_PA_S,
            "hydrodynamics": "fixed ideal no-slip sphere-plane lubrication 6*pi*eta*R^2*v/D",
            "free_parameters_per_independent_curve": [
                "lambda_D_nm",
                "zeta_magnitude_mV",
                "constant_baseline_pN",
            ],
            "fit_window_nm": [FIT_MIN_NM, FIT_MAX_NM],
            "distance_step_nm": 5.0,
            "objective": "ordinary least squares with analytical per-curve baseline; same objective for model comparison",
            "parameter_grid": audit["grid"],
        },
        "contact_zero": {
            "fixed_offsets_nm": CONTACT_OFFSETS_NM.tolist(),
            "interpretation": "D_true=D_nominal+offset",
            "profile": "discrete nuisance diagnostic; boundary hits are not physical D0 estimates",
        },
        "statistics": {
            "random_seed": RANDOM_SEED,
            "Moran_permutations": MORAN_PERMUTATIONS,
            "row_cluster_bootstrap_samples": BOOTSTRAP_SAMPLES,
            "row_signflip": "exact all 2^8 sign assignments",
            "multiple_testing": "Benjamini-Hochberg within explicitly named analysis families",
            "distance_bin_correlation": "AIC/BIC/profile/distance tests diagnostic only",
        },
        "numerical_audit": audit,
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "dependency_hashes": {
            str(Path(module.__file__).resolve().relative_to(ROOT)): sha256_file(
                Path(module.__file__).resolve()
            )
            for module in (raw, colloid, base)
        },
    }
    (RESULTS / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    artifacts = [
        path
        for path in RESULTS.rglob("*")
        if path.is_file() and path.name != "artifact_manifest.sha256"
    ]
    artifacts += [Path(__file__).resolve(), *input_files]
    create_manifest(artifacts, RESULTS / "artifact_manifest.sha256")
    print(f"Wrote {RESULTS}")
    print(
        f"primary fits=256, distance bins={data['distance_nm'].size}, "
        f"contact offsets={CONTACT_OFFSETS_NM.size}, core audit={audit['core_numerical_checks_pass']}"
    )


if __name__ == "__main__":
    main()

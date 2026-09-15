#!/usr/bin/env python3
"""Test velocity-distance coupling in the 12-09-26 glycerol data.

The input forces already use the speed-conditioned constant baselines and the
batch contact reconstruction.  Within each 8-map palindrome block, symmetric
same-speed maps are averaged.  The resulting distance-by-speed matrix is then
double-centred, which removes an arbitrary distance-only surface-force curve
and an arbitrary speed-only constant offset.  What remains is the interaction
that a hydrodynamic law must explain.

This is a descriptive four-block model check.  Correlated 5 nm bins are not
treated as independent experimental replicates and are not used for p-values.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.optimize import least_squares

import fit_glycerol_surface_forces as base


ROOT = Path(__file__).resolve().parents[1]
FORCE_RESULTS = ROOT / "analysis" / "glycerol_99p7_D3_force_pb_time_results"
TIME_RESULTS = ROOT / "analysis" / "glycerol_99p7_D3_time_speed_results"
DECAY_RESULTS = ROOT / "analysis" / "glycerol_99p7_D3_decay_v0_results"
OUT = ROOT / "analysis" / "glycerol_99p7_D3_hydrodynamic_scaling_results"
FIG = OUT / "figures"

SPEEDS = np.asarray((0.1, 0.3, 0.9, 2.7), dtype=np.float64)
SPEED_PAIRS = ((0.1, 0.3), (0.3, 0.9), (0.9, 2.7), (0.1, 2.7))
BLOCKS = (1, 2, 3, 4)
PRIMARY_MIN_NM = 20.0
PRIMARY_MAX_NM = 200.0
REFERENCE_DISTANCE_NM = 200.0
FIT_RANGES_NM = (
    (20.0, 200.0, "primary_20_200"),
    (20.0, 250.0, "sensitivity_20_250"),
    (20.0, 300.0, "sensitivity_20_300"),
    (30.0, 200.0, "sensitivity_30_200"),
    (50.0, 200.0, "sensitivity_50_200"),
)
MODEL_ORDER = (
    "additive_null",
    "theoretical_no_slip",
    "classical_fitted",
    "distance_shift_only",
    "velocity_power_only",
    "velocity_power_and_shift",
)
MODEL_LABELS = {
    "additive_null": "additive null",
    "theoretical_no_slip": "fixed no-slip",
    "classical_fitted": "fitted v/D",
    "distance_shift_only": "v/(D+D0)",
    "velocity_power_only": "v^alpha/D",
    "velocity_power_and_shift": "v^alpha/(D+D0)",
}
PAIR_COLORS = {
    (0.1, 0.3): "#007f5f",
    (0.3, 0.9): "#2b6cb0",
    (0.9, 2.7): "#c05621",
    (0.1, 2.7): "#5b3f8c",
}


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


def double_center(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2 or not np.all(np.isfinite(array)):
        raise ValueError("double_center requires a finite two-dimensional array")
    return (
        array
        - np.mean(array, axis=1, keepdims=True)
        - np.mean(array, axis=0, keepdims=True)
        + np.mean(array)
    )


def model_matrix(
    distances_nm: np.ndarray,
    coefficient: float,
    velocity_exponent: float,
    distance_shift_nm: float,
) -> np.ndarray:
    distances = np.asarray(distances_nm, dtype=np.float64)
    if np.any(distances + distance_shift_nm <= 0.0):
        raise ValueError("Shifted separation must stay positive")
    distance_basis = 1.0 / (distances + distance_shift_nm)
    velocity_basis = SPEEDS**velocity_exponent
    return coefficient * np.outer(
        distance_basis - np.mean(distance_basis),
        velocity_basis - np.mean(velocity_basis),
    )


def no_slip_coefficient_pN_nm_per_um_s() -> tuple[float, float]:
    viscosity_mPa_s = base.cheng_viscosity_mPa_s(0.997, base.TEMPERATURE_C)
    viscosity_Pa_s = viscosity_mPa_s * 1.0e-3
    coefficient = 6.0 * math.pi * viscosity_Pa_s * base.PROBE_RADIUS_M**2 * 1.0e15
    return viscosity_mPa_s, coefficient


def build_pair_means(
    chronology: list[dict[str, str]], force_rows: list[dict[str, str]]
) -> tuple[list[dict], dict[tuple[int, float, float], float]]:
    metadata = {int(row["acquisition_order"]): row for row in chronology}
    if len(metadata) != 32:
        raise RuntimeError("Expected a complete 32-map chronology")
    force_lookup = {
        (int(row["acquisition_order"]), float(row["distance_nm"])): float(row["map_median_force_pN"])
        for row in force_rows
    }
    distances = sorted(
        {
            float(row["distance_nm"])
            for row in force_rows
            if PRIMARY_MIN_NM <= float(row["distance_nm"]) <= 300.0
        }
    )
    if distances != list(np.arange(PRIMARY_MIN_NM, 300.0 + 0.1, 5.0)):
        raise RuntimeError("Expected complete 5 nm force bins from 20 to 300 nm")

    rows: list[dict] = []
    lookup: dict[tuple[int, float, float], float] = {}
    for block in BLOCKS:
        for speed in SPEEDS:
            maps = sorted(
                (
                    row
                    for row in chronology
                    if int(row["block"]) == block
                    and np.isclose(float(row["speed_um_per_s"]), speed)
                ),
                key=lambda row: int(row["position_in_block"]),
            )
            if len(maps) != 2:
                raise RuntimeError(f"Block {block}, speed {speed:g}: expected two maps")
            if int(maps[0]["position_in_block"]) + int(maps[1]["position_in_block"]) != 9:
                raise RuntimeError(f"Block {block}, speed {speed:g}: maps are not palindrome mates")
            early_order = int(maps[0]["acquisition_order"])
            late_order = int(maps[1]["acquisition_order"])
            for distance in distances:
                forces = np.asarray(
                    (force_lookup[(early_order, distance)], force_lookup[(late_order, distance)]),
                    dtype=np.float64,
                )
                if not np.all(np.isfinite(forces)):
                    raise RuntimeError("Non-finite map force in the analysis window")
                pair_mean = float(np.mean(forces))
                lookup[(block, distance, float(speed))] = pair_mean
                rows.append(
                    {
                        "block": block,
                        "speed_um_per_s": float(speed),
                        "distance_nm": distance,
                        "early_map_order": early_order,
                        "late_map_order": late_order,
                        "early_force_pN": float(forces[0]),
                        "late_force_pN": float(forces[1]),
                        "pair_mean_force_pN": pair_mean,
                        "pair_half_difference_pN": float(abs(forces[1] - forces[0]) / 2.0),
                    }
                )
    return rows, lookup


def pair_mean_closure(pair_rows: list[dict]) -> float:
    existing_path = DECAY_RESULTS / "blockwise_speed_pair_means.csv"
    if not existing_path.exists():
        return float("nan")
    existing = {
        (int(row["block"]), float(row["distance_nm"]), float(row["speed_um_per_s"])): float(
            row["pair_mean_force_pN"]
        )
        for row in read_csv(existing_path)
    }
    differences = []
    for row in pair_rows:
        key = (int(row["block"]), float(row["distance_nm"]), float(row["speed_um_per_s"]))
        if key in existing:
            differences.append(abs(float(row["pair_mean_force_pN"]) - existing[key]))
    if len(differences) != 4 * 4 * 4:
        raise RuntimeError("Existing v0 pair-mean closure did not cover 20/50/100/200 nm")
    return float(max(differences))


def force_matrix(
    pair_lookup: dict[tuple[int, float, float], float],
    block: int,
    minimum_nm: float,
    maximum_nm: float,
) -> tuple[np.ndarray, np.ndarray]:
    distances = np.arange(minimum_nm, maximum_nm + 0.1, 5.0, dtype=np.float64)
    matrix = np.asarray(
        [
            [pair_lookup[(block, float(distance), float(speed))] for speed in SPEEDS]
            for distance in distances
        ],
        dtype=np.float64,
    )
    if matrix.shape != (distances.size, SPEEDS.size) or not np.all(np.isfinite(matrix)):
        raise RuntimeError("Invalid force matrix")
    return distances, matrix


def analytic_coefficient(observed: np.ndarray, basis: np.ndarray) -> float:
    denominator = float(np.sum(basis**2))
    if denominator <= 0.0:
        raise RuntimeError("Degenerate interaction basis")
    return float(np.sum(observed * basis) / denominator)


def nonlinear_fit(
    observed: np.ndarray,
    distances_nm: np.ndarray,
    mode: str,
) -> tuple[float, float, float, bool, int, str]:
    if mode not in {"distance_shift_only", "velocity_power_only", "velocity_power_and_shift"}:
        raise ValueError(f"Unknown nonlinear model {mode}")

    def unpack(parameters: np.ndarray) -> tuple[float, float, float]:
        if mode == "distance_shift_only":
            coefficient, shift = parameters
            return float(coefficient), 1.0, float(shift)
        if mode == "velocity_power_only":
            coefficient, exponent = parameters
            return float(coefficient), float(exponent), 0.0
        coefficient, exponent, shift = parameters
        return float(coefficient), float(exponent), float(shift)

    def residual(parameters: np.ndarray) -> np.ndarray:
        coefficient, exponent, shift = unpack(parameters)
        return (model_matrix(distances_nm, coefficient, exponent, shift) - observed).ravel()

    if mode == "distance_shift_only":
        starts = [np.asarray((coefficient, shift)) for coefficient in (5.0e4, 1.0e5, 2.0e5) for shift in (0.0, 40.0, 120.0)]
        lower = np.asarray((0.0, -float(np.min(distances_nm)) + 0.1))
        upper = np.asarray((1.0e7, 500.0))
    elif mode == "velocity_power_only":
        starts = [np.asarray((coefficient, exponent)) for coefficient in (3.0e4, 1.0e5, 3.0e5) for exponent in (0.4, 0.7, 1.0, 1.3)]
        lower = np.asarray((0.0, 0.05))
        upper = np.asarray((1.0e7, 2.0))
    else:
        starts = [
            np.asarray((coefficient, exponent, shift))
            for coefficient in (5.0e4, 1.5e5, 4.0e5)
            for exponent in (0.4, 0.7, 1.0, 1.3)
            for shift in (0.0, 40.0, 120.0)
        ]
        lower = np.asarray((0.0, 0.05, -float(np.min(distances_nm)) + 0.1))
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
    coefficient, exponent, shift = unpack(solution.x)
    return coefficient, exponent, shift, bool(solution.success), int(solution.nfev), str(solution.message)


def fit_models(
    pair_lookup: dict[tuple[int, float, float], float], theoretical_coefficient: float
) -> tuple[list[dict], list[dict]]:
    fit_rows: list[dict] = []
    centred_rows: list[dict] = []
    for minimum, maximum, range_label in FIT_RANGES_NM:
        for block in BLOCKS:
            distances, raw = force_matrix(pair_lookup, block, minimum, maximum)
            observed = double_center(raw)
            interaction_rms = float(np.sqrt(np.mean(observed**2)))
            singular_values = np.linalg.svd(observed, compute_uv=False)
            first_fraction = float(singular_values[0] ** 2 / np.sum(singular_values**2))
            if range_label == "primary_20_200":
                for row_index, distance in enumerate(distances):
                    for column, speed in enumerate(SPEEDS):
                        centred_rows.append(
                            {
                                "block": block,
                                "distance_nm": float(distance),
                                "speed_um_per_s": float(speed),
                                "pair_mean_force_pN": float(raw[row_index, column]),
                                "double_centered_interaction_pN": float(observed[row_index, column]),
                            }
                        )

            candidates: list[tuple[str, float, float, float, bool, int, str]] = [
                ("additive_null", 0.0, 1.0, 0.0, True, 0, "closed form"),
                ("theoretical_no_slip", theoretical_coefficient, 1.0, 0.0, True, 0, "fixed physical coefficient"),
            ]
            classical_basis = model_matrix(distances, 1.0, 1.0, 0.0)
            classical_coefficient = analytic_coefficient(observed, classical_basis)
            candidates.append(("classical_fitted", classical_coefficient, 1.0, 0.0, True, 0, "closed form"))
            for model in ("distance_shift_only", "velocity_power_only", "velocity_power_and_shift"):
                coefficient, exponent, shift, success, nfev, message = nonlinear_fit(observed, distances, model)
                candidates.append((model, coefficient, exponent, shift, success, nfev, message))

            for model, coefficient, exponent, shift, success, nfev, message in candidates:
                predicted = model_matrix(distances, coefficient, exponent, shift)
                residual = observed - predicted
                rmse = float(np.sqrt(np.mean(residual**2)))
                fit_rows.append(
                    {
                        "fit_range": range_label,
                        "minimum_distance_nm": minimum,
                        "maximum_distance_nm": maximum,
                        "block": block,
                        "model": model,
                        "coefficient_pN_nm_per_um_s_to_alpha": coefficient,
                        "velocity_exponent_alpha": exponent,
                        "distance_shift_D0_nm": shift,
                        "interaction_rms_pN": interaction_rms,
                        "residual_rmse_pN": rmse,
                        "normalized_rmse": rmse / interaction_rms,
                        "first_singular_value_variance_fraction": first_fraction,
                        "coefficient_over_bulk_no_slip": coefficient / theoretical_coefficient if np.isclose(exponent, 1.0) else float("nan"),
                        "optimizer_success": success,
                        "optimizer_nfev": nfev,
                        "optimizer_message": message,
                    }
                )
    return fit_rows, centred_rows


def block_reproducibility(centred_rows: list[dict]) -> list[dict]:
    vectors: dict[int, np.ndarray] = {}
    for block in BLOCKS:
        selected = sorted(
            (row for row in centred_rows if int(row["block"]) == block),
            key=lambda row: (float(row["distance_nm"]), float(row["speed_um_per_s"])),
        )
        vector = np.asarray(
            [float(row["double_centered_interaction_pN"]) for row in selected],
            dtype=np.float64,
        )
        if vector.size != 37 * 4 or not np.all(np.isfinite(vector)):
            raise RuntimeError(f"Block {block}: incomplete primary interaction vector")
        vectors[block] = vector

    rows: list[dict] = []
    for first in BLOCKS:
        for second in BLOCKS:
            if second <= first:
                continue
            first_vector = vectors[first]
            second_vector = vectors[second]
            correlation = float(np.corrcoef(first_vector, second_vector)[0, 1])
            scale = float(np.dot(second_vector, first_vector) / np.dot(first_vector, first_vector))
            residual = second_vector - scale * first_vector
            rows.append(
                {
                    "first_block": first,
                    "second_block": second,
                    "pearson_shape_correlation": correlation,
                    "least_squares_scale_second_over_first": scale,
                    "scaled_shape_normalized_rmse": float(
                        np.sqrt(np.mean(residual**2)) / np.sqrt(np.mean(second_vector**2))
                    ),
                    "first_interaction_rms_pN": float(np.sqrt(np.mean(first_vector**2))),
                    "second_interaction_rms_pN": float(np.sqrt(np.mean(second_vector**2))),
                }
            )
    return rows


def distancewise_speed_slopes(
    pair_lookup: dict[tuple[int, float, float], float]
) -> list[dict]:
    design = np.column_stack((np.ones(SPEEDS.size), SPEEDS))
    rows: list[dict] = []
    for block in BLOCKS:
        for distance in np.arange(PRIMARY_MIN_NM, 300.0 + 0.1, 5.0):
            values = np.asarray(
                [pair_lookup[(block, float(distance), float(speed))] for speed in SPEEDS]
            )
            intercept, slope = np.linalg.lstsq(design, values, rcond=None)[0]
            residual = values - design @ np.asarray((intercept, slope))
            total = float(np.sum((values - np.mean(values)) ** 2))
            r2 = float("nan") if total <= 0.0 else 1.0 - float(np.sum(residual**2)) / total
            rows.append(
                {
                    "block": block,
                    "distance_nm": float(distance),
                    "speed_intercept_pN": float(intercept),
                    "speed_slope_pN_per_um_s": float(slope),
                    "distance_times_speed_slope_pN_nm_per_um_s": float(distance * slope),
                    "speed_linear_r2": r2,
                }
            )
    return rows


def velocity_secants(
    pair_lookup: dict[tuple[int, float, float], float]
) -> list[dict]:
    rows: list[dict] = []
    for block in BLOCKS:
        for speed_low, speed_high in SPEED_PAIRS:
            delta_speed = speed_high - speed_low
            reference_contrast = (
                pair_lookup[(block, REFERENCE_DISTANCE_NM, speed_high)]
                - pair_lookup[(block, REFERENCE_DISTANCE_NM, speed_low)]
            )
            for distance in np.arange(PRIMARY_MIN_NM, REFERENCE_DISTANCE_NM + 0.1, 5.0):
                contrast = (
                    pair_lookup[(block, float(distance), speed_high)]
                    - pair_lookup[(block, float(distance), speed_low)]
                )
                double_difference = contrast - reference_contrast
                inverse_distance_difference = 1.0 / distance - 1.0 / REFERENCE_DISTANCE_NM
                rows.append(
                    {
                        "block": block,
                        "speed_low_um_per_s": speed_low,
                        "speed_high_um_per_s": speed_high,
                        "delta_speed_um_per_s": delta_speed,
                        "distance_nm": float(distance),
                        "reference_distance_nm": REFERENCE_DISTANCE_NM,
                        "speed_contrast_pN": contrast,
                        "reference_speed_contrast_pN": reference_contrast,
                        "distance_double_difference_pN": double_difference,
                        "double_difference_per_delta_speed_pN_per_um_s": double_difference / delta_speed,
                        "classical_pair_coefficient_pN_nm_per_um_s": (
                            double_difference / delta_speed / inverse_distance_difference
                            if not np.isclose(distance, REFERENCE_DISTANCE_NM)
                            else float("nan")
                        ),
                    }
                )
    return rows


def synthetic_self_checks() -> dict[str, float]:
    distances = np.arange(20.0, 200.0 + 0.1, 5.0)
    surface = 600.0 * np.exp(-distances / 60.0)
    speed_offset = 17.0 + 31.0 * SPEEDS**2
    coefficient = 12345.0
    synthetic = surface[:, None] + speed_offset[None, :] + coefficient * np.outer(1.0 / distances, SPEEDS)
    observed = double_center(synthetic)
    expected = model_matrix(distances, coefficient, 1.0, 0.0)
    recovered = analytic_coefficient(observed, model_matrix(distances, 1.0, 1.0, 0.0))
    return {
        "synthetic_double_center_max_abs_error_pN": float(np.max(np.abs(observed - expected))),
        "synthetic_classical_coefficient_relative_error": abs(recovered - coefficient) / coefficient,
        "double_center_row_mean_max_abs_pN": float(np.max(np.abs(np.mean(observed, axis=1)))),
        "double_center_column_mean_max_abs_pN": float(np.max(np.abs(np.mean(observed, axis=0)))),
    }


def percentile_text(values: list[float], digits: int = 3) -> str:
    array = np.asarray(values, dtype=np.float64)
    return f"{np.median(array):.{digits}f} [{np.min(array):.{digits}f}, {np.max(array):.{digits}f}]"


def diagnostics_figure(fit_rows: list[dict], secant_rows: list[dict]) -> None:
    primary = [row for row in fit_rows if row["fit_range"] == "primary_20_200"]
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2))

    axis = axes[0]
    for speed_pair in SPEED_PAIRS[:3]:
        low, high = speed_pair
        distances = sorted(
            {
                float(row["distance_nm"])
                for row in secant_rows
                if np.isclose(float(row["speed_low_um_per_s"]), low)
                and np.isclose(float(row["speed_high_um_per_s"]), high)
            }
        )
        arrays = []
        for block in BLOCKS:
            lookup = {
                float(row["distance_nm"]): float(row["double_difference_per_delta_speed_pN_per_um_s"])
                for row in secant_rows
                if int(row["block"]) == block
                and np.isclose(float(row["speed_low_um_per_s"]), low)
                and np.isclose(float(row["speed_high_um_per_s"]), high)
            }
            arrays.append([lookup[distance] for distance in distances])
        values = np.asarray(arrays)
        median = np.median(values, axis=0)
        axis.plot(distances, median, lw=2.0, color=PAIR_COLORS[speed_pair], label=f"{low:g}-{high:g} um/s")
        axis.fill_between(distances, np.min(values, axis=0), np.max(values, axis=0), color=PAIR_COLORS[speed_pair], alpha=0.13)
    axis.axhline(0.0, color="0.35", lw=0.8)
    axis.set_xlabel("Separation D (nm)")
    axis.set_ylabel("Double difference / delta v (pN per um/s)")
    axis.set_title("Velocity intervals do not collapse")
    axis.legend(frameon=False)
    axis.grid(alpha=0.2)

    axis = axes[1]
    shown_models = (
        "theoretical_no_slip",
        "classical_fitted",
        "distance_shift_only",
        "velocity_power_only",
        "velocity_power_and_shift",
    )
    x = np.arange(len(shown_models))
    for block in BLOCKS:
        values = [
            float(next(row for row in primary if int(row["block"]) == block and row["model"] == model)["normalized_rmse"])
            for model in shown_models
        ]
        axis.plot(x, values, marker="o", lw=1.1, alpha=0.72, label=f"block {block}")
    axis.axhline(1.0, color="0.35", lw=0.8, ls="--", label="additive null")
    axis.set_yscale("log")
    axis.set_xticks(x, [MODEL_LABELS[model] for model in shown_models], rotation=24, ha="right")
    axis.set_ylabel("RMSE / RMS interaction")
    axis.set_title("20-200 nm interaction model error")
    axis.grid(axis="y", alpha=0.2, which="both")
    axis.legend(frameon=False, ncol=2, fontsize=8)

    fig.suptitle("12 Sep 2026 - 99.7 wt% glycerol - hydrodynamic scaling diagnostics", fontsize=14)
    fig.text(
        0.5,
        0.01,
        "Same-speed palindrome pair means; 200 nm reference on left. Bands are four-block ranges, not confidence intervals.",
        ha="center",
        fontsize=9,
        color="0.3",
    )
    fig.tight_layout(rect=(0.02, 0.06, 0.99, 0.93))
    for suffix in ("png", "svg"):
        fig.savefig(FIG / f"hydrodynamic_scaling_diagnostics.{suffix}", dpi=220)
    plt.close(fig)


def write_report(
    fit_rows: list[dict],
    secant_rows: list[dict],
    reproducibility_rows: list[dict],
    viscosity_mPa_s: float,
    theoretical_coefficient: float,
    closure_error: float,
    self_checks: dict[str, float],
) -> None:
    primary = [row for row in fit_rows if row["fit_range"] == "primary_20_200"]
    block_table = []
    for block in BLOCKS:
        classical = next(
            row for row in primary if int(row["block"]) == block and row["model"] == "classical_fitted"
        )
        flexible = next(
            row
            for row in primary
            if int(row["block"]) == block and row["model"] == "velocity_power_and_shift"
        )
        block_table.append(
            "| {block} | {svd:.5f} | {coef:.0f} | {ratio:.3f} | {classic:.3f} | {alpha:.3f} | {shift:.1f} | {flex:.3f} |".format(
                block=block,
                svd=float(classical["first_singular_value_variance_fraction"]),
                coef=float(classical["coefficient_pN_nm_per_um_s_to_alpha"]),
                ratio=float(classical["coefficient_over_bulk_no_slip"]),
                classic=float(classical["normalized_rmse"]),
                alpha=float(flexible["velocity_exponent_alpha"]),
                shift=float(flexible["distance_shift_D0_nm"]),
                flex=float(flexible["normalized_rmse"]),
            )
        )

    secant_summary = []
    for speed_low, speed_high in SPEED_PAIRS[:3]:
        selected = [
            row
            for row in secant_rows
            if np.isclose(float(row["distance_nm"]), 20.0)
            and np.isclose(float(row["speed_low_um_per_s"]), speed_low)
            and np.isclose(float(row["speed_high_um_per_s"]), speed_high)
        ]
        coefficients = [float(row["classical_pair_coefficient_pN_nm_per_um_s"]) for row in selected]
        secant_summary.append(
            "| {low:g}-{high:g} | {median:.0f} | {minimum:.0f}-{maximum:.0f} | {fraction:.3f} |".format(
                low=speed_low,
                high=speed_high,
                median=float(np.median(coefficients)),
                minimum=float(np.min(coefficients)),
                maximum=float(np.max(coefficients)),
                fraction=float(np.median(coefficients) / theoretical_coefficient),
            )
        )

    flexible_primary = [row for row in primary if row["model"] == "velocity_power_and_shift"]
    classical_primary = [row for row in primary if row["model"] == "classical_fitted"]
    shifted_primary = [row for row in primary if row["model"] == "distance_shift_only"]
    theoretical_primary = [row for row in primary if row["model"] == "theoretical_no_slip"]
    sensitivity = [
        row
        for row in fit_rows
        if row["model"] == "velocity_power_and_shift" and row["fit_range"] in {"sensitivity_20_250", "sensitivity_20_300"}
    ]

    lines = [
        "# 12-09-26 - 99.7 wt% glycerol - hydrodynamic scaling test",
        "",
        "## 结论",
        "",
        "**数据验证了 baseline 校正后仍存在很强的 velocity-distance coupling，因此纯粹的 `F_hyd=F(v)` 不能描述残余信号；但数据不验证把 nominal scanner speed 直接代入的严格 no-slip `F_hyd=6 pi eta R^2 v/D`。**",
        "",
        "在每个 8-map block 内，先平均同速度的 palindrome 对称 map，再把 distance-by-speed force matrix 双重中心化。这样会精确消去任意的 distance-only 项 `A(D)` 和任意的 speed-only 常数 `B(v)`；如果 hydrodynamics 只剩 `B(v)`，中心化后的 interaction 应为零。实际 interaction 的第一 singular component 在四个 block 中解释 99.90%-99.95% 的平方幅值，说明残余不是随机杂乱项，而是高度可复现、近似 separable 的 `v x D` 耦合。",
        f"四个 block 的 interaction 图样两两 Pearson correlation 为 {min(float(row['pearson_shape_correlation']) for row in reproducibility_rows):.4f}-{max(float(row['pearson_shape_correlation']) for row in reproducibility_rows):.4f}；这描述跨 block 的形状复现度，不把相关的 distance bins 当作独立样本做显著性检验。",
        "",
        "但该耦合的形状不是严格 `v/D`。20-200 nm 内，自由拟合振幅的 `K v/D` 仍留下 interaction RMS 的 25.6%-28.7%；经验式 `K v^alpha/(D+D0)` 只留下 3.3%-3.8%，四个 block 给出 `alpha=0.588-0.719`、`D0=33.5-40.6 nm`。后者只是 compact empirical descriptor，不是 shear-thinning exponent 或 slip length 的识别。",
        "",
        "有限距离 baseline 的物理含义也需要写清楚。若 baseline 对应间距窗口 `D_b`，理想 lubrication force 经过 constant subtraction 后成为",
        "",
        "```text",
        "F_corr(D,v) = F_surface(D) - <F_surface(D_b)> + K v [1/D - <1/D_b>]",
        "```",
        "",
        "所以 residual 本来就不应消失；baseline 只把 hydrodynamic zero 移到有限距离。四个速度使用不同 `D_b` 时，两个被减去的窗口平均都可产生 speed-only offset。下面的双中心化正是为了把该 offset 消掉，只检验剩余的 distance-dependent interaction。",
        "",
        "## 检验定义",
        "",
        "设同一 block 内、同速度对称 pair 的平均力为 `F_b(D,v)`。分析量为",
        "",
        "```text",
        "I_b(D,v) = F_b(D,v) - <F_b>_v - <F_b>_D + <F_b>_{D,v}",
        "```",
        "",
        "因此 `F_surface(D)`、constant baseline subtraction 留下的 speed-only offset，以及 block 常数都不进入 `I_b`。若 `F_hyd=K v/D`，则必须有",
        "",
        "```text",
        "I_b(D,v) = K [1/D - <1/D>_D] [v - <v>_v].",
        "```",
        "",
        "5 nm bins 来自同一批曲线，彼此相关；它们只定义曲线形状，不被当作独立重复。重复单位是四个 acquisition blocks，表中范围是 block range，不是 confidence interval。",
        "",
        "## 20-200 nm block 结果",
        "",
        f"bulk no-slip 参照使用仓库既有 Cheng correlation：99.7 wt%、{base.TEMPERATURE_C:.1f} C 时 `eta={viscosity_mPa_s:.3f} mPa s`，`R={base.PROBE_RADIUS_M*1e6:.6f} um`，因此 `K_bulk=6 pi eta R^2={theoretical_coefficient:.0f} pN nm/(um/s)`。",
        "",
        "| block | first SVD fraction | fitted K for v/D | K/K_bulk | v/D NRMSE | alpha | D0 (nm) | flexible NRMSE |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
        *block_table,
        "",
        f"`v/D` 的 block-median NRMSE 为 {np.median([float(row['normalized_rmse']) for row in classical_primary]):.3f}；允许一个 distance shift、但保持 velocity linear 时降至 {np.median([float(row['normalized_rmse']) for row in shifted_primary]):.3f}；同时允许 sublinear nominal-speed dependence 后降至 {np.median([float(row['normalized_rmse']) for row in flexible_primary]):.3f}。固定 bulk viscosity 的 no-slip 预测 NRMSE 为 {np.median([float(row['normalized_rmse']) for row in theoretical_primary]):.2f}，量级也不闭合。",
        "",
        "## 速度线性检验",
        "",
        "为排除 constant baseline 的影响，再取 20 nm 与 200 nm 的 distance double difference，并除以相邻速度差。若 force 对 nominal speed 线性，三个速度区间反推出的 `K` 应相同；实际随速度升高系统性下降：",
        "",
        "| speed interval (um/s) | block-median K (pN nm/(um/s)) | block range | median / K_bulk |",
        "|---:|---:|---:|---:|",
        *secant_summary,
        "",
        "这不是由一个 speed-only baseline 常数造成的，因为 distance double difference 会把该常数精确消掉。它说明至少以 nominal scanner speed 为横轴时，`F proportional to v` 不成立。",
        "",
        "## 物理解读边界",
        "",
        "1. **可以说**：残余信号包含强而可复现的 velocity-distance interaction；纯 `F(v)` 被数据排除，hydrodynamic-like distance dependence 是合理解释。",
        "2. **不能说**：本数据已验证 classical no-slip `v/D` 或测得 bulk viscosity。自由 `v/D` 振幅只有 bulk no-slip 系数的约 7%-8%，而且形状残差有系统性。",
        "3. 公式中的速度应是 instantaneous gap-closing speed `U_gap=-dD/dt`，不是自动等于 nominal scanner speed。高黏度下 cantilever deflection/relaxation 会使两者显著不同；当前本地只保留 derived CSV，raw JPK 不在 checkout，因而本包不能重建 `U_gap(D,t)`。",
        "4. `D0` 也可能混合 contact-zero bias、hydrodynamic compliance、slip、粗糙度和控制回路响应，不能直接命名为 slip length。`alpha<1` 可能主要反映 nominal-to-gap velocity mapping，不足以证明 glycerol shear thinning。",
        "5. pair averaging 只压低 block 内一阶时间漂移；四个 block 的幅值仍随 history 改变。距离 bins 的高相关性与仅四个 speeds 也限制了参数 identifiability。",
        "",
        "## 下一步判别实验",
        "",
        "决定性检验应从每条 raw curve 同时重建 `D(t)` 与 `U_gap(t)=-dD/dt`，再拟合 `F(D,t)=F_surface(D)+6 pi eta R^2 U_gap/D+C_curve`；同时保留 block/palindrome 结构，并检查 approach/retract 的 hydrodynamic sign reversal。若使用 bulk `eta` 后系数、速度线性和 `1/D` 三项同时闭合，才可称为验证 classical drainage law。",
        "",
        "## 数值与可复现性",
        "",
        f"- 重建的 20/50/100/200 nm pair means 与既有 v0 包最大差为 `{closure_error:.3e} pN`。",
        f"- synthetic additive-plus-`v/D` self-check 的 interaction 最大误差为 `{self_checks['synthetic_double_center_max_abs_error_pN']:.3e} pN`，K 相对恢复误差为 `{self_checks['synthetic_classical_coefficient_relative_error']:.3e}`。",
        f"- 将上限扩到 250/300 nm 后，经验式 alpha 的全部 block 范围为 `{min(float(row['velocity_exponent_alpha']) for row in sensitivity):.3f}-{max(float(row['velocity_exponent_alpha']) for row in sensitivity):.3f}`，D0 为 `{min(float(row['distance_shift_D0_nm']) for row in sensitivity):.1f}-{max(float(row['distance_shift_D0_nm']) for row in sensitivity):.1f} nm`；定性结论不变。",
        "- 所有 nonlinear fits 均使用多起点 bounded least-squares；CSV 保存 solver status。没有从相关 distance bins 构造 p-value 或 nominal CI。",
        "",
        "## 文件",
        "",
        "- `palindrome_pair_force_matrix.csv`: 20-300 nm block x speed pair means。",
        "- `double_centered_interaction_20_200nm.csv`: 主检验矩阵。",
        "- `block_interaction_reproducibility.csv`: 四个 block 的 pairwise shape correlation 与 scale。",
        "- `interaction_model_fits.csv`: 主范围与四个 fit-window sensitivities。",
        "- `distancewise_speed_slopes.csv`: 固定距离的 all-four-speed descriptive slopes。",
        "- `velocity_interval_double_differences.csv`: 以 200 nm 为 reference 的 baseline-invariant speed contrasts。",
        "- `figures/hydrodynamic_scaling_diagnostics.*`: 速度区间 collapse 与模型误差。",
        "- `provenance.json`, `artifact_manifest.sha256`: 输入 hash、定义、软件与产物校验。",
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
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    force_path = FORCE_RESULTS / "map_force_by_separation.csv"
    chronology_path = TIME_RESULTS / "map_time_speed_observables.csv"
    chronology = sorted(read_csv(chronology_path), key=lambda row: int(row["acquisition_order"]))
    force_rows = read_csv(force_path)

    pair_rows, pair_lookup = build_pair_means(chronology, force_rows)
    closure_error = pair_mean_closure(pair_rows)
    viscosity_mPa_s, theoretical_coefficient = no_slip_coefficient_pN_nm_per_um_s()
    fit_rows, centred_rows = fit_models(pair_lookup, theoretical_coefficient)
    reproducibility_rows = block_reproducibility(centred_rows)
    slope_rows = distancewise_speed_slopes(pair_lookup)
    secant_rows = velocity_secants(pair_lookup)
    self_checks = synthetic_self_checks()

    if closure_error > 1.0e-9:
        raise RuntimeError(f"Pair-mean closure failed: {closure_error:g} pN")
    if max(self_checks.values()) > 1.0e-9:
        raise RuntimeError(f"Synthetic self-check failed: {self_checks}")
    if not all(bool(row["optimizer_success"]) for row in fit_rows):
        raise RuntimeError("At least one interaction-model fit failed")

    write_csv(OUT / "palindrome_pair_force_matrix.csv", pair_rows)
    write_csv(OUT / "double_centered_interaction_20_200nm.csv", centred_rows)
    write_csv(OUT / "block_interaction_reproducibility.csv", reproducibility_rows)
    write_csv(OUT / "interaction_model_fits.csv", fit_rows)
    write_csv(OUT / "distancewise_speed_slopes.csv", slope_rows)
    write_csv(OUT / "velocity_interval_double_differences.csv", secant_rows)
    diagnostics_figure(fit_rows, secant_rows)
    write_report(
        fit_rows,
        secant_rows,
        reproducibility_rows,
        viscosity_mPa_s,
        theoretical_coefficient,
        closure_error,
        self_checks,
    )

    provenance = {
        "script": str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
        "inputs": {
            "force_csv": str(force_path.relative_to(ROOT)).replace("\\", "/"),
            "force_csv_sha256": sha256_file(force_path),
            "chronology_csv": str(chronology_path.relative_to(ROOT)).replace("\\", "/"),
            "chronology_csv_sha256": sha256_file(chronology_path),
            "existing_pair_means_csv": str((DECAY_RESULTS / "blockwise_speed_pair_means.csv").relative_to(ROOT)).replace("\\", "/"),
            "existing_pair_means_csv_sha256": sha256_file(DECAY_RESULTS / "blockwise_speed_pair_means.csv"),
        },
        "data_contract": {
            "statistical_replication": "four 8-map acquisition blocks",
            "within_block_unit": "arithmetic mean of the two symmetric same-speed palindrome maps",
            "distance_bins": "5 nm correlated curve coordinates; equal-weight shape samples, not independent replicates",
            "primary_fit_range_nm": [PRIMARY_MIN_NM, PRIMARY_MAX_NM],
            "speed_um_per_s": SPEEDS.tolist(),
            "speed_semantics": "nominal scanner speed, not reconstructed instantaneous gap-closing speed",
        },
        "interaction_definition": "double-centre each block distance-by-speed pair-mean force matrix to remove arbitrary A(D)+B(v)",
        "block_reproducibility": {
            "pairwise_shape_correlation_min": min(
                float(row["pearson_shape_correlation"]) for row in reproducibility_rows
            ),
            "pairwise_shape_correlation_max": max(
                float(row["pearson_shape_correlation"]) for row in reproducibility_rows
            ),
            "interpretation": "descriptive cross-block shape replication; no independent-bin p-value",
        },
        "classical_model": "I=K*(1/D-mean_D(1/D))*(v-mean_v(v))",
        "empirical_model": "I=K*(1/(D+D0)-mean_D)*(v**alpha-mean_v)",
        "bulk_reference": {
            "viscosity_model": "Cheng glycerol-water mass-fraction correlation implemented in fit_glycerol_surface_forces.py",
            "glycerol_mass_fraction": 0.997,
            "temperature_C": base.TEMPERATURE_C,
            "viscosity_mPa_s": viscosity_mPa_s,
            "probe_radius_m": base.PROBE_RADIUS_M,
            "coefficient_pN_nm_per_um_s": theoretical_coefficient,
        },
        "self_checks": {
            **self_checks,
            "existing_pair_mean_closure_max_abs_pN": closure_error,
        },
        "claim_boundary": "supports a reproducible velocity-distance interaction; does not identify hydrodynamics uniquely, validate nominal-speed no-slip v/D, infer viscosity, or infer slip length",
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
                "pair_rows": len(pair_rows),
                "fit_rows": len(fit_rows),
                "secant_rows": len(secant_rows),
                "pair_mean_closure_max_abs_pN": closure_error,
                "synthetic_checks": self_checks,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

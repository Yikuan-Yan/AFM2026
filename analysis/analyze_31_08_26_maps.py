#!/usr/bin/env python3
"""Raw reconstruction and QC for the five 31-08-26 JPK force maps.

The map files do not explicitly name the cantilever or liquid.  Cantilever
identity is inferred from the pair of stored JPK calibration multipliers by
comparison with the D4-D6 calibration files.  Force is then rebuilt from raw
vDeflection/measuredHeight using a hard-contact InvOLS consensus from all maps
and the independently calibrated spring constant of the inferred cantilever.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import sys
from pathlib import Path
from zipfile import ZipFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import (
    linregress,
    mannwhitneyu,
    pearsonr,
    spearmanr,
    theilslopes,
    ttest_1samp,
    ttest_ind,
    ttest_rel,
    wilcoxon,
)


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import analyze_27_08_26_palindrome_pilot as pilot  # noqa: E402
import analyze_velocity_systematics as events  # noqa: E402
import fit_glycerol_surface_forces as base  # noqa: E402


DATA_ROOT = ROOT / "31-08-26"
RESULTS = ROOT / "analysis" / "maps_31_08_26_results"
FIGURES = RESULTS / "figures"
CALIBRATION_SUMMARY = (
    ROOT / "analysis" / "palindrome_27_08_26_full_results" / "calibration_summary.csv"
)
CALIBRATION_RAW_ROOT = ROOT / "27-08-26" / "calibration"

ARCHIVE_NAME = "documents-export-2026-8-31.zip"
ARCHIVE_SHA256_BEFORE_DELETION = (
    "dc299da6d0104fa0047602eac1c9068b8058bc88bd4cfc0c59c21fc07eb50d47"
)
CONTACT_SPAN_NM = 40.0
CONTACT_CHECK_SPANS_NM = (30.0, 50.0)
CONTACT_R2_MIN = 0.995
CONTACT_SENSITIVITY_RANGE_NM_PER_V = (40.0, 120.0)
BIN_CENTERS_NM = np.arange(5.0, 505.0, 5.0)
TARGET_DISTANCES_NM = (20.0, 50.0, 100.0, 200.0)
CONTEXT_TEMPERATURE_C = 25.6
CONTEXT_PROBE_RADIUS_M = 4.546848945303745e-6


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


def finite_median(rows: list[dict], field: str) -> float:
    values = np.asarray([row[field] for row in rows], dtype=np.float64)
    values = values[np.isfinite(values)]
    return float(np.median(values)) if values.size else float("nan")


def relative_source(path: Path) -> str:
    return str(path.relative_to(ROOT))


def stored_calibration(path: Path) -> tuple[float, float]:
    with ZipFile(path) as archive:
        properties = base.parse_properties(archive.read("shared-data/header.properties"))
    sensitivity = float(
        properties[
            "lcd-info.1.conversion-set.conversion.distance.scaling.multiplier"
        ]
    )
    spring = float(
        properties["lcd-info.1.conversion-set.conversion.force.scaling.multiplier"]
    )
    return sensitivity, spring


def infer_cantilever(
    sources: list[base.SourceData],
) -> tuple[str, float, float, dict[str, dict[str, float | str]]]:
    stored_sensitivities = np.asarray(
        [source.stored_sensitivity_m_per_V for source in sources], dtype=np.float64
    )
    stored_springs = np.asarray(
        [source.stored_spring_constant_N_per_m for source in sources], dtype=np.float64
    )
    if not np.all(np.isfinite(stored_sensitivities)) or not np.all(
        np.isfinite(stored_springs)
    ):
        raise RuntimeError("new map stored-calibration metadata are not finite")
    if np.ptp(stored_sensitivities) > 1e-15 or np.ptp(stored_springs) > 1e-9:
        raise RuntimeError("the new maps do not share one stored calibration")

    new_sensitivity = float(np.median(stored_sensitivities))
    new_spring = float(np.median(stored_springs))
    fingerprints: dict[str, dict[str, float | str]] = {}
    for name in ("D4", "D5", "D6"):
        candidates = sorted((CALIBRATION_RAW_ROOT / name).glob("*.jpk-force"))
        if not candidates:
            raise FileNotFoundError(f"missing calibration fingerprint for {name}")
        sensitivity, spring = stored_calibration(candidates[0])
        log_distance = math.hypot(
            math.log(new_sensitivity / sensitivity), math.log(new_spring / spring)
        )
        fingerprints[name] = {
            "reference_file": relative_source(candidates[0]),
            "stored_sensitivity_nm_per_V": sensitivity * 1e9,
            "stored_spring_constant_N_per_m": spring,
            "log_ratio_distance_to_new_maps": log_distance,
        }
    inferred = min(
        fingerprints,
        key=lambda name: float(fingerprints[name]["log_ratio_distance_to_new_maps"]),
    )

    calibrated = [
        row for row in read_csv(CALIBRATION_SUMMARY) if row["cantilever"] == inferred
    ]
    if len(calibrated) != 1:
        raise RuntimeError(f"missing independent calibration summary for {inferred}")
    spring = float(calibrated[0]["spring_constant_N_per_m"])
    spring_sd = float(calibrated[0]["spring_constant_repeatability_sd_N_per_m"])
    return inferred, spring, spring_sd, fingerprints


def physical_coordinates_um(path: Path, point_count: int) -> np.ndarray:
    x_values: list[float] = []
    y_values: list[float] = []
    x_key = (
        "force-segment-header.environment.xy-scanner-position-map.xy-scanner."
        "tip-scanner.position.x"
    )
    y_key = (
        "force-segment-header.environment.xy-scanner-position-map.xy-scanner."
        "tip-scanner.position.y"
    )
    with ZipFile(path) as archive:
        for point in range(point_count):
            properties = base.parse_properties(
                archive.read(f"index/{point}/segments/0/segment-header.properties")
            )
            x_values.append(float(base._required(properties, x_key)))
            y_values.append(float(base._required(properties, y_key)))
    coordinates = np.column_stack(
        [
            np.asarray(x_values, dtype=np.float64) * 1e6,
            np.asarray(y_values, dtype=np.float64) * 1e6,
        ]
    )
    if not np.all(np.isfinite(coordinates)):
        raise RuntimeError(f"non-finite XY coordinates in {path}")
    return coordinates


def physical_bounds_um(path: Path, point_count: int) -> dict[str, float]:
    coordinates = physical_coordinates_um(path, point_count)
    x = coordinates[:, 0]
    y = coordinates[:, 1]
    return {
        "x_min_um": float(np.min(x)),
        "x_max_um": float(np.max(x)),
        "y_min_um": float(np.min(y)),
        "y_max_um": float(np.max(y)),
        "x_center_um": float((np.min(x) + np.max(x)) / 2.0),
        "y_center_um": float((np.min(y) + np.max(y)) / 2.0),
    }


def contact_analysis(
    sources: list[base.SourceData],
) -> tuple[list[dict], float]:
    rows: list[dict] = []
    retained_global: list[float] = []
    for source in sources:
        map_rows: list[dict] = []
        valid_values: list[float] = []
        for curve in source.curves:
            curve.far_field_fit = base.fit_far_field_drift(
                curve.measured_height_m, curve.deflection_V
            )
            corrected = curve.deflection_V - base.baseline_voltage(
                curve.measured_height_m, curve.deflection_V, curve.far_field_fit
            )
            primary = pilot.terminal_contact_fit(
                curve.measured_height_m, corrected, CONTACT_SPAN_NM
            )
            checks = {
                span: pilot.terminal_contact_fit(
                    curve.measured_height_m, corrected, span
                )
                for span in CONTACT_CHECK_SPANS_NM
            }
            curve.contact_fit = base.ContactFit(
                sensitivity_m_per_V=float(primary["sensitivity_m_per_V"]),
                start=int(primary["start"]),
                stop=int(primary["stop"]),
                slope_V_per_m=float(primary["slope_V_per_m"]),
                intercept_V=float(primary["intercept_V"]),
                r2=float(primary["r2"]),
            )
            sensitivity = float(primary["sensitivity_m_per_V"])
            valid = bool(
                float(primary["r2"]) >= CONTACT_R2_MIN
                and CONTACT_SENSITIVITY_RANGE_NM_PER_V[0] * 1e-9
                <= sensitivity
                <= CONTACT_SENSITIVITY_RANGE_NM_PER_V[1] * 1e-9
            )
            if valid:
                valid_values.append(sensitivity)
            row, column = base.map_pixel_from_index(
                int(curve.point_index), int(source.map_grid_i), source.map_back_and_forth
            )
            record = {
                "source": relative_source(source.path),
                "timestamp": source.timestamp.isoformat(),
                "point_index": int(curve.point_index),
                "row": row,
                "column": column,
                "contact_invOLS_40nm_nm_per_V": sensitivity * 1e9,
                "contact_invOLS_30nm_nm_per_V": float(
                    checks[30.0]["sensitivity_m_per_V"]
                )
                * 1e9,
                "contact_invOLS_50nm_nm_per_V": float(
                    checks[50.0]["sensitivity_m_per_V"]
                )
                * 1e9,
                "contact_fit_R2": float(primary["r2"]),
                "contact_points": int(primary["points"]),
                "valid_before_map_outlier_filter": valid,
                "retained_for_global_consensus": False,
            }
            rows.append(record)
            map_rows.append(record)

        values = np.asarray(valid_values, dtype=np.float64)
        if values.size < 56:
            raise RuntimeError(
                f"{source.path}: only {values.size} valid hard-contact fits"
            )
        center = float(np.median(values))
        mad = pilot.robust_mad(values)
        tolerance = max(3.5 * mad, 0.5e-9)
        retained = values[np.abs(values - center) <= tolerance]
        if retained.size < 56:
            raise RuntimeError(
                f"{source.path}: only {retained.size} contact fits after outlier filter"
            )
        for record in map_rows:
            value = float(record["contact_invOLS_40nm_nm_per_V"]) * 1e-9
            record["retained_for_global_consensus"] = bool(
                record["valid_before_map_outlier_filter"]
                and abs(value - center) <= tolerance
            )
        source.sensitivity_anchor_m_per_V = float(np.median(retained))
        source.sensitivity_anchor_mad_m_per_V = pilot.robust_mad(retained)
        source.sensitivity_valid_curves = int(retained.size)
        retained_global.extend(retained.tolist())

    global_sensitivity = float(np.median(np.asarray(retained_global)))
    if not np.isfinite(global_sensitivity):
        raise RuntimeError("global hard-contact sensitivity is non-finite")
    for source in sources:
        source.sensitivity_used_m_per_V = global_sensitivity
        source.sensitivity_method = (
            "global median of retained 40 nm hard-contact fits from both 31-08-26 maps"
        )
    return rows, global_sensitivity


def analyze_maps(
    sources: list[base.SourceData], spring_N_per_m: float
) -> tuple[list[dict], list[dict], list[dict], list[dict], dict[str, np.ndarray]]:
    base.SPRING_CONSTANT_N_PER_M = spring_N_per_m
    pixel_rows: list[dict] = []
    map_rows: list[dict] = []
    force_rows: list[dict] = []
    trend_rows: list[dict] = []
    matrices: dict[str, np.ndarray] = {}

    for map_order, source in enumerate(sources, start=1):
        key = relative_source(source.path)
        point_count = int(source.map_grid_i) * int(source.map_grid_j)
        if point_count != 64 or len(source.curves) != 64:
            raise RuntimeError(f"{source.path}: expected one complete 8x8 map")
        retracts, retract_skipped = events.load_branch(source.path, "retract")
        retract_by_point = {
            int(curve.point_index): curve
            for curve in retracts
            if curve.point_index is not None
        }
        line_matrix = np.full((64, BIN_CENTERS_NM.size), np.nan, dtype=np.float64)
        constant_matrix = np.full_like(line_matrix, np.nan)
        scanner_speeds: list[float] = []
        gap_speeds: list[float] = []
        terminal_loads: list[float] = []
        far_slopes: list[float] = []
        contact_heights: list[float] = []
        event_rows: list[dict] = []

        for acquisition_order, curve in enumerate(source.curves, start=1):
            if curve.point_index is None or curve.contact_fit is None:
                raise RuntimeError(f"{source.path}: incomplete prepared curve")
            point = int(curve.point_index)
            if point < 0 or point >= 64:
                raise RuntimeError(f"{source.path}: invalid point index {point}")
            sensitivity = float(source.sensitivity_used_m_per_V)
            far = curve.far_field_fit
            if far is None:
                raise RuntimeError(f"{source.path}: missing far-field fit")
            baseline = base.baseline_voltage(
                curve.measured_height_m, curve.deflection_V, far
            )
            corrected_V = curve.deflection_V - baseline
            delta_m = sensitivity * corrected_V
            fit = curve.contact_fit
            slope, intercept, _, _ = base.robust_line(
                curve.measured_height_m[fit.start : fit.stop],
                corrected_V[fit.start : fit.stop],
            )
            contact_height_m = float(-intercept / slope)
            distance_nm = (
                curve.measured_height_m + delta_m - contact_height_m
            ) * 1e9
            force_line_pN = spring_N_per_m * delta_m * 1e12
            reference_V = float(np.median(curve.deflection_V[: far.n_points]))
            force_constant_pN = (
                spring_N_per_m
                * sensitivity
                * (curve.deflection_V - reference_V)
                * 1e12
            )
            precontact = np.arange(distance_nm.size) < fit.start
            line_matrix[point] = pilot.bin_median(
                distance_nm, force_line_pN, precontact
            )
            constant_matrix[point] = pilot.bin_median(
                distance_nm, force_constant_pN, precontact
            )

            time_s = (
                np.arange(distance_nm.size, dtype=np.float64) + 0.5
            ) * curve.duration_s / distance_nm.size
            trim = max(5, int(math.ceil(0.10 * distance_nm.size)))
            scanner_slope, _, _, _ = base.robust_line(
                time_s[trim:-trim], curve.measured_height_m[trim:-trim]
            )
            scanner_speed = abs(scanner_slope) * 1e6
            gap_window = precontact & (distance_nm >= 20.0) & (distance_nm <= 200.0)
            if np.count_nonzero(gap_window) >= 40:
                gap_slope, _, _, _ = base.robust_line(
                    time_s[gap_window], distance_nm[gap_window] * 1e-9
                )
                gap_speed = abs(gap_slope) * 1e6
            else:
                gap_speed = float("nan")
            far_slope = (
                spring_N_per_m * sensitivity * far.slope_V_per_m * 1e5
            )
            terminal_load_nN = float(np.median(force_line_pN[-10:]) / 1e3)
            event = events.analyze_pair(
                curve, retract_by_point.get(point), sensitivity
            )
            row, column = base.map_pixel_from_index(
                point, int(source.map_grid_i), source.map_back_and_forth
            )
            scanner_speeds.append(scanner_speed)
            gap_speeds.append(gap_speed)
            terminal_loads.append(terminal_load_nN)
            far_slopes.append(far_slope)
            contact_heights.append(contact_height_m * 1e6)
            event_rows.append(event)
            pixel_rows.append(
                {
                    "source": key,
                    "timestamp": source.timestamp.isoformat(),
                    "map_order": map_order,
                    "point_index": point,
                    "row": row,
                    "column": column,
                    "pixel_acquisition_order": acquisition_order,
                    "global_InvOLS_used_nm_per_V": sensitivity * 1e9,
                    "local_contact_InvOLS_nm_per_V": fit.sensitivity_m_per_V * 1e9,
                    "contact_fit_R2": fit.r2,
                    "terminal_load_nN": terminal_load_nN,
                    "far_slope_pN_per_100nm": far_slope,
                    "far_noise_pN": spring_N_per_m
                    * sensitivity
                    * far.residual_mad_V
                    * 1e12,
                    "scanner_speed_um_per_s": scanner_speed,
                    "gap_speed_20_200nm_um_per_s": gap_speed,
                    "approach_snap_detected": event["approach_snap_detected"],
                    "approach_snap_distance_nm": event["approach_snap_distance_nm"],
                    "retract_free_baseline_valid": event[
                        "retract_free_baseline_valid"
                    ],
                    "retract_snapoff_detected": event["retract_snapoff_detected"],
                    "retract_pull_off_censored": event[
                        "retract_pull_off_censored"
                    ],
                    "retract_pull_off_force_nN": event[
                        "retract_pull_off_force_nN"
                    ],
                    "retract_detachment_piezo_travel_nm": event[
                        "retract_detachment_piezo_travel_nm"
                    ],
                }
            )

        matrices[key + "|line"] = line_matrix
        matrices[key + "|constant"] = constant_matrix
        for method, matrix in (
            ("linear_drift_corrected", line_matrix),
            ("far_constant_referenced", constant_matrix),
        ):
            for distance_index, distance_nm in enumerate(BIN_CENTERS_NM):
                values = matrix[:, distance_index]
                values = values[np.isfinite(values)]
                if not values.size:
                    continue
                force_rows.append(
                    {
                        "source": key,
                        "timestamp": source.timestamp.isoformat(),
                        "map_order": map_order,
                        "baseline_method": method,
                        "distance_nm": distance_nm,
                        "force_median_pN": float(np.median(values)),
                        "force_q25_pN": float(np.quantile(values, 0.25)),
                        "force_q75_pN": float(np.quantile(values, 0.75)),
                        "valid_pixels": int(values.size),
                    }
                )

        retract_duration = float(
            np.median([curve.duration_s for curve in retracts])
        )
        approach_duration = float(np.median([curve.duration_s for curve in source.curves]))
        protocol_duration_s = 64.0 * (approach_duration + retract_duration)
        bounds = physical_bounds_um(source.path, 64)
        map_row = {
            "source": key,
            "sha256": source.sha256,
            "timestamp": source.timestamp.isoformat(),
            "map_order": map_order,
            "grid_i": source.map_grid_i,
            "grid_j": source.map_grid_j,
            "field_u_um": source.map_ulength_m * 1e6,
            "field_v_um": source.map_vlength_m * 1e6,
            **bounds,
            "approach_duration_per_curve_s": approach_duration,
            "retract_duration_per_curve_s": retract_duration,
            "protocol_duration_no_XY_overhead_min": protocol_duration_s / 60.0,
            "first_to_64th_no_XY_overhead_min": protocol_duration_s
            * 63.0
            / 64.0
            / 60.0,
            "approach_curves": len(source.curves),
            "approach_skipped_curves": source.skipped_curves,
            "retract_curves": len(retracts),
            "retract_skipped_curves": retract_skipped,
            "stored_InvOLS_nm_per_V": source.stored_sensitivity_m_per_V * 1e9,
            "stored_spring_constant_N_per_m": source.stored_spring_constant_N_per_m,
            "map_contact_InvOLS_nm_per_V": source.sensitivity_anchor_m_per_V * 1e9,
            "map_contact_InvOLS_MAD_nm_per_V": source.sensitivity_anchor_mad_m_per_V
            * 1e9,
            "global_InvOLS_used_nm_per_V": source.sensitivity_used_m_per_V * 1e9,
            "calibrated_spring_constant_used_N_per_m": spring_N_per_m,
            "scanner_speed_median_um_per_s": float(np.median(scanner_speeds)),
            "gap_speed_20_200nm_median_um_per_s": float(
                np.nanmedian(gap_speeds)
            ),
            "terminal_load_median_nN": float(np.median(terminal_loads)),
            "far_slope_median_pN_per_100nm": float(np.median(far_slopes)),
            "far_slope_MAD_pN_per_100nm": pilot.robust_mad(
                np.asarray(far_slopes)
            ),
            "contact_height_median_um": float(np.median(contact_heights)),
            "approach_snap_detected_fraction": float(
                np.mean([bool(row["approach_snap_detected"]) for row in event_rows])
            ),
            "approach_snap_distance_median_nm": finite_median(
                event_rows, "approach_snap_distance_nm"
            ),
            "retract_free_baseline_valid_fraction": float(
                np.mean(
                    [bool(row["retract_free_baseline_valid"]) for row in event_rows]
                )
            ),
            "retract_snapoff_detected_fraction": float(
                np.mean([bool(row["retract_snapoff_detected"]) for row in event_rows])
            ),
            "retract_pull_off_censored_fraction": float(
                np.mean(
                    [bool(row["retract_pull_off_censored"]) for row in event_rows]
                )
            ),
            "retract_pull_off_force_median_nN": finite_median(
                event_rows, "retract_pull_off_force_nN"
            ),
            "retract_detachment_piezo_travel_median_nm": finite_median(
                event_rows, "retract_detachment_piezo_travel_nm"
            ),
        }
        map_rows.append(map_row)

        for target in TARGET_DISTANCES_NM:
            target_index = int(np.flatnonzero(np.isclose(BIN_CENTERS_NM, target))[0])
            values = line_matrix[:, target_index]
            valid = np.isfinite(values)
            if np.count_nonzero(valid) < 56:
                raise RuntimeError(f"{source.path}: too few pixels at {target} nm")
            x = np.arange(1, 65, dtype=np.float64)[valid]
            slope, _, _, _ = theilslopes(values[valid], x, alpha=0.95)
            trend_rows.append(
                {
                    "source": key,
                    "map_order": map_order,
                    "distance_nm": target,
                    "force_first_pN": values[0],
                    "force_64th_pN": values[63],
                    "endpoint_delta_64_minus_1_pN": values[63] - values[0],
                    "first_row_median_pN": float(np.nanmedian(values[:8])),
                    "last_row_median_pN": float(np.nanmedian(values[56:64])),
                    "last_minus_first_row_median_pN": float(
                        np.nanmedian(values[56:64]) - np.nanmedian(values[:8])
                    ),
                    "Theil_Sen_endpoint_change_pN": float(slope * 63.0),
                    "first_to_64th_no_XY_overhead_min": map_row[
                        "first_to_64th_no_XY_overhead_min"
                    ],
                }
            )

    return pixel_rows, map_rows, force_rows, trend_rows, matrices


def plot_force_distance(
    sources: list[base.SourceData], map_rows: list[dict], matrices: dict[str, np.ndarray]
) -> None:
    fig, axes = plt.subplots(1, len(sources), figsize=(22.0, 5.5), sharex=True)
    for ax, source, summary in zip(axes, sources, map_rows, strict=True):
        matrix = matrices[relative_source(source.path) + "|line"]
        for curve in matrix:
            ax.plot(BIN_CENTERS_NM, curve, color="#277da1", alpha=0.10, lw=0.8)
        median = np.nanmedian(matrix, axis=0)
        q25 = np.nanquantile(matrix, 0.25, axis=0)
        q75 = np.nanquantile(matrix, 0.75, axis=0)
        ax.fill_between(BIN_CENTERS_NM, q25, q75, color="#43aa8b", alpha=0.28)
        ax.plot(BIN_CENTERS_NM, median, color="#1d3557", lw=2.2)
        ax.axhline(0.0, color="0.35", lw=0.8)
        ax.set_xlim(10, 250)
        ax.grid(alpha=0.2)
        ax.set_title(
            f"Map {summary['map_order']}: {summary['field_u_um']:.1f}×"
            f"{summary['field_v_um']:.1f} µm\n"
            f"gap speed {summary['gap_speed_20_200nm_median_um_per_s']:.3f} µm/s"
        )
        ax.set_xlabel("Separation D (nm)")
    axes[0].set_ylabel("Force (pN), line-corrected")
    fig.suptitle("31-08-26 approach force-distance distributions")
    fig.tight_layout()
    fig.savefig(FIGURES / "force_distance_distributions.png", dpi=220)
    plt.close(fig)


def plot_force_maps(
    sources: list[base.SourceData], map_rows: list[dict], matrices: dict[str, np.ndarray]
) -> None:
    fig, axes = plt.subplots(
        len(sources),
        4,
        figsize=(17.5, 3.15 * len(sources)),
        layout="constrained",
    )
    for column, target in enumerate(TARGET_DISTANCES_NM):
        index = int(np.flatnonzero(np.isclose(BIN_CENTERS_NM, target))[0])
        maps = []
        for source in sources:
            vector = matrices[relative_source(source.path) + "|line"][:, index]
            physical = np.full((8, 8), np.nan, dtype=float)
            for point_index, value in enumerate(vector):
                row, physical_column = base.map_pixel_from_index(
                    point_index, 8, source.map_back_and_forth
                )
                physical[row, physical_column] = value
            maps.append(physical)
        combined = np.concatenate([values[np.isfinite(values)] for values in maps])
        vmin, vmax = np.quantile(combined, [0.02, 0.98])
        if np.isclose(vmin, vmax):
            vmin -= 1.0
            vmax += 1.0
        images = []
        for row, (values, summary) in enumerate(zip(maps, map_rows, strict=True)):
            image = axes[row, column].imshow(
                values, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax
            )
            images.append(image)
            axes[row, column].set_title(
                f"Map {summary['map_order']}, D={target:g} nm"
            )
            axes[row, column].set_xlabel("physical column")
            axes[row, column].set_ylabel("physical row")
        fig.colorbar(
            images[-1],
            ax=list(axes[:, column]),
            shrink=0.78,
            pad=0.02,
            label="Force (pN)",
        )
    fig.suptitle("Line-corrected force maps; color scale shared between maps at each D")
    fig.savefig(FIGURES / "force_maps_20_50_100_200nm.png", dpi=220)
    plt.close(fig)


def plot_baseline_comparison(
    sources: list[base.SourceData], map_rows: list[dict], matrices: dict[str, np.ndarray]
) -> None:
    fig, axes = plt.subplots(
        1, len(sources), figsize=(22.0, 5.5), sharex=True, sharey=True
    )
    for ax, source, summary in zip(axes, sources, map_rows, strict=True):
        key = relative_source(source.path)
        for suffix, label, color, linestyle in (
            ("|line", "far-linear", "#1d3557", "-"),
            ("|constant", "far-constant", "#e76f51", "--"),
        ):
            matrix = matrices[key + suffix]
            median = np.nanmedian(matrix, axis=0)
            q25 = np.nanquantile(matrix, 0.25, axis=0)
            q75 = np.nanquantile(matrix, 0.75, axis=0)
            ax.fill_between(BIN_CENTERS_NM, q25, q75, color=color, alpha=0.12)
            ax.plot(
                BIN_CENTERS_NM,
                median,
                color=color,
                linestyle=linestyle,
                linewidth=2.0,
                label=label,
            )
        ax.axhline(0.0, color="0.35", lw=0.8)
        ax.set_xlim(10, 250)
        ax.grid(alpha=0.2)
        ax.set_title(f"Map {summary['map_order']}")
        ax.set_xlabel("Separation D (nm)")
        ax.legend(frameon=False)
    axes[0].set_ylabel("Force (pN)")
    fig.suptitle("Sensitivity to far-field zero-force definition; median and IQR")
    fig.tight_layout()
    fig.savefig(FIGURES / "far_linear_vs_constant_baseline.png", dpi=220)
    plt.close(fig)


def plot_qc_maps(pixel_rows: list[dict], map_rows: list[dict]) -> None:
    fig, axes = plt.subplots(
        len(map_rows), 2, figsize=(9.4, 3.25 * len(map_rows))
    )
    for row_index, summary in enumerate(map_rows):
        rows = [item for item in pixel_rows if item["source"] == summary["source"]]
        sensitivity = np.full((8, 8), np.nan)
        far_slope = np.full((8, 8), np.nan)
        for item in rows:
            position = (int(item["row"]), int(item["column"]))
            sensitivity[position] = float(item["local_contact_InvOLS_nm_per_V"])
            far_slope[position] = float(item["far_slope_pN_per_100nm"])
        im0 = axes[row_index, 0].imshow(sensitivity, origin="lower", cmap="magma")
        im1 = axes[row_index, 1].imshow(far_slope, origin="lower", cmap="coolwarm")
        axes[row_index, 0].set_title(f"Map {summary['map_order']} local InvOLS")
        axes[row_index, 1].set_title(f"Map {summary['map_order']} far slope")
        fig.colorbar(im0, ax=axes[row_index, 0], label="nm/V", shrink=0.84)
        fig.colorbar(im1, ax=axes[row_index, 1], label="pN/100 nm", shrink=0.84)
        for ax in axes[row_index]:
            ax.set_xlabel("physical column")
            ax.set_ylabel("physical row")
    fig.suptitle("31-08-26 contact and far-field spatial QC")
    fig.tight_layout()
    fig.savefig(FIGURES / "contact_sensitivity_and_far_slope_maps.png", dpi=220)
    plt.close(fig)


def plot_acquisition_order(
    sources: list[base.SourceData], map_rows: list[dict], matrices: dict[str, np.ndarray]
) -> None:
    target_index = int(np.flatnonzero(np.isclose(BIN_CENTERS_NM, 50.0))[0])
    fig, axes = plt.subplots(1, len(sources), figsize=(22.0, 5.5), sharex=True)
    x = np.arange(1, 65, dtype=np.float64)
    for ax, source, summary in zip(axes, sources, map_rows, strict=True):
        values = matrices[relative_source(source.path) + "|line"][:, target_index]
        valid = np.isfinite(values)
        slope, intercept, _, _ = theilslopes(values[valid], x[valid], alpha=0.95)
        for boundary in range(8, 64, 8):
            ax.axvline(boundary + 0.5, color="0.82", lw=0.8)
        ax.plot(x, values, color="#8ecae6", lw=1.0)
        ax.scatter(x, values, color="#1261a0", s=27, zorder=2)
        ax.plot(x, slope * x + intercept, color="#e76f51", ls="--", lw=2.0)
        ax.set_title(
            f"Map {summary['map_order']}: Δtrend={slope*63:+.1f} pN, "
            f"≈{summary['first_to_64th_no_XY_overhead_min']:.2f} min"
        )
        ax.set_xlabel("Acquisition order")
        ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("Force at D=50 nm (pN)")
    fig.suptitle("Within-map force versus acquisition order (time and position confounded)")
    fig.tight_layout()
    fig.savefig(FIGURES / "force_50nm_vs_acquisition_order.png", dpi=220)
    plt.close(fig)


def plot_individual_acquisition_order_50nm(
    sources: list[base.SourceData],
    map_rows: list[dict],
    matrices: dict[str, np.ndarray],
) -> list[dict]:
    """Plot map2--4 separately with robust and ordinary linear fits."""
    target_index = int(np.flatnonzero(np.isclose(BIN_CENTERS_NM, 50.0))[0])
    x = np.arange(1, 65, dtype=np.float64)
    selected = list(zip(sources[1:], map_rows[1:], strict=True))
    vectors = [
        np.asarray(
            matrices[relative_source(source.path) + "|line"][:, target_index],
            dtype=np.float64,
        )
        for source, _ in selected
    ]
    if len(vectors) != 4 or any(vector.shape != (64,) for vector in vectors):
        raise RuntimeError("individual acquisition plots require map2--5 with 64 points")
    if any(not np.all(np.isfinite(vector)) for vector in vectors):
        raise RuntimeError("map2--5 D=50 nm vectors must be fully finite")

    pooled = np.concatenate(vectors)
    span = float(np.max(pooled) - np.min(pooled))
    padding = max(0.06 * span, 10.0)
    common_ylim = (float(np.min(pooled) - padding), float(np.max(pooled) + padding))

    fit_rows: list[dict] = []
    colors = ("#2a6fbb", "#5a3d9a", "#00897b", "#d17a22")
    for (source, summary), force_pN, point_color in zip(
        selected, vectors, colors, strict=True
    ):
        theil_slope, theil_intercept, theil_low, theil_high = theilslopes(
            force_pN, x, alpha=0.95
        )
        ordinary = linregress(x, force_pN)
        map_order = int(summary["map_order"])
        duration_min = float(summary["first_to_64th_no_XY_overhead_min"])
        fit_rows.append(
            {
                "source": relative_source(source.path),
                "map_order": map_order,
                "distance_nm": 50.0,
                "points": 64,
                "force_median_pN": float(np.median(force_pN)),
                "force_first_pN": float(force_pN[0]),
                "force_64th_pN": float(force_pN[-1]),
                "raw_64_minus_1_pN": float(force_pN[-1] - force_pN[0]),
                "Theil_Sen_slope_pN_per_sample": float(theil_slope),
                "Theil_Sen_slope_95CI_low_pN_per_sample": float(theil_low),
                "Theil_Sen_slope_95CI_high_pN_per_sample": float(theil_high),
                "Theil_Sen_intercept_pN": float(theil_intercept),
                "Theil_Sen_endpoint_change_pN": float(theil_slope * 63.0),
                "Theil_Sen_endpoint_95CI_low_pN": float(theil_low * 63.0),
                "Theil_Sen_endpoint_95CI_high_pN": float(theil_high * 63.0),
                "OLS_slope_pN_per_sample": float(ordinary.slope),
                "OLS_intercept_pN": float(ordinary.intercept),
                "OLS_endpoint_change_pN": float(ordinary.slope * 63.0),
                "OLS_r_squared": float(ordinary.rvalue**2),
                "OLS_slope_p_value": float(ordinary.pvalue),
                "first_to_64th_no_XY_overhead_min": duration_min,
            }
        )

        figure, axis = plt.subplots(figsize=(11.0, 6.3))
        for raster_row in range(8):
            if raster_row % 2 == 0:
                axis.axvspan(
                    8 * raster_row + 0.5,
                    8 * (raster_row + 1) + 0.5,
                    color="#457b9d",
                    alpha=0.045,
                    linewidth=0,
                )
        for boundary in range(8, 64, 8):
            axis.axvline(boundary + 0.5, color="0.78", linewidth=0.8)

        axis.plot(x, force_pN, color=point_color, alpha=0.48, linewidth=1.15)
        axis.scatter(
            x,
            force_pN,
            s=34,
            color=point_color,
            edgecolor="white",
            linewidth=0.45,
            zorder=3,
            label="64 force curves",
        )
        axis.plot(
            x,
            theil_slope * x + theil_intercept,
            color="#d1495b",
            linewidth=2.4,
            label="Theil–Sen robust linear fit",
        )
        axis.plot(
            x,
            ordinary.slope * x + ordinary.intercept,
            color="0.20",
            linestyle="--",
            linewidth=1.4,
            label="OLS linear fit",
        )
        axis.axhline(
            float(np.median(force_pN)),
            color="0.45",
            linestyle=":",
            linewidth=1.2,
            label="64-point median",
        )
        axis.set_xlim(0.5, 64.5)
        axis.set_ylim(*common_ylim)
        axis.set_xticks(np.arange(1, 65, 4))
        axis.set_xlabel("Acquisition order within the 8×8 map")
        axis.set_ylabel("Force at D = 50 nm (pN), far-linear corrected")
        axis.grid(axis="y", alpha=0.22)
        axis.legend(frameon=False, ncol=2, loc="upper left")
        axis.set_title(
            f"Map {map_order}: 50 nm force versus acquisition order\n"
            f"Theil–Sen slope {theil_slope:+.2f} pN/sample "
            f"(95% CI {theil_low:+.2f} to {theil_high:+.2f}); "
            f"fitted 1→64 change {theil_slope * 63.0:+.1f} pN"
        )
        figure.text(
            0.5,
            0.012,
            "Alternating bands are 8-point raster rows; acquisition order mixes elapsed time and spatial position.",
            ha="center",
            fontsize=9,
            color="0.35",
        )
        figure.tight_layout(rect=(0.0, 0.055, 1.0, 1.0))
        figure.savefig(
            FIGURES / f"map{map_order}_force_50nm_vs_acquisition_order.png",
            dpi=240,
        )
        plt.close(figure)

    return fit_rows


def physical_matrix(source: base.SourceData, vector: np.ndarray) -> np.ndarray:
    matrix = np.full((8, 8), np.nan, dtype=np.float64)
    for point_index, value in enumerate(vector):
        row, column = base.map_pixel_from_index(
            point_index, 8, source.map_back_and_forth
        )
        matrix[row, column] = value
    return matrix


def plot_paired_time_maps(
    sources: list[base.SourceData], matrices: dict[str, np.ndarray]
) -> None:
    fig, axes = plt.subplots(3, 4, figsize=(16.8, 10.4), layout="constrained")
    map2 = matrices[relative_source(sources[1].path) + "|line"]
    map3 = matrices[relative_source(sources[2].path) + "|line"]
    for column, distance_nm in enumerate(TARGET_DISTANCES_NM):
        index = int(np.flatnonzero(np.isclose(BIN_CENTERS_NM, distance_nm))[0])
        before = physical_matrix(sources[1], map2[:, index])
        after = physical_matrix(sources[2], map3[:, index])
        difference = after - before
        combined = np.concatenate([before.ravel(), after.ravel()])
        vmin, vmax = np.nanquantile(combined, [0.02, 0.98])
        limit = float(np.nanquantile(np.abs(difference), 0.98))
        limit = max(limit, 1.0)
        image_before = axes[0, column].imshow(
            before, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax
        )
        axes[1, column].imshow(
            after, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax
        )
        image_difference = axes[2, column].imshow(
            difference,
            origin="lower",
            cmap="coolwarm",
            vmin=-limit,
            vmax=limit,
        )
        axes[0, column].set_title(f"D={distance_nm:g} nm")
        fig.colorbar(
            image_before,
            ax=[axes[0, column], axes[1, column]],
            shrink=0.75,
            pad=0.02,
            label="Force (pN)",
        )
        fig.colorbar(
            image_difference,
            ax=axes[2, column],
            shrink=0.75,
            pad=0.02,
            label="Map3−Map2 (pN)",
        )
        for row in range(3):
            axes[row, column].set_xlabel("physical column")
            axes[row, column].set_ylabel("physical row")
    axes[0, 0].text(
        -0.42, 0.5, "Map 2", rotation=90, va="center", transform=axes[0, 0].transAxes
    )
    axes[1, 0].text(
        -0.42, 0.5, "Map 3", rotation=90, va="center", transform=axes[1, 0].transAxes
    )
    axes[2, 0].text(
        -0.42,
        0.5,
        "Paired difference",
        rotation=90,
        va="center",
        transform=axes[2, 0].transAxes,
    )
    fig.suptitle("Same-location time comparison: line-corrected force")
    fig.savefig(FIGURES / "map2_map3_paired_force_maps.png", dpi=220)
    plt.close(fig)


def plot_paired_time_statistics(
    sources: list[base.SourceData],
    matrices: dict[str, np.ndarray],
    paired_rows: list[dict],
) -> None:
    map2 = matrices[relative_source(sources[1].path) + "|line"]
    map3 = matrices[relative_source(sources[2].path) + "|line"]
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 9.5))
    colors = ("#264653", "#e76f51")
    for matrix, label, color in zip((map2, map3), ("Map 2", "Map 3"), colors):
        median = np.nanmedian(matrix, axis=0)
        q25 = np.nanquantile(matrix, 0.25, axis=0)
        q75 = np.nanquantile(matrix, 0.75, axis=0)
        axes[0, 0].fill_between(BIN_CENTERS_NM, q25, q75, color=color, alpha=0.14)
        axes[0, 0].plot(BIN_CENTERS_NM, median, color=color, lw=2.0, label=label)
    difference = map3 - map2
    difference_median = np.nanmedian(difference, axis=0)
    difference_q25 = np.nanquantile(difference, 0.25, axis=0)
    difference_q75 = np.nanquantile(difference, 0.75, axis=0)
    axes[0, 1].fill_between(
        BIN_CENTERS_NM, difference_q25, difference_q75, color="#457b9d", alpha=0.22
    )
    axes[0, 1].plot(BIN_CENTERS_NM, difference_median, color="#1d3557", lw=2.0)
    x = np.arange(1, 65)
    index50 = int(np.flatnonzero(np.isclose(BIN_CENTERS_NM, 50.0))[0])
    axes[1, 0].plot(x, map2[:, index50], "o-", ms=3, lw=0.8, color=colors[0], label="Map 2")
    axes[1, 0].plot(x, map3[:, index50], "o-", ms=3, lw=0.8, color=colors[1], label="Map 3")
    delta50 = difference[:, index50]
    slope, intercept, _, _ = theilslopes(delta50, x, alpha=0.95)
    axes[1, 1].scatter(x, delta50, s=24, color="#1261a0")
    axes[1, 1].plot(x, slope * x + intercept, "--", color="#e76f51", lw=2.0)
    axes[1, 1].axhline(0.0, color="0.35", lw=0.8)
    axes[0, 0].set_title("Same-location F–D distributions")
    axes[0, 1].set_title("Paired Map3−Map2 difference")
    axes[1, 0].set_title("D=50 nm by acquisition order")
    axes[1, 1].set_title(f"D=50 nm paired difference; Δtrend={slope*63:+.1f} pN")
    for ax in axes[0]:
        ax.set_xlim(10, 250)
        ax.set_xlabel("Separation D (nm)")
        ax.set_ylabel("Force (pN)")
        ax.axhline(0.0, color="0.35", lw=0.8)
        ax.grid(alpha=0.2)
    for ax in axes[1]:
        ax.set_xlabel("Acquisition order")
        ax.set_ylabel("Force (pN)")
        ax.grid(alpha=0.2)
    axes[0, 0].legend(frameon=False)
    axes[1, 0].legend(frameon=False)
    time_lag_min = float(paired_rows[0]["start_time_lag_min"])
    figure.suptitle(
        f"Map 2 versus Map 3: same speed and pixels, {time_lag_min:.2f} min apart"
    )
    figure.tight_layout()
    figure.savefig(FIGURES / "map2_map3_paired_time_statistics.png", dpi=220)
    plt.close(figure)


def plot_same_location_time_maps(
    sources: list[base.SourceData], matrices: dict[str, np.ndarray]
) -> None:
    time_sources = sources[1:]
    time_matrices = [
        matrices[relative_source(source.path) + "|line"] for source in time_sources
    ]
    difference_count = len(time_sources) - 1
    row_count = len(time_sources) + difference_count
    figure, axes = plt.subplots(
        row_count,
        4,
        figsize=(16.8, 3.05 * row_count),
        layout="constrained",
    )
    row_labels = [f"Map {map_order}" for map_order in range(2, len(sources) + 1)]
    row_labels += [
        f"Map{map_order + 1}−Map{map_order}"
        for map_order in range(2, len(sources))
    ]
    for column, distance_nm in enumerate(TARGET_DISTANCES_NM):
        index = int(np.flatnonzero(np.isclose(BIN_CENTERS_NM, distance_nm))[0])
        physical_maps = [
            physical_matrix(source, matrix[:, index])
            for source, matrix in zip(time_sources, time_matrices, strict=True)
        ]
        differences = [
            physical_maps[index + 1] - physical_maps[index]
            for index in range(difference_count)
        ]
        combined = np.concatenate([matrix.ravel() for matrix in physical_maps])
        vmin, vmax = np.nanquantile(combined, [0.02, 0.98])
        difference_limit = max(
            float(
                np.nanquantile(
                    np.abs(np.concatenate([matrix.ravel() for matrix in differences])),
                    0.98,
                )
            ),
            1.0,
        )
        force_image = None
        for row, matrix in enumerate(physical_maps):
            force_image = axes[row, column].imshow(
                matrix, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax
            )
        difference_image = None
        for row_offset, matrix in enumerate(differences, start=len(time_sources)):
            difference_image = axes[row_offset, column].imshow(
                matrix,
                origin="lower",
                cmap="coolwarm",
                vmin=-difference_limit,
                vmax=difference_limit,
            )
        axes[0, column].set_title(f"D={distance_nm:g} nm")
        figure.colorbar(
            force_image,
            ax=list(axes[: len(time_sources), column]),
            shrink=0.72,
            pad=0.02,
            label="Force (pN)",
        )
        figure.colorbar(
            difference_image,
            ax=list(axes[len(time_sources) :, column]),
            shrink=0.72,
            pad=0.02,
            label="Paired difference (pN)",
        )
        for row in range(row_count):
            axes[row, column].set_xlabel("physical column")
            axes[row, column].set_ylabel("physical row")
    for row, label in enumerate(row_labels):
        axes[row, 0].set_ylabel(f"{label}\nphysical row")
    figure.suptitle("Same-location four-map time sequence: line-corrected force")
    figure.savefig(FIGURES / "map2_map5_time_force_maps.png", dpi=220)
    plt.close(figure)


def plot_same_location_time_statistics(
    sources: list[base.SourceData],
    matrices: dict[str, np.ndarray],
    time_trend_rows: list[dict],
) -> None:
    time_sources = sources[1:]
    time_matrices = [
        matrices[relative_source(source.path) + "|line"] for source in time_sources
    ]
    relative_times = np.asarray(
        [
            (source.timestamp - time_sources[0].timestamp).total_seconds() / 60.0
            for source in time_sources
        ]
    )
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 9.8))
    colors = ("#264653", "#e76f51", "#2a9d8f", "#7b2cbf")
    map_orders = tuple(range(2, len(sources) + 1))
    for map_order, matrix, color in zip(map_orders, time_matrices, colors, strict=True):
        median = np.nanmedian(matrix, axis=0)
        q25 = np.nanquantile(matrix, 0.25, axis=0)
        q75 = np.nanquantile(matrix, 0.75, axis=0)
        axes[0, 0].fill_between(BIN_CENTERS_NM, q25, q75, color=color, alpha=0.10)
        axes[0, 0].plot(
            BIN_CENTERS_NM, median, color=color, lw=2.0, label=f"Map {map_order}"
        )
    plotted_pairs = (
        (0, 1, "Map3−Map2", "#457b9d"),
        (1, 2, "Map4−Map3", "#e76f51"),
        (2, 3, "Map5−Map4", "#7b2cbf"),
        (0, 3, "Map5−Map2", "#2a9d8f"),
    )
    for before, after, label, color in plotted_pairs:
        difference = time_matrices[after] - time_matrices[before]
        median = np.nanmedian(difference, axis=0)
        q25 = np.nanquantile(difference, 0.25, axis=0)
        q75 = np.nanquantile(difference, 0.75, axis=0)
        axes[0, 1].fill_between(BIN_CENTERS_NM, q25, q75, color=color, alpha=0.09)
        axes[0, 1].plot(BIN_CENTERS_NM, median, color=color, lw=1.8, label=label)
    index50 = int(np.flatnonzero(np.isclose(BIN_CENTERS_NM, 50.0))[0])
    order = np.arange(1, 65)
    for map_order, matrix, color in zip(map_orders, time_matrices, colors, strict=True):
        axes[1, 0].plot(
            order,
            matrix[:, index50],
            "o-",
            ms=2.8,
            lw=0.75,
            color=color,
            label=f"Map {map_order}",
        )
    primary = [
        row
        for row in time_trend_rows
        if row["baseline_method"] == "linear_drift_corrected"
    ]
    distances = np.asarray([float(row["distance_nm"]) for row in primary])
    medians = np.asarray(
        [float(row["fitted_full_interval_change_median_pN"]) for row in primary]
    )
    lower = medians - np.asarray(
        [float(row["fitted_full_interval_change_q25_pN"]) for row in primary]
    )
    upper = np.asarray(
        [float(row["fitted_full_interval_change_q75_pN"]) for row in primary]
    ) - medians
    axes[1, 1].errorbar(
        distances,
        medians,
        yerr=np.vstack([lower, upper]),
        fmt="o-",
        capsize=4,
        color="#1d3557",
    )
    axes[0, 0].set_title("Same-location F–D distributions")
    axes[0, 1].set_title("Paired force differences")
    axes[1, 0].set_title("D=50 nm by acquisition order")
    axes[1, 1].set_title(
        f"Per-pixel fitted change over {relative_times[-1]:.2f} min; median and IQR"
    )
    for axis in axes[0]:
        axis.set_xlim(10, 250)
        axis.set_xlabel("Separation D (nm)")
        axis.set_ylabel("Force (pN)")
        axis.axhline(0.0, color="0.35", lw=0.8)
        axis.grid(alpha=0.2)
        axis.legend(frameon=False)
    axes[1, 0].set_xlabel("Acquisition order")
    axes[1, 0].set_ylabel("Force at 50 nm (pN)")
    axes[1, 0].grid(alpha=0.2)
    axes[1, 0].legend(frameon=False)
    axes[1, 1].set_xlabel("Separation D (nm)")
    axes[1, 1].set_ylabel("Fitted Map5−Map2 change (pN)")
    axes[1, 1].axhline(0.0, color="0.35", lw=0.8)
    axes[1, 1].grid(alpha=0.2)
    figure.suptitle("Map 2–5: same speed and pixels, four time points")
    figure.tight_layout()
    figure.savefig(FIGURES / "map2_map5_time_statistics.png", dpi=220)
    plt.close(figure)


def plot_absolute_force_time_statistics(
    sources: list[base.SourceData], matrices: dict[str, np.ndarray]
) -> list[dict]:
    """Plot mean per-pixel absolute force and its spatial spread versus time."""
    time_sources = sources[1:]
    if len(time_sources) != 4:
        raise RuntimeError("absolute-force time plot requires maps 2--5")
    relative_times_min = np.asarray(
        [
            (source.timestamp - time_sources[0].timestamp).total_seconds() / 60.0
            for source in time_sources
        ],
        dtype=np.float64,
    )
    if relative_times_min.shape != (4,) or not np.all(np.isfinite(relative_times_min)):
        raise RuntimeError("invalid map2--5 relative times")

    rows: list[dict] = []
    by_distance: dict[float, dict[str, np.ndarray]] = {}
    for distance_nm in TARGET_DISTANCES_NM:
        distance_index = int(
            np.flatnonzero(np.isclose(BIN_CENTERS_NM, distance_nm))[0]
        )
        means: list[float] = []
        standard_deviations: list[float] = []
        variances: list[float] = []
        for map_order, source, relative_time_min in zip(
            range(2, 6), time_sources, relative_times_min, strict=True
        ):
            signed_force_pN = np.asarray(
                matrices[relative_source(source.path) + "|line"][:, distance_index],
                dtype=np.float64,
            )
            if signed_force_pN.shape != (64,) or not np.all(
                np.isfinite(signed_force_pN)
            ):
                raise RuntimeError(
                    f"map{map_order} D={distance_nm:g} nm requires 64 finite forces"
                )
            absolute_force_pN = np.abs(signed_force_pN)
            mean_pN = float(np.mean(absolute_force_pN))
            std_pN = float(np.std(absolute_force_pN, ddof=1))
            variance_pN2 = float(np.var(absolute_force_pN, ddof=1))
            if not np.isclose(variance_pN2, std_pN**2, rtol=1e-13, atol=1e-12):
                raise RuntimeError("sample variance is inconsistent with sample SD")
            means.append(mean_pN)
            standard_deviations.append(std_pN)
            variances.append(variance_pN2)
            rows.append(
                {
                    "source": relative_source(source.path),
                    "map_order": map_order,
                    "relative_time_min": float(relative_time_min),
                    "distance_nm": distance_nm,
                    "pixels": 64,
                    "baseline_method": "linear_drift_corrected",
                    "absolute_force_definition": "absolute value per pixel before aggregation",
                    "mean_absolute_force_pN": mean_pN,
                    "median_absolute_force_pN": float(
                        np.median(absolute_force_pN)
                    ),
                    "sample_std_absolute_force_pN": std_pN,
                    "sample_variance_absolute_force_pN2": variance_pN2,
                    "sem_absolute_force_pN": float(std_pN / np.sqrt(64.0)),
                }
            )
        by_distance[distance_nm] = {
            "mean": np.asarray(means, dtype=np.float64),
            "std": np.asarray(standard_deviations, dtype=np.float64),
            "variance": np.asarray(variances, dtype=np.float64),
        }

    colors = {
        20.0: "#d1495b",
        50.0: "#e07a1f",
        100.0: "#2a9d8f",
        200.0: "#386cb0",
    }
    figure, axes = plt.subplots(
        3,
        1,
        figsize=(10.8, 12.2),
        sharex=True,
        gridspec_kw={"height_ratios": (1.55, 1.0, 1.0)},
    )
    for distance_nm in TARGET_DISTANCES_NM:
        statistics = by_distance[distance_nm]
        color = colors[distance_nm]
        mean = statistics["mean"]
        std = statistics["std"]
        variance = statistics["variance"]
        label = f"D = {distance_nm:g} nm"
        axes[0].plot(
            relative_times_min,
            mean,
            "o-",
            color=color,
            linewidth=2.2,
            markersize=6,
            label=label,
        )
        axes[0].fill_between(
            relative_times_min,
            np.maximum(mean - std, 0.0),
            mean + std,
            color=color,
            alpha=0.13,
            linewidth=0,
        )
        axes[1].plot(
            relative_times_min,
            std,
            "o-",
            color=color,
            linewidth=2.0,
            markersize=5,
            label=label,
        )
        axes[2].plot(
            relative_times_min,
            variance,
            "o-",
            color=color,
            linewidth=2.0,
            markersize=5,
            label=label,
        )

    axes[0].set_ylabel(r"Mean $|F_i|$ (pN)")
    axes[0].set_title(
        "Absolute force versus experimental time\n"
        "Lines: 64-pixel mean; shaded regions: ±1 sample SD"
    )
    axes[0].legend(frameon=False, ncol=2)
    axes[1].set_ylabel(r"SD of $|F_i|$ (pN)")
    axes[1].set_title("Spatial sample standard deviation (ddof = 1)")
    axes[2].set_ylabel(r"Variance of $|F_i|$ (pN$^2$)")
    axes[2].set_title("Spatial sample variance (ddof = 1)")
    axes[2].set_xlabel("Time from start of Map 2 (min)")
    axes[2].set_xticks(relative_times_min)
    axes[2].set_xticklabels([f"{value:.2f}" for value in relative_times_min])
    for axis in axes:
        axis.set_xlim(relative_times_min[0] - 2.5, relative_times_min[-1] + 2.5)
        axis.set_ylim(bottom=0.0)
        axis.grid(alpha=0.24)
    figure.tight_layout()
    figure.savefig(FIGURES / "absolute_force_vs_time_mean_std_variance.png", dpi=240)
    plt.close(figure)
    return rows


def force_value(
    rows: list[dict], source: str, method: str, distance_nm: float
) -> dict:
    selected = [
        row
        for row in rows
        if row["source"] == source
        and row["baseline_method"] == method
        and np.isclose(
            float(row["distance_nm"]), distance_nm, rtol=0.0, atol=1e-9
        )
    ]
    if len(selected) != 1:
        raise RuntimeError(
            f"missing force row for {source}, {method}, {distance_nm} nm"
        )
    return selected[0]


def distribution_comparisons(
    sources: list[base.SourceData], matrices: dict[str, np.ndarray]
) -> list[dict]:
    """Describe map-to-map pixel separation without claiming independent replicates."""

    rows: list[dict] = []
    for method, suffix in (
        ("linear_drift_corrected", "|line"),
        ("far_constant_referenced", "|constant"),
    ):
        first = matrices[relative_source(sources[0].path) + suffix]
        second = matrices[relative_source(sources[1].path) + suffix]
        for distance_nm in TARGET_DISTANCES_NM:
            index = int(np.flatnonzero(np.isclose(BIN_CENTERS_NM, distance_nm))[0])
            a = first[:, index]
            b = second[:, index]
            a = a[np.isfinite(a)]
            b = b[np.isfinite(b)]
            welch = ttest_ind(a, b, equal_var=False)
            mann_whitney = mannwhitneyu(a, b, alternative="two-sided")
            cliff_delta = float(
                (
                    np.count_nonzero(b[:, None] > a[None, :])
                    - np.count_nonzero(b[:, None] < a[None, :])
                )
                / (a.size * b.size)
            )
            pooled_sd = math.sqrt(
                (
                    (a.size - 1) * float(np.var(a, ddof=1))
                    + (b.size - 1) * float(np.var(b, ddof=1))
                )
                / (a.size + b.size - 2)
            )
            cohen_d = (
                float((np.mean(b) - np.mean(a)) / pooled_sd)
                if pooled_sd > np.finfo(np.float64).tiny
                else float("nan")
            )
            rows.append(
                {
                    "baseline_method": method,
                    "distance_nm": distance_nm,
                    "map1_pixels": int(a.size),
                    "map2_pixels": int(b.size),
                    "map1_median_pN": float(np.median(a)),
                    "map2_median_pN": float(np.median(b)),
                    "map2_minus_map1_median_pN": float(
                        np.median(b) - np.median(a)
                    ),
                    "Welch_t_statistic_map1_minus_map2": float(welch.statistic),
                    "Welch_t_p_value_naive_pixels": float(welch.pvalue),
                    "Mann_Whitney_U_statistic": float(mann_whitney.statistic),
                    "Mann_Whitney_p_value_naive_pixels": float(
                        mann_whitney.pvalue
                    ),
                    "Cliffs_delta_map2_greater_map1": cliff_delta,
                    "Cohen_d_map2_minus_map1": cohen_d,
                }
            )
    return rows


def paired_time_comparisons(
    sources: list[base.SourceData], matrices: dict[str, np.ndarray]
) -> tuple[list[dict], dict[str, float | bool]]:
    """Paired map3-map2 comparison at exactly matched physical pixels."""

    if len(sources) != 3:
        raise RuntimeError("paired time comparison requires exactly three maps")
    map2_coordinates = physical_coordinates_um(sources[1].path, 64)
    map3_coordinates = physical_coordinates_um(sources[2].path, 64)
    coordinate_difference_nm = np.abs(map3_coordinates - map2_coordinates) * 1e3
    max_coordinate_difference_nm = float(np.max(coordinate_difference_nm))
    coordinates_match = bool(max_coordinate_difference_nm <= 1e-6)
    if not coordinates_match:
        raise RuntimeError(
            "map2/map3 point coordinates do not match: "
            f"max difference {max_coordinate_difference_nm:g} nm"
        )
    time_lag_min = float(
        (sources[2].timestamp - sources[1].timestamp).total_seconds() / 60.0
    )
    metadata = {
        "coordinates_match_pointwise": coordinates_match,
        "maximum_pointwise_coordinate_difference_nm": max_coordinate_difference_nm,
        "start_time_lag_min": time_lag_min,
    }
    rows: list[dict] = []
    acquisition_order = np.arange(1, 65, dtype=np.float64)
    for method, suffix in (
        ("linear_drift_corrected", "|line"),
        ("far_constant_referenced", "|constant"),
    ):
        map2 = matrices[relative_source(sources[1].path) + suffix]
        map3 = matrices[relative_source(sources[2].path) + suffix]
        for distance_nm in TARGET_DISTANCES_NM:
            index = int(np.flatnonzero(np.isclose(BIN_CENTERS_NM, distance_nm))[0])
            before = map2[:, index]
            after = map3[:, index]
            valid = np.isfinite(before) & np.isfinite(after)
            before = before[valid]
            after = after[valid]
            order = acquisition_order[valid]
            if before.size < 56:
                raise RuntimeError(
                    f"too few paired map2/map3 pixels at {distance_nm:g} nm"
                )
            difference = after - before
            full_difference = np.full(64, np.nan, dtype=np.float64)
            full_difference[valid] = difference
            physical_difference = physical_matrix(sources[1], full_difference)
            row_medians = np.nanmedian(physical_difference, axis=1)
            paired_t = ttest_rel(after, before)
            signed_rank = wilcoxon(after, before, alternative="two-sided")
            row_t = ttest_1samp(row_medians, 0.0)
            row_signed_rank = wilcoxon(
                row_medians, np.zeros_like(row_medians), alternative="two-sided"
            )
            pearson = pearsonr(before, after)
            spearman = spearmanr(before, after)
            time_slope, _, time_low, time_high = theilslopes(
                difference, order, alpha=0.95
            )
            map_relation_slope, map_relation_intercept, _, _ = theilslopes(
                after, before, alpha=0.95
            )
            difference_sd = float(np.std(difference, ddof=1))
            rows.append(
                {
                    "baseline_method": method,
                    "distance_nm": distance_nm,
                    "paired_pixels": int(before.size),
                    "map2_median_pN": float(np.median(before)),
                    "map3_median_pN": float(np.median(after)),
                    "paired_difference_map3_minus_map2_mean_pN": float(
                        np.mean(difference)
                    ),
                    "paired_difference_map3_minus_map2_median_pN": float(
                        np.median(difference)
                    ),
                    "paired_difference_q25_pN": float(
                        np.quantile(difference, 0.25)
                    ),
                    "paired_difference_q75_pN": float(
                        np.quantile(difference, 0.75)
                    ),
                    "paired_difference_sd_pN": difference_sd,
                    "paired_Cohen_dz": float(np.mean(difference) / difference_sd)
                    if difference_sd > np.finfo(np.float64).tiny
                    else float("nan"),
                    "fraction_map3_greater_map2": float(np.mean(difference > 0.0)),
                    "paired_t_statistic": float(paired_t.statistic),
                    "paired_t_p_value": float(paired_t.pvalue),
                    "Wilcoxon_statistic": float(signed_rank.statistic),
                    "Wilcoxon_p_value": float(signed_rank.pvalue),
                    "physical_row_median_difference_mean_pN": float(
                        np.mean(row_medians)
                    ),
                    "physical_row_median_difference_sd_pN": float(
                        np.std(row_medians, ddof=1)
                    ),
                    "row_block_t_p_value": float(row_t.pvalue),
                    "row_block_Wilcoxon_p_value": float(row_signed_rank.pvalue),
                    "Pearson_r_map2_map3": float(pearson.statistic),
                    "Pearson_p_value": float(pearson.pvalue),
                    "Spearman_rho_map2_map3": float(spearman.statistic),
                    "Spearman_p_value": float(spearman.pvalue),
                    "Theil_Sen_delta_endpoint_change_pN": float(time_slope * 63.0),
                    "Theil_Sen_delta_endpoint_low_pN": float(time_low * 63.0),
                    "Theil_Sen_delta_endpoint_high_pN": float(time_high * 63.0),
                    "Theil_Sen_map3_vs_map2_slope": float(map_relation_slope),
                    "Theil_Sen_map3_vs_map2_intercept_pN": float(
                        map_relation_intercept
                    ),
                    "start_time_lag_min": time_lag_min,
                }
            )
    return rows, metadata


def baseline_coupling(
    sources: list[base.SourceData],
    pixel_rows: list[dict],
    matrices: dict[str, np.ndarray],
) -> list[dict]:
    """Quantify coupling of reconstructed force to fitted baseline/contact QC."""

    rows: list[dict] = []
    for map_order, source in enumerate(sources, start=1):
        selected = sorted(
            (row for row in pixel_rows if int(row["map_order"]) == map_order),
            key=lambda row: int(row["point_index"]),
        )
        if [int(row["point_index"]) for row in selected] != list(range(64)):
            raise RuntimeError(f"map {map_order}: pixel QC order is incomplete")
        metrics = {
            "far_slope_pN_per_100nm": np.asarray(
                [float(row["far_slope_pN_per_100nm"]) for row in selected]
            ),
            "local_contact_InvOLS_nm_per_V": np.asarray(
                [float(row["local_contact_InvOLS_nm_per_V"]) for row in selected]
            ),
            "acquisition_order": np.asarray(
                [float(row["pixel_acquisition_order"]) for row in selected]
            ),
            "physical_row": np.asarray([float(row["row"]) for row in selected]),
            "physical_column": np.asarray(
                [float(row["column"]) for row in selected]
            ),
        }
        for method, suffix in (
            ("linear_drift_corrected", "|line"),
            ("far_constant_referenced", "|constant"),
        ):
            matrix = matrices[relative_source(source.path) + suffix]
            for distance_nm in TARGET_DISTANCES_NM:
                index = int(
                    np.flatnonzero(np.isclose(BIN_CENTERS_NM, distance_nm))[0]
                )
                force = matrix[:, index]
                for metric_name, metric in metrics.items():
                    valid = np.isfinite(force) & np.isfinite(metric)
                    correlation = spearmanr(force[valid], metric[valid])
                    rows.append(
                        {
                            "map_order": map_order,
                            "baseline_method": method,
                            "distance_nm": distance_nm,
                            "QC_metric": metric_name,
                            "pixels": int(np.count_nonzero(valid)),
                            "Spearman_rho": float(correlation.statistic),
                            "Spearman_p_value_naive_pixels": float(
                                correlation.pvalue
                            ),
                        }
                    )
    return rows


def paired_qc_comparisons(pixel_rows: list[dict]) -> list[dict]:
    map2 = sorted(
        (row for row in pixel_rows if int(row["map_order"]) == 2),
        key=lambda row: int(row["point_index"]),
    )
    map3 = sorted(
        (row for row in pixel_rows if int(row["map_order"]) == 3),
        key=lambda row: int(row["point_index"]),
    )
    if len(map2) != 64 or len(map3) != 64:
        raise RuntimeError("paired QC comparison requires 64 pixels in map2/map3")
    metrics = (
        ("local_contact_InvOLS", "local_contact_InvOLS_nm_per_V", "nm/V"),
        ("far_slope", "far_slope_pN_per_100nm", "pN/100 nm"),
        ("far_noise", "far_noise_pN", "pN"),
        ("terminal_load", "terminal_load_nN", "nN"),
        ("gap_speed", "gap_speed_20_200nm_um_per_s", "um/s"),
        ("retract_pull_off_force", "retract_pull_off_force_nN", "nN"),
        (
            "retract_detachment_travel",
            "retract_detachment_piezo_travel_nm",
            "nm",
        ),
    )
    rows: list[dict] = []
    for metric, field, unit in metrics:
        before = np.asarray([float(row[field]) for row in map2], dtype=np.float64)
        after = np.asarray([float(row[field]) for row in map3], dtype=np.float64)
        valid = np.isfinite(before) & np.isfinite(after)
        before = before[valid]
        after = after[valid]
        difference = after - before
        signed_rank = wilcoxon(difference, alternative="two-sided")
        correlation = spearmanr(before, after)
        rows.append(
            {
                "metric": metric,
                "unit": unit,
                "paired_pixels": int(before.size),
                "map2_median": float(np.median(before)),
                "map3_median": float(np.median(after)),
                "paired_difference_map3_minus_map2_median": float(
                    np.median(difference)
                ),
                "paired_difference_q25": float(np.quantile(difference, 0.25)),
                "paired_difference_q75": float(np.quantile(difference, 0.75)),
                "Wilcoxon_p_value": float(signed_rank.pvalue),
                "Spearman_rho_map2_map3": float(correlation.statistic),
                "Spearman_p_value": float(correlation.pvalue),
            }
        )
    return rows


def same_location_time_comparisons(
    sources: list[base.SourceData], matrices: dict[str, np.ndarray]
) -> tuple[list[dict], list[dict], dict]:
    """Pairwise and four-timepoint statistics for same-location maps 2--5."""

    if len(sources) != 5:
        raise RuntimeError("four-map time comparison requires five total maps")
    time_sources = sources[1:]
    reference_coordinates = physical_coordinates_um(time_sources[0].path, 64)
    maximum_coordinate_difference_nm = 0.0
    for source in time_sources[1:]:
        coordinates = physical_coordinates_um(source.path, 64)
        maximum_coordinate_difference_nm = max(
            maximum_coordinate_difference_nm,
            float(np.max(np.abs(coordinates - reference_coordinates)) * 1e3),
        )
    if maximum_coordinate_difference_nm > 1e-6:
        raise RuntimeError(
            "map2/map3/map4/map5 coordinates do not match: "
            f"max difference {maximum_coordinate_difference_nm:g} nm"
        )
    relative_times_min = np.asarray(
        [
            (source.timestamp - time_sources[0].timestamp).total_seconds() / 60.0
            for source in time_sources
        ],
        dtype=np.float64,
    )
    metadata = {
        "maps": [2, 3, 4, 5],
        "coordinates_match_pointwise": True,
        "maximum_pointwise_coordinate_difference_nm": maximum_coordinate_difference_nm,
        "relative_start_times_min": relative_times_min.tolist(),
        "adjacent_start_time_lags_min": np.diff(relative_times_min).tolist(),
        "total_start_time_lag_min": float(relative_times_min[-1]),
    }

    pair_rows: list[dict] = []
    pairs = tuple(
        (before_index, after_index)
        for before_index in range(1, len(sources) - 1)
        for after_index in range(before_index + 1, len(sources))
    )
    acquisition_order = np.arange(1, 65, dtype=np.float64)
    for before_index, after_index in pairs:
        before_map = before_index + 1
        after_map = after_index + 1
        lag_min = float(relative_times_min[after_index - 1] - relative_times_min[before_index - 1])
        for method, suffix in (
            ("linear_drift_corrected", "|line"),
            ("far_constant_referenced", "|constant"),
        ):
            before_matrix = matrices[
                relative_source(sources[before_index].path) + suffix
            ]
            after_matrix = matrices[
                relative_source(sources[after_index].path) + suffix
            ]
            for distance_nm in TARGET_DISTANCES_NM:
                index = int(
                    np.flatnonzero(np.isclose(BIN_CENTERS_NM, distance_nm))[0]
                )
                before_full = before_matrix[:, index]
                after_full = after_matrix[:, index]
                valid = np.isfinite(before_full) & np.isfinite(after_full)
                before = before_full[valid]
                after = after_full[valid]
                order = acquisition_order[valid]
                difference = after - before
                full_difference = np.full(64, np.nan, dtype=np.float64)
                full_difference[valid] = difference
                row_medians = np.nanmedian(
                    physical_matrix(sources[before_index], full_difference), axis=1
                )
                paired_t = ttest_rel(after, before)
                signed_rank = wilcoxon(difference, alternative="two-sided")
                row_signed_rank = wilcoxon(
                    row_medians,
                    np.zeros_like(row_medians),
                    alternative="two-sided",
                )
                pearson = pearsonr(before, after)
                spearman = spearmanr(before, after)
                order_slope, _, order_low, order_high = theilslopes(
                    difference, order, alpha=0.95
                )
                pair_rows.append(
                    {
                        "before_map": before_map,
                        "after_map": after_map,
                        "pair": f"map{before_map}_to_map{after_map}",
                        "baseline_method": method,
                        "distance_nm": distance_nm,
                        "time_lag_min": lag_min,
                        "paired_pixels": int(before.size),
                        "before_median_pN": float(np.median(before)),
                        "after_median_pN": float(np.median(after)),
                        "paired_difference_after_minus_before_mean_pN": float(
                            np.mean(difference)
                        ),
                        "paired_difference_after_minus_before_median_pN": float(
                            np.median(difference)
                        ),
                        "paired_difference_q25_pN": float(
                            np.quantile(difference, 0.25)
                        ),
                        "paired_difference_q75_pN": float(
                            np.quantile(difference, 0.75)
                        ),
                        "fraction_after_greater_before": float(
                            np.mean(difference > 0.0)
                        ),
                        "paired_t_p_value": float(paired_t.pvalue),
                        "Wilcoxon_p_value": float(signed_rank.pvalue),
                        "row_block_Wilcoxon_p_value": float(
                            row_signed_rank.pvalue
                        ),
                        "Pearson_r_before_after": float(pearson.statistic),
                        "Spearman_rho_before_after": float(spearman.statistic),
                        "Theil_Sen_delta_endpoint_change_pN": float(
                            order_slope * 63.0
                        ),
                        "Theil_Sen_delta_endpoint_low_pN": float(
                            order_low * 63.0
                        ),
                        "Theil_Sen_delta_endpoint_high_pN": float(
                            order_high * 63.0
                        ),
                    }
                )

    trend_rows: list[dict] = []
    centered_time = relative_times_min - np.mean(relative_times_min)
    time_denominator = float(np.sum(centered_time**2))
    for method, suffix in (
        ("linear_drift_corrected", "|line"),
        ("far_constant_referenced", "|constant"),
    ):
        time_matrices = np.asarray(
            [matrices[relative_source(source.path) + suffix] for source in time_sources]
        )
        for distance_nm in TARGET_DISTANCES_NM:
            index = int(np.flatnonzero(np.isclose(BIN_CENTERS_NM, distance_nm))[0])
            values = time_matrices[:, :, index]
            if not np.all(np.isfinite(values)):
                raise RuntimeError(f"non-finite four-map time data at {distance_nm:g} nm")
            slopes = np.sum(centered_time[:, None] * values, axis=0) / time_denominator
            full_change = slopes * float(relative_times_min[-1])
            physical_change = physical_matrix(time_sources[0], full_change)
            row_medians = np.nanmedian(physical_change, axis=1)
            signed_rank = wilcoxon(full_change, alternative="two-sided")
            row_signed_rank = wilcoxon(
                row_medians, np.zeros_like(row_medians), alternative="two-sided"
            )
            decreasing = np.all(np.diff(values, axis=0) < 0.0, axis=0)
            increasing = np.all(np.diff(values, axis=0) > 0.0, axis=0)
            trend_row = {
                    "baseline_method": method,
                    "distance_nm": distance_nm,
                    "map2_median_pN": float(np.median(values[0])),
                    "map3_median_pN": float(np.median(values[1])),
                    "map4_median_pN": float(np.median(values[2])),
                    "map5_median_pN": float(np.median(values[3])),
                    "pixel_time_slope_median_pN_per_min": float(np.median(slopes)),
                    "pixel_time_slope_q25_pN_per_min": float(
                        np.quantile(slopes, 0.25)
                    ),
                    "pixel_time_slope_q75_pN_per_min": float(
                        np.quantile(slopes, 0.75)
                    ),
                    "fitted_full_interval_change_median_pN": float(
                        np.median(full_change)
                    ),
                    "fitted_full_interval_change_q25_pN": float(
                        np.quantile(full_change, 0.25)
                    ),
                    "fitted_full_interval_change_q75_pN": float(
                        np.quantile(full_change, 0.75)
                    ),
                    "fraction_negative_time_slope": float(np.mean(slopes < 0.0)),
                    "fraction_strict_monotonic_decrease": float(
                        np.mean(decreasing)
                    ),
                    "fraction_strict_monotonic_increase": float(
                        np.mean(increasing)
                    ),
                    "pixel_slope_Wilcoxon_p_value": float(signed_rank.pvalue),
                    "row_block_slope_Wilcoxon_p_value": float(
                        row_signed_rank.pvalue
                    ),
                    "total_time_span_min": float(relative_times_min[-1]),
                }
            trend_rows.append(trend_row)
    return pair_rows, trend_rows, metadata


def same_location_paired_qc_comparisons(pixel_rows: list[dict]) -> list[dict]:
    maps = {
        map_order: sorted(
            (row for row in pixel_rows if int(row["map_order"]) == map_order),
            key=lambda row: int(row["point_index"]),
        )
        for map_order in (2, 3, 4, 5)
    }
    if any(len(rows) != 64 for rows in maps.values()):
        raise RuntimeError("four-map paired QC requires 64 pixels per map")
    metrics = (
        ("local_contact_InvOLS", "local_contact_InvOLS_nm_per_V", "nm/V"),
        ("far_slope", "far_slope_pN_per_100nm", "pN/100 nm"),
        ("far_noise", "far_noise_pN", "pN"),
        ("terminal_load", "terminal_load_nN", "nN"),
        ("gap_speed", "gap_speed_20_200nm_um_per_s", "um/s"),
        ("retract_pull_off_force", "retract_pull_off_force_nN", "nN"),
        ("retract_detachment_travel", "retract_detachment_piezo_travel_nm", "nm"),
    )
    rows: list[dict] = []
    pairs = tuple(
        (before_map, after_map)
        for before_map in (2, 3, 4)
        for after_map in range(before_map + 1, 6)
    )
    for before_map, after_map in pairs:
        for metric, field, unit in metrics:
            before = np.asarray(
                [float(row[field]) for row in maps[before_map]], dtype=np.float64
            )
            after = np.asarray(
                [float(row[field]) for row in maps[after_map]], dtype=np.float64
            )
            valid = np.isfinite(before) & np.isfinite(after)
            before = before[valid]
            after = after[valid]
            difference = after - before
            signed_rank = wilcoxon(difference, alternative="two-sided")
            correlation = spearmanr(before, after)
            rows.append(
                {
                    "before_map": before_map,
                    "after_map": after_map,
                    "pair": f"map{before_map}_to_map{after_map}",
                    "metric": metric,
                    "unit": unit,
                    "paired_pixels": int(before.size),
                    "before_median": float(np.median(before)),
                    "after_median": float(np.median(after)),
                    "paired_difference_after_minus_before_median": float(
                        np.median(difference)
                    ),
                    "paired_difference_q25": float(np.quantile(difference, 0.25)),
                    "paired_difference_q75": float(np.quantile(difference, 0.75)),
                    "Wilcoxon_p_value": float(signed_rank.pvalue),
                    "Spearman_rho_before_after": float(correlation.statistic),
                }
            )
    return rows


def render_report_legacy_two_maps(
    inferred_cantilever: str,
    spring_N_per_m: float,
    spring_sd_N_per_m: float,
    fingerprints: dict[str, dict[str, float | str]],
    global_sensitivity: float,
    map_rows: list[dict],
    force_rows: list[dict],
    trend_rows: list[dict],
    contact_rows: list[dict],
    comparison_rows: list[dict],
    coupling_rows: list[dict],
) -> str:
    center_distance_um = math.hypot(
        float(map_rows[1]["x_center_um"]) - float(map_rows[0]["x_center_um"]),
        float(map_rows[1]["y_center_um"]) - float(map_rows[0]["y_center_um"]),
    )
    water_viscosity_Pa_s = (
        base.cheng_viscosity_mPa_s(0.0, CONTEXT_TEMPERATURE_C) * 1e-3
    )
    delta_gap_speed_m_per_s = (
        float(map_rows[1]["gap_speed_20_200nm_median_um_per_s"])
        - float(map_rows[0]["gap_speed_20_200nm_median_um_per_s"])
    ) * 1e-6
    lines = [
        "# 31-08-26 两张 force map：raw reconstruction 与基础比较",
        "",
        "## 直接结论",
        "",
        f"两张文件都是完整8×8 map。raw metadata显示approach/retract单段分别为2.0/2.0 s和1.0/1.0 s；重建后的20–200 nm median gap speed见下表。两图的XY范围不重叠，因此它们不是同一区域的same-pixel velocity pair，不能用paired t-test或把差值直接归因于速度。",
        "",
        f"文件没有显式cantilever名称。stored `(InvOLS, k)`与D4-D6原始标定文件逐一比较后，最近指纹为 **{inferred_cantilever}**；因此force采用25.6 °C空气独立标定 `k={spring_N_per_m:.9f} N/m`，repeatability SD `{spring_sd_N_per_m:.6f} N/m`。这是metadata inference，不是文件内明确标签。map文件写入的`k=0.350609 N/m`未用于force。",
        "",
        f"两张map的全部hard-contact curve在40 nm末端窗口重新拟合，共同得到global liquid InvOLS **{global_sensitivity*1e9:.3f} nm/V**；文件写入的`78.275 nm/V`未用于force。两图local contact InvOLS差异较大，作为主要QC保留。",
        "",
        "液体组成、温度、球/平面材料和是否同一probe没有编码在map中；本报告只做校准、force reconstruction与map/systematic QC，不做hydrodynamic subtraction、PB、zeta potential或Debye length拟合。",
        "",
        "## Map inventory",
        "",
        "| map | time | field | XY center (µm) | approach/gap speed (µm/s) | local/global InvOLS (nm/V) | terminal load (nN) | far slope (pN/100 nm) |",
        "|---:|:---|:---|:---|:---|:---|---:|---:|",
    ]
    for row in map_rows:
        lines.append(
            f"| {row['map_order']} | {str(row['timestamp'])[11:19]} | "
            f"{row['field_u_um']:.1f}×{row['field_v_um']:.1f} µm | "
            f"({row['x_center_um']:.2f}, {row['y_center_um']:.2f}) | "
            f"{row['scanner_speed_median_um_per_s']:.3f} / "
            f"{row['gap_speed_20_200nm_median_um_per_s']:.3f} | "
            f"{row['map_contact_InvOLS_nm_per_V']:.3f} / "
            f"{row['global_InvOLS_used_nm_per_V']:.3f} | "
            f"{row['terminal_load_median_nN']:.3f} | "
            f"{row['far_slope_median_pN_per_100nm']:.2f} |"
        )
    lines += [
        "",
        f"两张map中心相距约{center_distance_um:.2f} µm，且各自半宽不足1.8 µm，因此没有空间重叠。map2不是map1的zoom-in。",
        "",
        "## Sensitivity QC",
        "",
        "contact window从30到50 nm时，每张map的median InvOLS只缓慢变化；map间差异则持续存在，因此不是单一contact-window选择造成的。括号内为40 nm primary fit最终进入global consensus的curve数。",
        "",
        "| map | 30 nm window | 40 nm window | 50 nm window | retained / 64 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for map_row in map_rows:
        selected = [
            row for row in contact_rows if row["source"] == map_row["source"]
        ]
        medians = [
            float(
                np.median(
                    np.asarray(
                        [
                            row[f"contact_invOLS_{span}nm_nm_per_V"]
                            for row in selected
                        ],
                        dtype=np.float64,
                    )
                )
            )
            for span in (30, 40, 50)
        ]
        retained = sum(bool(row["retained_for_global_consensus"]) for row in selected)
        lines.append(
            f"| {map_row['map_order']} | {medians[0]:.3f} | {medians[1]:.3f} | "
            f"{medians[2]:.3f} | {retained} / 64 |"
        )
    lines += [
        "",
        "## Measured force slices",
        "",
        "以下均为64-pixel median，单位pN。Primary为逐curve far-field linear baseline subtraction；constant branch只减远场常数而保留斜率。",
        "",
        "| D (nm) | map1 linear | map2 linear | map1 constant | map2 constant |",
        "|---:|---:|---:|---:|---:|",
    ]
    for distance in TARGET_DISTANCES_NM:
        values = []
        for method in ("linear_drift_corrected", "far_constant_referenced"):
            for row in map_rows:
                values.append(
                    float(force_value(force_rows, row["source"], method, distance)["force_median_pN"])
                )
        lines.append(
            f"| {distance:.0f} | {values[0]:.2f} | {values[1]:.2f} | "
            f"{values[2]:.2f} | {values[3]:.2f} |"
        )
    lines += [
        "",
        "这些是finite-speed、不同区域的measured force，不是equilibrium force。map间差值同时包含速度、位置、21 min时间间隔、contact response变化和表面/探针history。",
        "",
        "在20和50 nm，map2−map1 median分别为+249.42和+91.22 pN（far-linear）；改用far-constant后为+236.24和+96.21 pN。因此近场map差值对这两种zero-force定义较稳健。相反，100和200 nm的符号或大小明显依赖baseline方法，不能作为稳健远场力结论。",
        "",
        "## Pixel-distribution separation（描述性）",
        "",
        "下表按64个pixel做Welch t-test与Mann–Whitney test，并给出Cliff's delta。由于pixel存在空间相关、两张map又不在同一区域，p值是假定pixel独立的naive数值，只说明当前两张图的pixel分布是否重叠，不能检验纯velocity effect。",
        "",
        "| D (nm) | map2−map1 median (pN) | Welch p | Mann–Whitney p | Cliff's delta |",
        "|---:|---:|---:|---:|---:|",
    ]
    primary_comparisons = [
        row
        for row in comparison_rows
        if row["baseline_method"] == "linear_drift_corrected"
    ]
    for row in primary_comparisons:
        lines.append(
            f"| {row['distance_nm']:.0f} | "
            f"{row['map2_minus_map1_median_pN']:+.2f} | "
            f"{row['Welch_t_p_value_naive_pixels']:.3g} | "
            f"{row['Mann_Whitney_p_value_naive_pixels']:.3g} | "
            f"{row['Cliffs_delta_map2_greater_map1']:+.3f} |"
        )
    lines += [
        "",
        "20和50 nm的两组pixel分布明显分离；100和200 nm在far-linear定义下不分离。这个结论是map-level描述，不是速度因果结论。",
        "",
        "## Conditional hydrodynamic scale check",
        "",
        "若沿用此前给定的25.6 °C纯水与sphere radius 4.547 µm，并取no-slip sphere-plane lubrication近似 `ΔF_hyd=6πηR²Δv/D`，两图实测gap-speed差对应：",
        "",
        "| D (nm) | expected ΔF_hyd (pN) | measured Δmedian, far-linear (pN) | measured / hyd |",
        "|---:|---:|---:|---:|",
    ]
    for distance_nm in TARGET_DISTANCES_NM:
        hydrodynamic_pN = (
            6.0
            * math.pi
            * water_viscosity_Pa_s
            * CONTEXT_PROBE_RADIUS_M**2
            * delta_gap_speed_m_per_s
            / (distance_nm * 1e-9)
            * 1e12
        )
        selected = [
            row
            for row in primary_comparisons
            if float(row["distance_nm"]) == distance_nm
        ]
        if len(selected) != 1:
            raise RuntimeError("missing hydrodynamic scale-comparison row")
        measured = float(selected[0]["map2_minus_map1_median_pN"])
        lines.append(
            f"| {distance_nm:.0f} | {hydrodynamic_pN:.2f} | "
            f"{measured:+.2f} | {measured / hydrodynamic_pN:+.1f}× |"
        )
    lines += [
        "",
        "在20–50 nm，measured map difference约为该理想hydrodynamic increment的13–15倍；因此即使液体和probe假设正确，也不能把map差直接解释成hydrodynamic force。100–200 nm又受baseline定义强烈影响，不适合用这个比值下结论。这里仅作量级检查，没有从数据中拟合hydrodynamic coefficient。",
        "",
        "## Baseline coupling",
        "",
        "far-linear subtraction后，pixel force与该pixel拟合出的far slope高度相关；far-constant branch基本消除这一耦合。下表为D=50 nm的Spearman rho：",
        "",
        "| map | far-linear rho(F, far slope) | far-constant rho(F, far slope) |",
        "|---:|---:|---:|",
    ]
    for map_row in map_rows:
        values = []
        for method in ("linear_drift_corrected", "far_constant_referenced"):
            selected = [
                row
                for row in coupling_rows
                if int(row["map_order"]) == int(map_row["map_order"])
                and row["baseline_method"] == method
                and float(row["distance_nm"]) == 50.0
                and row["QC_metric"] == "far_slope_pN_per_100nm"
            ]
            if len(selected) != 1:
                raise RuntimeError("missing D=50 nm baseline-coupling row")
            values.append(float(selected[0]["Spearman_rho"]))
        lines.append(
            f"| {map_row['map_order']} | {values[0]:+.3f} | {values[1]:+.3f} |"
        )
    lines += [
        "",
        "这不是独立的物理相关性：far slope本身进入了linear subtraction算子，因此高rho主要是在提示pixel-level远场图样被baseline外推控制。近场map median在两种baseline下接近，但100–200 nm的pixel图样和absolute zero必须谨慎解释。",
        "",
        "## Within-map acquisition-order change",
        "",
        "`Theil–Sen endpoint change`是64个取样序号的robust线性趋势乘63；仍混合时间与raster位置。",
        "",
        "| map | D (nm) | F64−F1 (pN) | last-row−first-row median (pN) | Theil–Sen endpoint change (pN) |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in trend_rows:
        lines.append(
            f"| {row['map_order']} | {row['distance_nm']:.0f} | "
            f"{row['endpoint_delta_64_minus_1_pN']:+.2f} | "
            f"{row['last_minus_first_row_median_pN']:+.2f} | "
            f"{row['Theil_Sen_endpoint_change_pN']:+.2f} |"
        )
    lines += [
        "",
        "## Approach/retract QC",
        "",
        "| map | approach snap fraction | retract free baseline | retract snap-off fraction | censored pull-off fraction | median exact pull-off (nN) |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in map_rows:
        lines.append(
            f"| {row['map_order']} | {row['approach_snap_detected_fraction']:.1%} | "
            f"{row['retract_free_baseline_valid_fraction']:.1%} | "
            f"{row['retract_snapoff_detected_fraction']:.1%} | "
            f"{row['retract_pull_off_censored_fraction']:.1%} | "
            f"{row['retract_pull_off_force_median_nN']:.3f} |"
        )
    lines += [
        "",
        "## Calibration fingerprint",
        "",
        "| candidate | stored InvOLS reference (nm/V) | stored k reference (N/m) | log-ratio distance |",
        "|:---|---:|---:|---:|",
    ]
    for name, record in fingerprints.items():
        lines.append(
            f"| {name} | {float(record['stored_sensitivity_nm_per_V']):.3f} | "
            f"{float(record['stored_spring_constant_N_per_m']):.6f} | "
            f"{float(record['log_ratio_distance_to_new_maps']):.4f} |"
        )
    lines += [
        "",
        "## Claim boundary",
        "",
        "- Force scale is calibrated under the D5 metadata-inference assumption. If the experimental cantilever was not D5, every force must be rescaled by the correct calibrated `k/0.311340173`. Distance also changes slightly if the hard-contact consensus is redefined.",
        "- Common InvOLS follows the standing rule of calculating sensitivity from all supplied map hard contacts rather than using the JPK value. Local map medians are retained as QC, not used to rescale each map independently.",
        "- The two maps have different speeds and non-overlapping areas; neither a t-test nor their force difference identifies velocity or hydrodynamic force.",
        "- Far-linear and far-constant outputs are both retained because far-field subtraction changes the finite-distance zero-force definition.",
        "",
        "## Outputs",
        "",
        "- `map_inventory_QC.csv`, `contact_sensitivity_fits.csv`, `pixel_QC.csv`.",
        "- `map_force_curves.csv`, `pixel_force_slices.csv`, `within_map_acquisition_trends.csv`, `pixel_force_curves.npz`.",
        "- `map_distribution_comparison.csv`, `baseline_coupling.csv`.",
        "- `figures/`: F–D distributions, 20/50/100/200 nm force maps, contact/far-slope QC and 50 nm acquisition-order trends.",
        "- `provenance.json`, `artifact_manifest.sha256`.",
        "",
    ]
    return "\n".join(lines)


def render_report_three_total_maps(
    inferred_cantilever: str,
    spring_N_per_m: float,
    spring_sd_N_per_m: float,
    fingerprints: dict[str, dict[str, float | str]],
    global_sensitivity: float,
    map_rows: list[dict],
    force_rows: list[dict],
    trend_rows: list[dict],
    contact_rows: list[dict],
    comparison_rows: list[dict],
    coupling_rows: list[dict],
    paired_rows: list[dict],
    paired_metadata: dict[str, float | bool],
    paired_qc_rows: list[dict],
) -> str:
    if len(map_rows) != 3:
        raise RuntimeError("three-map report requires exactly three map summaries")
    map12_distance_um = math.hypot(
        float(map_rows[1]["x_center_um"]) - float(map_rows[0]["x_center_um"]),
        float(map_rows[1]["y_center_um"]) - float(map_rows[0]["y_center_um"]),
    )
    map23_distance_nm = math.hypot(
        float(map_rows[2]["x_center_um"]) - float(map_rows[1]["x_center_um"]),
        float(map_rows[2]["y_center_um"]) - float(map_rows[1]["y_center_um"]),
    ) * 1e3
    lines = [
        "# 31-08-26 三张 force map：raw reconstruction、速度与时间比较",
        "",
        "## 直接结论",
        "",
        f"三张文件按timestamp依次定义为map1、map2、map3。map2与map3的64个XY坐标逐点完全相同（最大坐标差 `{float(paired_metadata['maximum_pointwise_coordinate_difference_nm']):.3g} nm`），field和scan order相同，gap speed也近似相同；两次开始时间相隔 **{float(paired_metadata['start_time_lag_min']):.2f} min**。因此map2/3是本批数据中可识别时间效应的same-pixel pair。",
        "",
        f"文件没有显式cantilever名称；stored calibration fingerprint最近的是 **{inferred_cantilever}**。force暂按25.6 °C空气独立标定 `k={spring_N_per_m:.9f} N/m`（repeatability SD `{spring_sd_N_per_m:.6f} N/m`）计算。三张map全部hard-contact共同得到global liquid InvOLS **{global_sensitivity*1e9:.3f} nm/V**；文件写入的InvOLS和force conversion均未用于最终force。",
        "",
        "## Map inventory",
        "",
        "| map | start time | field | XY center (µm) | scanner/gap speed (µm/s) | local/global InvOLS (nm/V) | terminal load (nN) | far slope (pN/100 nm) |",
        "|---:|:---|:---|:---|:---|:---|---:|---:|",
    ]
    for row in map_rows:
        lines.append(
            f"| {row['map_order']} | {str(row['timestamp'])[11:19]} | "
            f"{row['field_u_um']:.1f}×{row['field_v_um']:.1f} µm | "
            f"({row['x_center_um']:.2f}, {row['y_center_um']:.2f}) | "
            f"{row['scanner_speed_median_um_per_s']:.3f} / "
            f"{row['gap_speed_20_200nm_median_um_per_s']:.3f} | "
            f"{row['map_contact_InvOLS_nm_per_V']:.3f} / "
            f"{row['global_InvOLS_used_nm_per_V']:.3f} | "
            f"{row['terminal_load_median_nN']:.3f} | "
            f"{row['far_slope_median_pN_per_100nm']:.2f} |"
        )
    lines += [
        "",
        f"map1与map2中心相距 `{map12_distance_um:.2f} µm`，速度和区域同时改变，不能用于纯velocity inference。map2与map3中心差 `{map23_distance_nm:.3g} nm`，是严格same-location repeat。",
        "",
        "## Sensitivity QC",
        "",
        "| map | 30 nm window | 40 nm window | 50 nm window | retained / 64 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for map_row in map_rows:
        selected = [row for row in contact_rows if row["source"] == map_row["source"]]
        medians = []
        for span in (30, 40, 50):
            values = np.asarray(
                [row[f"contact_invOLS_{span}nm_nm_per_V"] for row in selected],
                dtype=np.float64,
            )
            medians.append(float(np.median(values)))
        retained = sum(bool(row["retained_for_global_consensus"]) for row in selected)
        lines.append(
            f"| {map_row['map_order']} | {medians[0]:.3f} | {medians[1]:.3f} | "
            f"{medians[2]:.3f} | {retained} / 64 |"
        )
    lines += [
        "",
        "local InvOLS只作为contact-state QC；三图仍使用同一个global InvOLS，避免把真实的map/time差异通过逐图rescaling人为消掉。",
        "",
        "## Measured force slices",
        "",
        "单位pN。每项为64-pixel median；linear是逐curve远场直线扣除，constant只减远场常数。",
        "",
        "| D (nm) | map1 linear | map2 linear | map3 linear | map1 constant | map2 constant | map3 constant |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for distance_nm in TARGET_DISTANCES_NM:
        values = []
        for method in ("linear_drift_corrected", "far_constant_referenced"):
            for map_row in map_rows:
                values.append(
                    float(
                        force_value(
                            force_rows, map_row["source"], method, distance_nm
                        )["force_median_pN"]
                    )
                )
        lines.append(
            f"| {distance_nm:.0f} | "
            + " | ".join(f"{value:.2f}" for value in values)
            + " |"
        )
    lines += [
        "",
        "## Map2 → Map3 paired time effect",
        "",
        "下面的差值是每个physical pixel先做 `map3−map2` 再统计，不是两个map median相减。pixel-paired检验利用了同位置重复；考虑到8×8 raster仍有空间/时间自相关，同时给出把每个physical row压缩成一个block后的较保守p值。",
        "",
        "### Primary: far-linear baseline",
        "",
        "| D (nm) | map2 median | map3 median | paired Δ median [IQR] (pN) | map3>map2 | paired t p | Wilcoxon p | row-block p | Spearman map2/map3 |",
        "|---:|---:|---:|:---|---:|---:|---:|---:|---:|",
    ]
    primary_paired = [
        row for row in paired_rows if row["baseline_method"] == "linear_drift_corrected"
    ]
    for row in primary_paired:
        lines.append(
            f"| {row['distance_nm']:.0f} | {row['map2_median_pN']:.2f} | "
            f"{row['map3_median_pN']:.2f} | "
            f"{row['paired_difference_map3_minus_map2_median_pN']:+.2f} "
            f"[{row['paired_difference_q25_pN']:+.2f}, {row['paired_difference_q75_pN']:+.2f}] | "
            f"{row['fraction_map3_greater_map2']:.1%} | "
            f"{row['paired_t_p_value']:.3g} | {row['Wilcoxon_p_value']:.3g} | "
            f"{row['row_block_Wilcoxon_p_value']:.3g} | "
            f"{row['Spearman_rho_map2_map3']:+.3f} |"
        )
    lines += [
        "",
        "### Baseline sensitivity of the paired time effect",
        "",
        "| D (nm) | paired Δ median, linear (pN) | paired Δ median, constant (pN) | constant Wilcoxon p |",
        "|---:|---:|---:|---:|",
    ]
    constant_paired = [
        row for row in paired_rows if row["baseline_method"] == "far_constant_referenced"
    ]
    for linear_row, constant_row in zip(primary_paired, constant_paired, strict=True):
        lines.append(
            f"| {linear_row['distance_nm']:.0f} | "
            f"{linear_row['paired_difference_map3_minus_map2_median_pN']:+.2f} | "
            f"{constant_row['paired_difference_map3_minus_map2_median_pN']:+.2f} | "
            f"{constant_row['Wilcoxon_p_value']:.3g} |"
        )
    qc_by_metric = {row["metric"]: row for row in paired_qc_rows}
    speed_delta = float(
        qc_by_metric["gap_speed"]["paired_difference_map3_minus_map2_median"]
    )
    speed_relative = speed_delta / float(qc_by_metric["gap_speed"]["map2_median"])
    water_viscosity_Pa_s = (
        base.cheng_viscosity_mPa_s(0.0, CONTEXT_TEMPERATURE_C) * 1e-3
    )
    hydrodynamic_speed_increment_pN = {
        distance_nm: 6.0
        * math.pi
        * water_viscosity_Pa_s
        * CONTEXT_PROBE_RADIUS_M**2
        * (speed_delta * 1e-6)
        / (distance_nm * 1e-9)
        * 1e12
        for distance_nm in (20.0, 50.0)
    }
    sensitivity_delta = float(
        qc_by_metric["local_contact_InvOLS"][
            "paired_difference_map3_minus_map2_median"
        ]
    )
    pull_off_delta = float(
        qc_by_metric["retract_pull_off_force"][
            "paired_difference_map3_minus_map2_median"
        ]
    )
    lines += [
        "",
        "若linear与constant给出同方向、相近量级，time effect对zero-force定义较稳健；若两者差别大，则变化主要与far-field slope/baseline有关，而不能解释为真实interaction force变化。Spearman相关衡量两次是否保留相同spatial pattern：高相关更像原空间图样上的幅值/offset变化，低相关则表示图样本身重排。",
        "",
        "## Time-effect interpretation",
        "",
        f"- map2/3 paired gap-speed median只改变 `{speed_delta:+.4f} µm/s`（`{speed_relative:+.2%}`），而terminal load、local contact InvOLS与far-slope median也基本稳定。按25.6 °C纯水、no-slip sphere-plane估算，这个speed差在20/50 nm只对应 `{hydrodynamic_speed_increment_pN[20.0]:.3f}/{hydrodynamic_speed_increment_pN[50.0]:.3f} pN`。local InvOLS的paired median变化仅 `{sensitivity_delta:+.3f} nm/V`，Wilcoxon `p={float(qc_by_metric['local_contact_InvOLS']['Wilcoxon_p_value']):.3g}`。因此20–100 nm的几十pN变化不能合理归因于speed或重新计算sensitivity。",
        f"- 20、50、100 nm的paired median均随时间下降，linear/constant分别约 `−32/−35`、`−28/−31`、`−12/−11 pN`；方向和量级对baseline定义相对稳定。200 nm则转为正值且baseline敏感，主要反映far-field zero/slope变化。",
        "- 这些下降不是uniform offset：20 nm的map2/map3 spatial Spearman仍较高，但50–200 nm只中低相关；每个pixel的difference IQR很宽并跨过零。50 nm两个map的整体median相差约−73 pN，但same-pixel paired median只有约−28 pN，说明distribution shape和spatial pattern也随时间改变。",
        "- naive paired-pixel test在50 nm处接近显著或显著，但8个physical-row block的Wilcoxon结果不显著；因此数据表明测得force distribution在27 min后并不完全可重复，支持存在time/history/repeatability effect，尚不支持一个全区域统一、精确为−30 pN的time correction。",
        f"- retract pull-off的paired median变化 `{pull_off_delta:+.3f} nN`，但Wilcoxon `p={float(qc_by_metric['retract_pull_off_force']['Wilcoxon_p_value']):.3g}`；当前没有证据表明adhesion发生统一变化。",
        "",
        "### Paired QC variables",
        "",
        "| metric | map2 median | map3 median | paired Δ median [IQR] | Wilcoxon p | Spearman |",
        "|:---|---:|---:|:---|---:|---:|",
    ]
    for row in paired_qc_rows:
        lines.append(
            f"| {row['metric']} ({row['unit']}) | {row['map2_median']:.4g} | "
            f"{row['map3_median']:.4g} | "
            f"{row['paired_difference_map3_minus_map2_median']:+.4g} "
            f"[{row['paired_difference_q25']:+.4g}, {row['paired_difference_q75']:+.4g}] | "
            f"{row['Wilcoxon_p_value']:.3g} | "
            f"{row['Spearman_rho_map2_map3']:+.3f} |"
        )
    lines += [
        "",
        "## Within-map acquisition-order trend",
        "",
        "Theil–Sen endpoint change是该map内部64个取样序号的robust斜率乘63；row与时间仍然部分混杂。",
        "",
        "| map | D (nm) | F64−F1 (pN) | row8−row1 median (pN) | Theil–Sen endpoint (pN) |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in trend_rows:
        lines.append(
            f"| {row['map_order']} | {row['distance_nm']:.0f} | "
            f"{row['endpoint_delta_64_minus_1_pN']:+.2f} | "
            f"{row['last_minus_first_row_median_pN']:+.2f} | "
            f"{row['Theil_Sen_endpoint_change_pN']:+.2f} |"
        )
    lines += [
        "",
        "## Approach/retract contact-state QC",
        "",
        "| map | approach snap fraction | retract snap-off fraction | pull-off median (nN) | detachment travel median (nm) |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in map_rows:
        lines.append(
            f"| {row['map_order']} | {row['approach_snap_detected_fraction']:.1%} | "
            f"{row['retract_snapoff_detected_fraction']:.1%} | "
            f"{row['retract_pull_off_force_median_nN']:.3f} | "
            f"{row['retract_detachment_piezo_travel_median_nm']:.2f} |"
        )
    lines += [
        "",
        "snap detection是algorithmic QC flag；特别是异常长距离检测不能直接解释成physical snap-in position。pull-off和contact InvOLS若在map2/3间同步变化，则支持probe/contact state随时间改变。",
        "",
        "## Calibration fingerprint",
        "",
        "| candidate | stored InvOLS reference (nm/V) | stored k reference (N/m) | log-ratio distance |",
        "|:---|---:|---:|---:|",
    ]
    for name, record in fingerprints.items():
        lines.append(
            f"| {name} | {float(record['stored_sensitivity_nm_per_V']):.3f} | "
            f"{float(record['stored_spring_constant_N_per_m']):.6f} | "
            f"{float(record['log_ratio_distance_to_new_maps']):.4f} |"
        )
    lines += [
        "",
        "## Claim boundary",
        "",
        "- map2/map3支持same-location、same-speed条件下的约27 min时间差比较；它仍是一次map pair，不能单独建立通用time law。",
        "- paired pixel p值仍受spatial/serial correlation影响；row-block结果是更保守的敏感性检查，而不是完整的spatial model。",
        "- Force scale以D5 metadata-inference为条件；若实际cantilever不是D5，应按正确calibrated k整体rescale force。",
        "- far-linear与far-constant都保留，因为zero-force定义会显著改变100–200 nm结果。",
        "- map1/map2同时改变区域和速度，不能把差值唯一归因于hydrodynamics。",
        "",
        "## Outputs",
        "",
        "- `map_inventory_QC.csv`, `contact_sensitivity_fits.csv`, `pixel_QC.csv`.",
        "- `map_force_curves.csv`, `pixel_force_slices.csv`, `within_map_acquisition_trends.csv`, `pixel_force_curves.npz`.",
        "- `map2_map3_paired_time_comparison.csv`, `map2_map3_paired_QC.csv`, `map_distribution_comparison.csv`, `baseline_coupling.csv`.",
        "- `figures/`: 三图F–D/map/QC，map2/3 paired force maps与time statistics.",
        "- `provenance.json`, `artifact_manifest.sha256`.",
        "",
    ]
    return "\n".join(lines)


def render_report(
    inferred_cantilever: str,
    spring_N_per_m: float,
    spring_sd_N_per_m: float,
    fingerprints: dict[str, dict[str, float | str]],
    global_sensitivity: float,
    map_rows: list[dict],
    force_rows: list[dict],
    within_map_trend_rows: list[dict],
    contact_rows: list[dict],
    time_pair_rows: list[dict],
    time_trend_rows: list[dict],
    time_metadata: dict,
    paired_qc_rows: list[dict],
    acquisition_fit_rows: list[dict],
) -> str:
    if len(map_rows) != 5:
        raise RuntimeError("five-map report requires five map summaries")
    relative_times = [float(value) for value in time_metadata["relative_start_times_min"]]
    relative_time_text = "、".join(f"{value:.2f}" for value in relative_times)
    lines = [
        "# 31-08-26 五张 force map：map2–5 同地点四时点比较",
        "",
        "## 直接结论",
        "",
        f"全部文件按timestamp依次定义为map1–map5。map2–map5的64个XY坐标逐点完全相同（最大差 `{float(time_metadata['maximum_pointwise_coordinate_difference_nm']):.3g} nm`），相对开始时间为 **{relative_time_text} min**。四图field、scan order和nominal speed相同，因此构成same-location四时点时间序列。",
        "",
        f"Force scale暂按metadata fingerprint最近的 **{inferred_cantilever}**、25.6 °C空气独立标定 `k={spring_N_per_m:.9f} N/m`（repeatability SD `{spring_sd_N_per_m:.6f} N/m`）。五张map全部hard-contact共同得到global liquid InvOLS **{global_sensitivity*1e9:.3f} nm/V**；文件写入的InvOLS与force conversion未用于最终force。",
        "",
        "## Map inventory",
        "",
        "| map | start time | relative time (min) | field | XY center (µm) | scanner/gap speed (µm/s) | local/global InvOLS (nm/V) | terminal load (nN) | far slope (pN/100 nm) |",
        "|---:|:---|---:|:---|:---|:---|:---|---:|---:|",
    ]
    for row in map_rows:
        relative = "—" if int(row["map_order"]) == 1 else f"{relative_times[int(row['map_order']) - 2]:.2f}"
        lines.append(
            f"| {row['map_order']} | {str(row['timestamp'])[11:19]} | {relative} | "
            f"{row['field_u_um']:.1f}×{row['field_v_um']:.1f} µm | "
            f"({row['x_center_um']:.2f}, {row['y_center_um']:.2f}) | "
            f"{row['scanner_speed_median_um_per_s']:.3f} / "
            f"{row['gap_speed_20_200nm_median_um_per_s']:.3f} | "
            f"{row['map_contact_InvOLS_nm_per_V']:.3f} / "
            f"{row['global_InvOLS_used_nm_per_V']:.3f} | "
            f"{row['terminal_load_median_nN']:.3f} | "
            f"{row['far_slope_median_pN_per_100nm']:.2f} |"
        )
    lines += [
        "",
        "map1与后四图位置、field和速度不同，只保留为背景；时间结论只来自map2–5。",
        "",
        "## Sensitivity QC",
        "",
        "| map | 30 nm contact | 40 nm contact | 50 nm contact | retained / 64 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for map_row in map_rows:
        selected = [row for row in contact_rows if row["source"] == map_row["source"]]
        medians = [
            float(
                np.median(
                    np.asarray(
                        [row[f"contact_invOLS_{span}nm_nm_per_V"] for row in selected]
                    )
                )
            )
            for span in (30, 40, 50)
        ]
        retained = sum(bool(row["retained_for_global_consensus"]) for row in selected)
        lines.append(
            f"| {map_row['map_order']} | {medians[0]:.3f} | {medians[1]:.3f} | "
            f"{medians[2]:.3f} | {retained} / 64 |"
        )
    lines += [
        "",
        "## Measured force slices",
        "",
        "以下为far-linear、64-pixel median，单位pN。",
        "",
        "| D (nm) | map1 | map2 | map3 | map4 | map5 |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for distance_nm in TARGET_DISTANCES_NM:
        values = [
            float(
                force_value(
                    force_rows,
                    row["source"],
                    "linear_drift_corrected",
                    distance_nm,
                )["force_median_pN"]
            )
            for row in map_rows
        ]
        lines.append(
            f"| {distance_nm:.0f} | " + " | ".join(f"{value:.2f}" for value in values) + " |"
        )
    lines += [
        "",
        "## 四时点pixel-level时间趋势",
        "",
        f"对每个physical pixel分别用map2–5四个点拟合force对时间的直线，再把slope乘总时间 `{relative_times[-1]:.2f} min`。这比直接拟合四个map median更能保留same-pixel结构；但每个pixel只有四个时间点，仍只能视为linear trend diagnostic。",
        "",
        "### Primary: far-linear baseline",
        "",
        "| D (nm) | map2 → map3 → map4 → map5 medians (pN) | fitted Map5−Map2 median [IQR] (pN) | negative slopes | strict decrease | strict increase | pixel p | row-block p |",
        "|---:|:---|:---|---:|---:|---:|---:|---:|",
    ]
    primary_trends = [
        row for row in time_trend_rows if row["baseline_method"] == "linear_drift_corrected"
    ]
    constant_trends = [
        row for row in time_trend_rows if row["baseline_method"] == "far_constant_referenced"
    ]
    for row in primary_trends:
        lines.append(
            f"| {row['distance_nm']:.0f} | {row['map2_median_pN']:.2f} → "
            f"{row['map3_median_pN']:.2f} → {row['map4_median_pN']:.2f} → "
            f"{row['map5_median_pN']:.2f} | "
            f"{row['fitted_full_interval_change_median_pN']:+.2f} "
            f"[{row['fitted_full_interval_change_q25_pN']:+.2f}, "
            f"{row['fitted_full_interval_change_q75_pN']:+.2f}] | "
            f"{row['fraction_negative_time_slope']:.1%} | "
            f"{row['fraction_strict_monotonic_decrease']:.1%} | "
            f"{row['fraction_strict_monotonic_increase']:.1%} | "
            f"{row['pixel_slope_Wilcoxon_p_value']:.3g} | "
            f"{row['row_block_slope_Wilcoxon_p_value']:.3g} |"
        )
    lines += [
        "",
        "### Baseline sensitivity of full-interval trend",
        "",
        "| D (nm) | fitted change linear (pN) | fitted change constant (pN) | constant row-block p |",
        "|---:|---:|---:|---:|",
    ]
    for linear_row, constant_row in zip(primary_trends, constant_trends, strict=True):
        lines.append(
            f"| {linear_row['distance_nm']:.0f} | "
            f"{linear_row['fitted_full_interval_change_median_pN']:+.2f} | "
            f"{constant_row['fitted_full_interval_change_median_pN']:+.2f} | "
            f"{constant_row['row_block_slope_Wilcoxon_p_value']:.3g} |"
        )
    lines += [
        "",
        "## 相邻与全间隔paired changes",
        "",
        "| pair | D (nm) | paired Δ median [IQR] (pN) | Wilcoxon p | row-block p | spatial Spearman |",
        "|:---|---:|:---|---:|---:|---:|",
    ]
    primary_pairs = [
        row for row in time_pair_rows if row["baseline_method"] == "linear_drift_corrected"
    ]
    for row in primary_pairs:
        lines.append(
            f"| map{row['before_map']}→map{row['after_map']} | {row['distance_nm']:.0f} | "
            f"{row['paired_difference_after_minus_before_median_pN']:+.2f} "
            f"[{row['paired_difference_q25_pN']:+.2f}, {row['paired_difference_q75_pN']:+.2f}] | "
            f"{row['Wilcoxon_p_value']:.3g} | {row['row_block_Wilcoxon_p_value']:.3g} | "
            f"{row['Spearman_rho_before_after']:+.3f} |"
        )
    primary_by_distance = {
        float(row["distance_nm"]): row for row in primary_trends
    }
    constant_by_distance = {
        float(row["distance_nm"]): row for row in constant_trends
    }
    pair_lookup = {
        (
            int(row["before_map"]),
            int(row["after_map"]),
            float(row["distance_nm"]),
        ): row
        for row in primary_pairs
    }
    qc_full = {
        row["metric"]: row
        for row in paired_qc_rows
        if int(row["before_map"]) == 2 and int(row["after_map"]) == 5
    }
    speed_delta_um_per_s = float(
        qc_full["gap_speed"]["paired_difference_after_minus_before_median"]
    )
    water_viscosity_Pa_s = (
        base.cheng_viscosity_mPa_s(0.0, CONTEXT_TEMPERATURE_C) * 1e-3
    )
    hydrodynamic_increment = {
        distance_nm: 6.0
        * math.pi
        * water_viscosity_Pa_s
        * CONTEXT_PROBE_RADIUS_M**2
        * (speed_delta_um_per_s * 1e-6)
        / (distance_nm * 1e-9)
        * 1e12
        for distance_nm in (20.0, 50.0)
    }
    lines += [
        "",
        "## Interpretation of the four-map sequence",
        "",
        f"- **20 nm存在最清楚的time/history dependence。** map median为 `{primary_by_distance[20.0]['map2_median_pN']:.1f} → {primary_by_distance[20.0]['map3_median_pN']:.1f} → {primary_by_distance[20.0]['map4_median_pN']:.1f} → {primary_by_distance[20.0]['map5_median_pN']:.1f} pN`；same-pixel full-interval change为 `{primary_by_distance[20.0]['fitted_full_interval_change_median_pN']:+.1f} pN`（linear）或 `{constant_by_distance[20.0]['fitted_full_interval_change_median_pN']:+.1f} pN`（constant）。linear和constant的row-block p分别为 `{primary_by_distance[20.0]['row_block_slope_Wilcoxon_p_value']:.3g}`、`{constant_by_distance[20.0]['row_block_slope_Wilcoxon_p_value']:.3g}`。",
        f"- **50 nm的精确量级更受baseline与空间相关影响。** full-interval paired trend为 `{primary_by_distance[50.0]['fitted_full_interval_change_median_pN']:+.1f} pN`（linear）和 `{constant_by_distance[50.0]['fitted_full_interval_change_median_pN']:+.1f} pN`（constant）；negative slope pixel占 `{primary_by_distance[50.0]['fraction_negative_time_slope']:.1%}`。linear row-block p为 `{primary_by_distance[50.0]['row_block_slope_Wilcoxon_p_value']:.3g}`，constant为 `{constant_by_distance[50.0]['row_block_slope_Wilcoxon_p_value']:.3g}`。因此不能把任一单一拟合值当作所有pixel的唯一校正。",
        f"- **100 nm结论依赖zero-force定义。** linear给出 `{primary_by_distance[100.0]['fitted_full_interval_change_median_pN']:+.1f} pN`且row-block不显著；constant给出 `{constant_by_distance[100.0]['fitted_full_interval_change_median_pN']:+.1f} pN`且row-block p=`{constant_by_distance[100.0]['row_block_slope_Wilcoxon_p_value']:.3g}`。200 nm主要受far-field baseline控制，不宜作为absolute interaction-force time trend。",
        f"- **这不是简单统一的linear/exponential relaxation。** 相邻map在20 nm的paired median依次变化 `{pair_lookup[(2, 3, 20.0)]['paired_difference_after_minus_before_median_pN']:+.1f}`、`{pair_lookup[(3, 4, 20.0)]['paired_difference_after_minus_before_median_pN']:+.1f}`、`{pair_lookup[(4, 5, 20.0)]['paired_difference_after_minus_before_median_pN']:+.1f} pN`；50 nm依次为 `{pair_lookup[(2, 3, 50.0)]['paired_difference_after_minus_before_median_pN']:+.1f}`、`{pair_lookup[(3, 4, 50.0)]['paired_difference_after_minus_before_median_pN']:+.1f}`、`{pair_lookup[(4, 5, 50.0)]['paired_difference_after_minus_before_median_pN']:+.1f} pN`。四个时点仍不足以唯一确定time law。",
        f"- speed、load和sensitivity不是主要解释。map2→5 paired gap-speed只改变 `{speed_delta_um_per_s:+.4f} µm/s`；按纯水no-slip估算，在20/50 nm只对应 `{hydrodynamic_increment[20.0]:.3f}/{hydrodynamic_increment[50.0]:.3f} pN`。local InvOLS paired median变化 `{float(qc_full['local_contact_InvOLS']['paired_difference_after_minus_before_median']):+.3f} nm/V`，terminal load变化 `{float(qc_full['terminal_load']['paired_difference_after_minus_before_median']):+.4f} nN`。",
        f"- 最稳妥的结论是：20–50 nm force随约{relative_times[-1]:.1f} min实验history变化，但具有空间异质性，不能用单一offset或单一time slope校正所有pixel。",
        "",
    ]
    lines += [
        "",
        "## Within-map acquisition-order trend",
        "",
        "| map | D (nm) | Theil–Sen endpoint change (pN) |",
        "|---:|---:|---:|",
    ]
    for row in within_map_trend_rows:
        lines.append(
            f"| {row['map_order']} | {row['distance_nm']:.0f} | "
            f"{row['Theil_Sen_endpoint_change_pN']:+.2f} |"
        )
    lines += [
        "",
        "### Map2–5 at 50 nm: separate 64-point fits",
        "",
        "Theil–Sen为主拟合，OLS只作对照。斜率的sample单位是一条force curve；1→64拟合变化为斜率×63。它描述acquisition-order association，不能单独解释为纯时间drift，因为raster位置与时间混杂。",
        "",
        "| map | Theil–Sen slope (pN/sample) | 95% slope CI | fitted 1→64 change (pN) | OLS change (pN) | OLS R² |",
        "|---:|---:|:---|---:|---:|---:|",
    ]
    for row in acquisition_fit_rows:
        lines.append(
            f"| {row['map_order']} | {row['Theil_Sen_slope_pN_per_sample']:+.3f} | "
            f"[{row['Theil_Sen_slope_95CI_low_pN_per_sample']:+.3f}, "
            f"{row['Theil_Sen_slope_95CI_high_pN_per_sample']:+.3f}] | "
            f"{row['Theil_Sen_endpoint_change_pN']:+.2f} | "
            f"{row['OLS_endpoint_change_pN']:+.2f} | "
            f"{row['OLS_r_squared']:.3f} |"
        )
    lines += [
        "",
        "## QC pair comparisons",
        "",
        "以下只列map2→map5全间隔QC；完整六对结果在CSV中。",
        "",
        "| metric | map2 | map5 | paired Δ median [IQR] | Wilcoxon p |",
        "|:---|---:|---:|:---|---:|",
    ]
    for row in paired_qc_rows:
        if int(row["before_map"]) != 2 or int(row["after_map"]) != 5:
            continue
        lines.append(
            f"| {row['metric']} ({row['unit']}) | {row['before_median']:.4g} | "
            f"{row['after_median']:.4g} | "
            f"{row['paired_difference_after_minus_before_median']:+.4g} "
            f"[{row['paired_difference_q25']:+.4g}, {row['paired_difference_q75']:+.4g}] | "
            f"{row['Wilcoxon_p_value']:.3g} |"
        )
    lines += [
        "",
        "## Calibration fingerprint",
        "",
        "| candidate | stored InvOLS reference (nm/V) | stored k reference (N/m) | log-ratio distance |",
        "|:---|---:|---:|---:|",
    ]
    for name, record in fingerprints.items():
        lines.append(
            f"| {name} | {float(record['stored_sensitivity_nm_per_V']):.3f} | "
            f"{float(record['stored_spring_constant_N_per_m']):.6f} | "
            f"{float(record['log_ratio_distance_to_new_maps']):.4f} |"
        )
    lines += [
        "",
        "## Claim boundary",
        "",
        "- map2–5是same-location、same-speed四时点序列，可判断repeatability/chronology是否存在，但四个时点仍不足以唯一确定linear、exponential或其他time law。",
        "- 64个pixel具有spatial/serial correlation；pixel p值偏乐观，8个physical-row block结果作为保守敏感性检查。",
        "- far-linear与far-constant都保留；两者不一致时不能把变化解释为absolute interaction-force变化。",
        "- Force scale仍以D5 metadata inference为条件；若实际cantilever不同，force需按正确calibrated k整体rescale。",
        "- map1不参与时间拟合，也不能与map2唯一识别velocity effect。",
        "",
        "## Outputs",
        "",
        "- `map2_map5_pairwise_time_comparison.csv`, `map2_map5_time_trends.csv`, `map2_map5_paired_QC.csv`.",
        "- `map_inventory_QC.csv`, `contact_sensitivity_fits.csv`, `pixel_QC.csv`, `map_force_curves.csv`, `pixel_force_curves.npz`.",
        "- `map2_map5_force_50nm_acquisition_fits.csv` and the four separate `figures/map[2-5]_force_50nm_vs_acquisition_order.png` plots.",
        "- `map2_map5_absolute_force_time_statistics.csv` and `figures/absolute_force_vs_time_mean_std_variance.png`; statistics use per-pixel `|F|`, with sample SD/variance (`ddof=1`).",
        "- `figures/map2_map5_time_force_maps.png`, `figures/map2_map5_time_statistics.png`及五图总览/QC。",
        "- `provenance.json`, `artifact_manifest.sha256`.",
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
    paths = sorted(DATA_ROOT.glob("*.jpk-force-map"))
    if len(paths) != 5:
        raise RuntimeError(f"expected exactly five 31-08-26 maps, got {len(paths)}")
    sources = [base.load_source(path.resolve(), 0) for path in paths]
    sources.sort(key=lambda source: source.timestamp)

    inferred, spring, spring_sd, fingerprints = infer_cantilever(sources)
    contact_rows, global_sensitivity = contact_analysis(sources)
    pixel_rows, map_rows, force_rows, trend_rows, matrices = analyze_maps(
        sources, spring
    )
    comparison_rows = distribution_comparisons(sources, matrices)
    coupling_rows = baseline_coupling(sources, pixel_rows, matrices)
    time_pair_rows, time_trend_rows, time_metadata = same_location_time_comparisons(
        sources, matrices
    )
    paired_qc_rows = same_location_paired_qc_comparisons(pixel_rows)
    acquisition_fit_rows = plot_individual_acquisition_order_50nm(
        sources, map_rows, matrices
    )
    absolute_force_time_rows = plot_absolute_force_time_statistics(
        sources, matrices
    )

    write_csv(RESULTS / "contact_sensitivity_fits.csv", contact_rows)
    write_csv(RESULTS / "pixel_QC.csv", pixel_rows)
    write_csv(RESULTS / "map_inventory_QC.csv", map_rows)
    write_csv(RESULTS / "map_force_curves.csv", force_rows)
    write_csv(RESULTS / "within_map_acquisition_trends.csv", trend_rows)
    write_csv(RESULTS / "map_distribution_comparison.csv", comparison_rows)
    write_csv(RESULTS / "baseline_coupling.csv", coupling_rows)
    write_csv(
        RESULTS / "map2_map5_pairwise_time_comparison.csv", time_pair_rows
    )
    write_csv(RESULTS / "map2_map5_time_trends.csv", time_trend_rows)
    write_csv(RESULTS / "map2_map5_paired_QC.csv", paired_qc_rows)
    write_csv(
        RESULTS / "map2_map5_force_50nm_acquisition_fits.csv",
        acquisition_fit_rows,
    )
    write_csv(
        RESULTS / "map2_map5_absolute_force_time_statistics.csv",
        absolute_force_time_rows,
    )

    pixel_slice_rows: list[dict] = []
    for source, summary in zip(sources, map_rows, strict=True):
        key = relative_source(source.path)
        for target in TARGET_DISTANCES_NM:
            index = int(np.flatnonzero(np.isclose(BIN_CENTERS_NM, target))[0])
            for point in range(64):
                row, column = base.map_pixel_from_index(point, 8, True)
                pixel_slice_rows.append(
                    {
                        "source": key,
                        "map_order": summary["map_order"],
                        "point_index": point,
                        "row": row,
                        "column": column,
                        "distance_nm": target,
                        "force_linear_drift_corrected_pN": matrices[
                            key + "|line"
                        ][point, index],
                        "force_far_constant_referenced_pN": matrices[
                            key + "|constant"
                        ][point, index],
                    }
                )
    write_csv(RESULTS / "pixel_force_slices.csv", pixel_slice_rows)
    np.savez_compressed(
        RESULTS / "pixel_force_curves.npz",
        distance_nm=BIN_CENTERS_NM,
        source=np.asarray([relative_source(source.path) for source in sources]),
        force_linear_drift_corrected_pN=np.asarray(
            [matrices[relative_source(source.path) + "|line"] for source in sources]
        ),
        force_far_constant_referenced_pN=np.asarray(
            [
                matrices[relative_source(source.path) + "|constant"]
                for source in sources
            ]
        ),
    )

    plot_force_distance(sources, map_rows, matrices)
    plot_force_maps(sources, map_rows, matrices)
    plot_baseline_comparison(sources, map_rows, matrices)
    plot_qc_maps(pixel_rows, map_rows)
    plot_acquisition_order(sources, map_rows, matrices)
    plot_same_location_time_maps(sources, matrices)
    plot_same_location_time_statistics(sources, matrices, time_trend_rows)

    report = render_report(
        inferred,
        spring,
        spring_sd,
        fingerprints,
        global_sensitivity,
        map_rows,
        force_rows,
        trend_rows,
        contact_rows,
        time_pair_rows,
        time_trend_rows,
        time_metadata,
        paired_qc_rows,
        acquisition_fit_rows,
    )
    (RESULTS / "REPORT.md").write_text(report, encoding="utf-8")

    provenance = {
        "analysis": "31-08-26 five-map raw reconstruction, QC, and same-location map2/map3/map4/map5 time sequence",
        "claim_status": "calibrated_force_conditional_on_D5_metadata_identity_inference",
        "archive": {
            "name": ARCHIVE_NAME,
            "sha256_before_deletion": ARCHIVE_SHA256_BEFORE_DELETION,
            "zip_crc_test": "passed before extraction",
            "deleted_after_verified_extraction": not (ROOT / ARCHIVE_NAME).exists(),
        },
        "input_hashes": {
            relative_source(source.path): source.sha256 for source in sources
        },
        "map_environment": "not encoded in source map metadata",
        "cantilever_name_encoded_in_maps": False,
        "cantilever_identity_inference": {
            "selected": inferred,
            "method": "minimum Euclidean distance between log ratios of stored InvOLS and stored spring constant versus D4-D6 calibration-file fingerprints",
            "fingerprints": fingerprints,
        },
        "independent_calibration_summary": relative_source(CALIBRATION_SUMMARY),
        "independent_calibrated_spring_constant_N_per_m": spring,
        "spring_repeatability_sd_N_per_m": spring_sd,
        "contact_primary_span_nm": CONTACT_SPAN_NM,
        "contact_check_spans_nm": CONTACT_CHECK_SPANS_NM,
        "contact_R2_min": CONTACT_R2_MIN,
        "global_map_hard_contact_InvOLS_nm_per_V": global_sensitivity * 1e9,
        "force_definition": "F=k*InvOLS*(vDeflection - per-curve far baseline)",
        "distance_definition": "D=measuredHeight + InvOLS*corrected_vDeflection - fitted_contact_height",
        "force_baselines": [
            "per-curve far-field linear subtraction",
            "per-curve far-field constant reference",
        ],
        "distribution_comparison": {
            "tests": ["Welch t-test", "Mann-Whitney U"],
            "effect_sizes": ["Cliff's delta", "Cohen's d"],
            "claim_limit": "descriptive naive-pixel comparison; spatial correlation and non-overlapping map areas preclude a velocity-causal test",
        },
        "baseline_coupling": "Spearman correlation of force with far-slope, local contact InvOLS, acquisition order, physical row, and physical column",
        "paired_time_comparison": {
            **time_metadata,
            "time_series": "map2, map3, map4, map5",
            "tests": [
                "pairwise paired t-test and Wilcoxon signed-rank over matched pixels",
                "pairwise Wilcoxon signed-rank over eight physical-row median differences",
                "per-pixel four-timepoint linear slope with pixel and row-block Wilcoxon tests",
            ],
            "claim_limit": "four same-location time points; time-law functional form remains weakly identifiable and spatial/serial correlation remains",
        },
        "paired_QC_comparison": {
            "pairs": [
                "map2 versus map3",
                "map2 versus map4",
                "map2 versus map5",
                "map3 versus map4",
                "map3 versus map5",
                "map4 versus map5",
            ],
            "test": "Wilcoxon signed-rank over matched pixels",
            "metrics": sorted(set(row["metric"] for row in paired_qc_rows)),
        },
        "within_map_50nm_linear_fits": {
            "maps": [2, 3, 4, 5],
            "primary_fit": "Theil-Sen robust slope with scipy.stats.theilslopes 95% slope confidence interval",
            "comparison_fit": "ordinary least squares via scipy.stats.linregress",
            "sample_order": "1 through 64 in stored map acquisition order",
            "claim_limit": "acquisition order confounds elapsed time with raster position",
        },
        "absolute_force_time_statistics": {
            "maps": [2, 3, 4, 5],
            "distances_nm": list(TARGET_DISTANCES_NM),
            "force_definition": "absolute value of each pixel force before aggregation",
            "center": "arithmetic mean over 64 pixels",
            "spread": "sample standard deviation and sample variance over 64 pixels with ddof=1",
            "baseline_method": "per-curve far-field linear subtraction",
            "claim_limit": "spatial distribution width, not uncertainty of the mean; absolute value folds positive and negative far-field noise together",
        },
        "conditional_hydrodynamic_scale_check": {
            "temperature_C_from_user_context": CONTEXT_TEMPERATURE_C,
            "liquid_assumption": "pure water",
            "probe_radius_m_from_prior_sphere_characterization": CONTEXT_PROBE_RADIUS_M,
            "boundary_condition": "no slip",
            "formula": "delta_F=6*pi*eta*R^2*delta_v/D",
            "claim_limit": "order-of-magnitude comparison only; liquid and probe are not encoded in the map files",
        },
        "bin_centers_nm": BIN_CENTERS_NM.tolist(),
        "target_distances_nm": TARGET_DISTANCES_NM,
        "randomness": "none",
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "dependency_hashes": {
            relative_source(ROOT / "analysis" / "fit_glycerol_surface_forces.py"): sha256_file(
                ROOT / "analysis" / "fit_glycerol_surface_forces.py"
            ),
            relative_source(ROOT / "analysis" / "analyze_velocity_systematics.py"): sha256_file(
                ROOT / "analysis" / "analyze_velocity_systematics.py"
            ),
            relative_source(
                ROOT / "analysis" / "analyze_27_08_26_palindrome_pilot.py"
            ): sha256_file(
                ROOT / "analysis" / "analyze_27_08_26_palindrome_pilot.py"
            ),
        },
    }
    (RESULTS / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    for stale_relative_path in (
        "map2_map3_paired_time_comparison.csv",
        "map2_map3_paired_QC.csv",
        "map2_map3_map4_pairwise_time_comparison.csv",
        "map2_map3_map4_time_trends.csv",
        "map2_map3_map4_paired_QC.csv",
        "map2_map3_map4_force_50nm_acquisition_fits.csv",
        "figures/map2_map3_paired_force_maps.png",
        "figures/map2_map3_paired_time_statistics.png",
        "figures/map2_map3_map4_time_force_maps.png",
        "figures/map2_map3_map4_time_statistics.png",
    ):
        (RESULTS / stale_relative_path).unlink(missing_ok=True)
    artifacts = [
        path
        for path in RESULTS.rglob("*")
        if path.is_file() and path.name != "artifact_manifest.sha256"
    ]
    artifacts += [Path(__file__).resolve(), CALIBRATION_SUMMARY]
    create_manifest(artifacts, RESULTS / "artifact_manifest.sha256")

    print(f"Wrote {RESULTS}")
    print(
        f"cantilever={inferred}, k={spring:.9f} N/m, "
        f"global InvOLS={global_sensitivity*1e9:.6f} nm/V"
    )
    for row in map_rows:
        print(
            f"map{row['map_order']}: scanner={row['scanner_speed_median_um_per_s']:.4f} "
            f"um/s, gap={row['gap_speed_20_200nm_median_um_per_s']:.4f} um/s, "
            f"local InvOLS={row['map_contact_InvOLS_nm_per_V']:.3f} nm/V"
        )


if __name__ == "__main__":
    main()

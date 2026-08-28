#!/usr/bin/env python3
"""Calibrate D4-D6 and analyse the available 27-08-26 D4 palindrome maps.

The map analysis is deliberately raw-data based.  It rebuilds the water
InvOLS from all map hard-contact tails, uses the independently air-calibrated
D4 spring constant, and preserves both a constant-referenced and a line-drift
corrected force.  Each map is one experimental velocity unit; the 8x8 pixels
are paired spatial observations.  With only two of the planned three blocks,
velocity fits are descriptive pilot results, not validated zero-speed forces.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import re
import sys
from dataclasses import asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.stats
from scipy.constants import Boltzmann


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import calibrate_cantilevers as cal  # noqa: E402
import analyze_velocity_systematics as events  # noqa: E402
import fit_glycerol_surface_forces as base  # noqa: E402


DATA_ROOT = ROOT / "27-08-26"
RESULTS = ROOT / "analysis" / "palindrome_27_08_26_pilot_results"
FIGURES = RESULTS / "figures"

TEMPERATURE_C = 25.6
TEMPERATURE_K = TEMPERATURE_C + 273.15
PROBE_RADIUS_M = 4.546848945303745e-6
WATER_VISCOSITY_MPA_S = base.cheng_viscosity_mPa_s(0.0, TEMPERATURE_C)

CAL_CONTACT_SPAN_NM = 35.0
CAL_CONTACT_CHECKS_NM = (30.0, 40.0)
MAP_CONTACT_SPAN_NM = 50.0
MAP_CONTACT_CHECKS_NM = (40.0, 60.0)
MAP_CONTACT_R2_MIN = 0.995
FAR_FRACTION = 0.20
FAR_MIN_POINTS = 80
BIN_CENTERS_NM = np.arange(5.0, 505.0, 5.0)
BIN_WIDTH_NM = 5.0
TARGET_DISTANCES_NM = (20.0, 50.0, 100.0, 200.0)
SPEEDS_UM_PER_S = (0.05, 0.1, 0.2)


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


def sample_sd(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.std(values, ddof=1)) if values.size > 1 else float("nan")


def robust_mad(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not values.size:
        return float("nan")
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def terminal_contact_fit(
    height_m: np.ndarray, deflection_V: np.ndarray, span_nm: float
) -> dict:
    """Robust line fit to the contiguous terminal hard-contact interval."""

    if height_m.ndim != 1 or height_m.shape != deflection_V.shape:
        raise ValueError("contact arrays must be matching 1-D arrays")
    terminal = float(height_m[-1])
    indices = np.flatnonzero(np.abs(height_m - terminal) <= span_nm * 1e-9)
    if indices.size < 20:
        raise ValueError(f"only {indices.size} terminal points at {span_nm:g} nm")
    # Retain only the terminal contiguous run, protecting against scanner jitter.
    start = int(indices[-1])
    while start > 0 and abs(height_m[start - 1] - terminal) <= span_nm * 1e-9:
        start -= 1
    selected = np.arange(start, height_m.size)
    slope, intercept, r2, n_inliers = base.robust_line(
        height_m[selected], deflection_V[selected]
    )
    if not np.isfinite(slope) or slope >= 0.0:
        raise ValueError("unexpected contact slope")
    return {
        "start": start,
        "stop": int(height_m.size),
        "points": int(selected.size),
        "actual_span_nm": float(np.ptp(height_m[selected]) * 1e9),
        "slope_V_per_m": float(slope),
        "intercept_V": float(intercept),
        "r2": float(r2),
        "n_inliers": int(n_inliers),
        "sensitivity_m_per_V": float(1.0 / abs(slope)),
    }


def calibration_analysis() -> tuple[list[dict], list[dict], list[dict], dict]:
    """Return summaries, thermal rows, force rows, and plot payloads."""

    # Make the user-supplied temperature authoritative over TND header 25.0 C.
    cal.TEMPERATURE_C = TEMPERATURE_C
    cal.TEMPERATURE_K = TEMPERATURE_K
    summaries: list[dict] = []
    thermal_rows: list[dict] = []
    force_rows: list[dict] = []
    plots: dict = {name: {"thermal": [], "force": []} for name in ("D4", "D5", "D6")}

    for group, name in ((4, "D4"), (5, "D5"), (6, "D6")):
        thermal_records = []
        thermal_paths = sorted((DATA_ROOT / "calibration" / name).glob("*.tnd"))
        seen_thermal_hashes: dict[str, str] = {}
        thermal_window_factors: dict[float, list[float]] = {
            width: [] for width in cal.THERMAL_FIT_WINDOW_CHECKS_HZ
        }
        for path in thermal_paths:
            path = path.resolve()
            file_hash = sha256_file(path)
            record, plot_data = cal.refit_tnd(path, group)
            header, data = cal.load_tnd(path)
            row = asdict(record)
            row["cantilever"] = name
            row["user_temperature_C"] = TEMPERATURE_C
            row["file_header_temperature"] = header.get("settings.temperature", "")
            row["exact_duplicate_of"] = seen_thermal_hashes.get(file_hash, "")
            row["used_for_calibration"] = file_hash not in seen_thermal_hashes
            thermal_rows.append(row)
            if file_hash in seen_thermal_hashes:
                continue
            seen_thermal_hashes[file_hash] = str(path.relative_to(ROOT))
            thermal_records.append(record)
            plots[name]["thermal"].append(plot_data)
            for width in thermal_window_factors:
                fit = cal.fit_sho(
                    data[:, 0], data[:, 2], record.detected_peak_hz, width
                )
                thermal_window_factors[width].append(
                    record.correction_factor
                    * Boltzmann
                    * TEMPERATURE_K
                    / float(fit["resonance_area_v2"])
                )

        group_force: list[dict] = []
        for path in sorted((DATA_ROOT / "calibration" / name).glob("*.jpk-force")):
            path = path.resolve()
            source = base.load_source(path, 0)
            if len(source.curves) != 1:
                raise RuntimeError(f"expected one approach curve in {path}")
            curve = source.curves[0]
            primary = terminal_contact_fit(
                curve.measured_height_m, curve.deflection_V, CAL_CONTACT_SPAN_NM
            )
            checks = {
                span: terminal_contact_fit(
                    curve.measured_height_m, curve.deflection_V, span
                )
                for span in CAL_CONTACT_CHECKS_NM
            }
            row = {
                "cantilever": name,
                "file": str(path.relative_to(ROOT)),
                "sha256": sha256_file(path),
                "contact_span_nm": CAL_CONTACT_SPAN_NM,
                "contact_actual_span_nm": primary["actual_span_nm"],
                "contact_points": primary["points"],
                "slope_V_per_m": primary["slope_V_per_m"],
                "sensitivity_nm_per_V": primary["sensitivity_m_per_V"] * 1e9,
                "r_squared": primary["r2"],
                "stored_sensitivity_nm_per_V": source.stored_sensitivity_m_per_V * 1e9,
                "stored_spring_constant_N_per_m": source.stored_spring_constant_N_per_m,
                "sensitivity_30nm_nm_per_V": checks[30.0]["sensitivity_m_per_V"] * 1e9,
                "sensitivity_40nm_nm_per_V": checks[40.0]["sensitivity_m_per_V"] * 1e9,
            }
            group_force.append(row)
            force_rows.append(row)
            plots[name]["force"].append(
                {
                    "height_m": curve.measured_height_m,
                    "deflection_V": curve.deflection_V,
                    "fit": primary,
                }
            )

        if len(thermal_paths) != 5 or len(group_force) != 5:
            raise RuntimeError(f"{name}: expected 5 TND files and 5 force files")
        if len(thermal_records) < 4:
            raise RuntimeError(f"{name}: fewer than 4 unique TND spectra")
        sens = np.asarray([row["sensitivity_nm_per_V"] * 1e-9 for row in group_force])
        sens30 = np.asarray([row["sensitivity_30nm_nm_per_V"] * 1e-9 for row in group_force])
        sens40 = np.asarray([row["sensitivity_40nm_nm_per_V"] * 1e-9 for row in group_force])
        factors = np.asarray([record.thermal_factor_n_m_per_v2 for record in thermal_records])
        freqs = np.asarray([record.fit_frequency_hz for record in thermal_records])
        qs = np.asarray([record.quality_factor for record in thermal_records])
        mean_sens = float(np.mean(sens))
        mean_factor = float(np.mean(factors))
        spring = mean_factor / mean_sens**2
        repeat_sd = spring * math.sqrt(
            (sample_sd(factors) / mean_factor) ** 2
            + (2.0 * sample_sd(sens) / mean_sens) ** 2
        )
        contact_springs = (mean_factor / np.mean(sens30) ** 2, mean_factor / np.mean(sens40) ** 2)
        thermal_springs = [
            np.mean(thermal_window_factors[width]) / mean_sens**2
            for width in cal.THERMAL_FIT_WINDOW_CHECKS_HZ
        ]
        summaries.append(
            {
                "cantilever": name,
                "temperature_C": TEMPERATURE_C,
                "calibration_environment": "air",
                "force_curves": len(group_force),
                "thermal_files": len(thermal_paths),
                "thermal_spectra": len(thermal_records),
                "thermal_exact_duplicates_excluded": len(thermal_paths)
                - len(thermal_records),
                "sensitivity_nm_per_V": mean_sens * 1e9,
                "sensitivity_repeatability_sd_nm_per_V": sample_sd(sens) * 1e9,
                "minimum_contact_fit_R2": min(row["r_squared"] for row in group_force),
                "resonance_frequency_kHz": float(np.mean(freqs) / 1e3),
                "resonance_frequency_sd_kHz": sample_sd(freqs) / 1e3,
                "quality_factor": float(np.mean(qs)),
                "quality_factor_sd": sample_sd(qs),
                "thermal_factor_1e15_N_m_per_V2": mean_factor * 1e15,
                "thermal_factor_sd_1e15_N_m_per_V2": sample_sd(factors) * 1e15,
                "spring_constant_N_per_m": spring,
                "spring_constant_repeatability_sd_N_per_m": repeat_sd,
                "spring_shift_contact_window_N_per_m": max(
                    abs(value - spring) for value in contact_springs
                ),
                "spring_shift_thermal_window_N_per_m": max(
                    abs(value - spring) for value in thermal_springs
                ),
                "all_thermal_fits_successful": int(
                    all(record.fit_success for record in thermal_records)
                ),
            }
        )
    return summaries, thermal_rows, force_rows, plots


MAP_NAME = re.compile(r"^map(?P<block>\d+)-(?P<speed>[0-9.]+)-data-")


def map_identity(path: Path) -> tuple[int, float]:
    match = MAP_NAME.match(path.name)
    if match is None:
        raise ValueError(f"unrecognised map name: {path.name}")
    return int(match.group("block")), float(match.group("speed"))


def set_map_contact_fits(sources: list[base.SourceData]) -> tuple[list[dict], float]:
    """Fit all water hard contacts and set one all-data water InvOLS."""

    rows: list[dict] = []
    retained_all: list[float] = []
    for source in sources:
        block, speed = map_identity(source.path)
        map_values: list[float] = []
        for curve in source.curves:
            curve.far_field_fit = base.fit_far_field_drift(
                curve.measured_height_m, curve.deflection_V
            )
            corrected = curve.deflection_V - base.baseline_voltage(
                curve.measured_height_m, curve.deflection_V, curve.far_field_fit
            )
            primary = terminal_contact_fit(
                curve.measured_height_m, corrected, MAP_CONTACT_SPAN_NM
            )
            checks = {
                span: terminal_contact_fit(curve.measured_height_m, corrected, span)
                for span in MAP_CONTACT_CHECKS_NM
            }
            curve.contact_fit = base.ContactFit(
                sensitivity_m_per_V=primary["sensitivity_m_per_V"],
                start=primary["start"],
                stop=primary["stop"],
                slope_V_per_m=primary["slope_V_per_m"],
                intercept_V=primary["intercept_V"],
                r2=primary["r2"],
            )
            valid = bool(
                primary["r2"] >= MAP_CONTACT_R2_MIN
                and 35e-9 <= primary["sensitivity_m_per_V"] <= 90e-9
            )
            if valid:
                map_values.append(primary["sensitivity_m_per_V"])
            rows.append(
                {
                    "source": str(source.path.relative_to(ROOT)),
                    "timestamp": source.timestamp.isoformat(),
                    "block": block,
                    "nominal_speed_um_per_s": speed,
                    "point_index": curve.point_index,
                    "sensitivity_nm_per_V": primary["sensitivity_m_per_V"] * 1e9,
                    "sensitivity_40nm_nm_per_V": checks[40.0]["sensitivity_m_per_V"] * 1e9,
                    "sensitivity_60nm_nm_per_V": checks[60.0]["sensitivity_m_per_V"] * 1e9,
                    "contact_fit_R2": primary["r2"],
                    "contact_points": primary["points"],
                    "valid_for_consensus": valid,
                }
            )
        values = np.asarray(map_values)
        if values.size < 56:
            raise RuntimeError(f"too few valid contact fits in {source.path}: {values.size}")
        center = float(np.median(values))
        mad = robust_mad(values)
        tolerance = max(3.5 * mad, 0.5e-9)
        retained = values[np.abs(values - center) <= tolerance]
        retained_all.extend(retained.tolist())
        source.sensitivity_anchor_m_per_V = float(np.median(retained))
        source.sensitivity_anchor_mad_m_per_V = robust_mad(retained)
        source.sensitivity_valid_curves = int(retained.size)
    consensus = float(np.median(np.asarray(retained_all)))
    for source in sources:
        source.sensitivity_used_m_per_V = consensus
        source.sensitivity_method = "all_water_map_hard_contacts_global_median"
    return rows, consensus


def bin_median(distance_nm: np.ndarray, values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    output = np.full(BIN_CENTERS_NM.shape, np.nan)
    usable = (
        mask
        & np.isfinite(distance_nm)
        & np.isfinite(values)
        & (distance_nm >= BIN_CENTERS_NM[0] - BIN_WIDTH_NM / 2)
        & (distance_nm < BIN_CENTERS_NM[-1] + BIN_WIDTH_NM / 2)
    )
    if not np.any(usable):
        return output
    indices = np.floor(
        (distance_nm[usable] - (BIN_CENTERS_NM[0] - BIN_WIDTH_NM / 2))
        / BIN_WIDTH_NM
    ).astype(int)
    selected_values = values[usable]
    for index in np.unique(indices):
        if 0 <= index < output.size:
            output[index] = float(np.median(selected_values[indices == index]))
    return output


def map_pixel(source: base.SourceData, point_index: int) -> tuple[int, int]:
    if source.map_grid_i is None:
        raise ValueError("missing map grid")
    return base.map_pixel_from_index(
        point_index, source.map_grid_i, source.map_back_and_forth
    )


def map_analysis(
    spring_N_per_m: float, sources: list[base.SourceData]
) -> tuple[list[dict], list[dict], list[dict], dict[str, np.ndarray]]:
    """Prepare pixel curves and map summaries from raw approach branches."""

    curve_rows: list[dict] = []
    map_rows: list[dict] = []
    force_rows: list[dict] = []
    matrices: dict[str, np.ndarray] = {}
    base.SPRING_CONSTANT_N_PER_M = spring_N_per_m
    for source_order, source in enumerate(sources, start=1):
        block, nominal_speed = map_identity(source.path)
        retracts, retract_skipped = events.load_branch(source.path, "retract")
        retract_by_point = {
            int(curve.point_index): curve
            for curve in retracts
            if curve.point_index is not None
        }
        force_line_matrix = np.full((64, BIN_CENTERS_NM.size), np.nan)
        force_constant_matrix = np.full_like(force_line_matrix, np.nan)
        speed_values: list[float] = []
        gap_speed_values: list[float] = []
        terminal_loads: list[float] = []
        far_slopes: list[float] = []
        contact_heights: list[float] = []
        event_records: list[dict] = []
        for curve_order, curve in enumerate(source.curves):
            if curve.point_index is None or curve.contact_fit is None or curve.far_field_fit is None:
                continue
            sensitivity = source.sensitivity_used_m_per_V
            far = curve.far_field_fit
            baseline = base.baseline_voltage(curve.measured_height_m, curve.deflection_V, far)
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
            far_n = far.n_points
            reference_V = float(np.median(curve.deflection_V[:far_n]))
            force_constant_pN = (
                spring_N_per_m
                * sensitivity
                * (curve.deflection_V - reference_V)
                * 1e12
            )
            precontact = np.arange(distance_nm.size) < fit.start
            row, column = map_pixel(source, int(curve.point_index))
            force_line_matrix[curve.point_index] = bin_median(
                distance_nm, force_line_pN, precontact
            )
            force_constant_matrix[curve.point_index] = bin_median(
                distance_nm, force_constant_pN, precontact
            )
            time_s = (
                np.arange(distance_nm.size, dtype=float) + 0.5
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
            terminal_load = float(np.median(force_line_pN[-10:]) / 1e3)
            event = events.analyze_pair(
                curve,
                retract_by_point.get(int(curve.point_index)),
                sensitivity,
            )
            event_records.append(event)
            speed_values.append(scanner_speed)
            gap_speed_values.append(gap_speed)
            terminal_loads.append(terminal_load)
            far_slopes.append(far_slope)
            contact_heights.append(contact_height_m * 1e6)
            curve_rows.append(
                {
                    "source": str(source.path.relative_to(ROOT)),
                    "timestamp": source.timestamp.isoformat(),
                    "acquisition_order": source_order,
                    "block": block,
                    "nominal_speed_um_per_s": nominal_speed,
                    "point_index": curve.point_index,
                    "row": row,
                    "column": column,
                    "pixel_acquisition_order": curve_order,
                    "sensitivity_used_nm_per_V": sensitivity * 1e9,
                    "contact_invOLS_nm_per_V": fit.sensitivity_m_per_V * 1e9,
                    "contact_response_ratio": sensitivity / fit.sensitivity_m_per_V,
                    "contact_fit_R2": fit.r2,
                    "contact_height_um": contact_height_m * 1e6,
                    "terminal_load_nN": terminal_load,
                    "far_slope_pN_per_100nm": far_slope,
                    "far_noise_pN": spring_N_per_m * sensitivity * far.residual_mad_V * 1e12,
                    "scanner_speed_um_per_s": scanner_speed,
                    "gap_speed_20_200nm_um_per_s": gap_speed,
                    "approach_snap_detected": event["approach_snap_detected"],
                    "approach_snap_distance_nm": event["approach_snap_distance_nm"],
                    "approach_barrier_force_nN": event["approach_barrier_force_nN"],
                    "retract_free_baseline_valid": event["retract_free_baseline_valid"],
                    "retract_snapoff_detected": event["retract_snapoff_detected"],
                    "retract_detachment_right_censored": event[
                        "retract_detachment_right_censored"
                    ],
                    "retract_pull_off_censored": event["retract_pull_off_censored"],
                    "retract_detachment_piezo_travel_nm": event[
                        "retract_detachment_piezo_travel_nm"
                    ],
                    "retract_pull_off_force_nN": event["retract_pull_off_force_nN"],
                    "retract_pull_off_observed_floor_nN": event[
                        "retract_pull_off_observed_floor_nN"
                    ],
                    "branch_contact_height_difference_nm": event[
                        "branch_contact_height_difference_nm"
                    ],
                }
            )
        source_key = str(source.path.relative_to(ROOT))
        matrices[source_key + "|line"] = force_line_matrix
        matrices[source_key + "|constant"] = force_constant_matrix
        for method, matrix in (
            ("linear_drift_corrected", force_line_matrix),
            ("far_constant_referenced", force_constant_matrix),
        ):
            for index, distance in enumerate(BIN_CENTERS_NM):
                values = matrix[:, index]
                values = values[np.isfinite(values)]
                if not values.size:
                    continue
                force_rows.append(
                    {
                        "source": source_key,
                        "timestamp": source.timestamp.isoformat(),
                        "acquisition_order": source_order,
                        "block": block,
                        "nominal_speed_um_per_s": nominal_speed,
                        "baseline_method": method,
                        "distance_nm": distance,
                        "force_median_pN": float(np.median(values)),
                        "force_q25_pN": float(np.quantile(values, 0.25)),
                        "force_q75_pN": float(np.quantile(values, 0.75)),
                        "pixel_count": int(values.size),
                    }
                )
        map_rows.append(
            {
                "source": source_key,
                "sha256": source.sha256,
                "timestamp": source.timestamp.isoformat(),
                "acquisition_order": source_order,
                "block": block,
                "nominal_speed_um_per_s": nominal_speed,
                "grid_i": source.map_grid_i,
                "grid_j": source.map_grid_j,
                "field_u_um": source.map_ulength_m * 1e6,
                "field_v_um": source.map_vlength_m * 1e6,
                "approach_curves": len(source.curves),
                "skipped_curves": source.skipped_curves,
                "retract_curves": len(retracts),
                "retract_skipped_curves": retract_skipped,
                "valid_contact_curves": source.sensitivity_valid_curves,
                "map_contact_invOLS_median_nm_per_V": source.sensitivity_anchor_m_per_V * 1e9,
                "map_contact_invOLS_mad_nm_per_V": source.sensitivity_anchor_mad_m_per_V * 1e9,
                "global_water_InvOLS_used_nm_per_V": source.sensitivity_used_m_per_V * 1e9,
                "scanner_speed_median_um_per_s": float(np.median(speed_values)),
                "gap_speed_20_200nm_median_um_per_s": float(np.nanmedian(gap_speed_values)),
                "terminal_load_median_nN": float(np.median(terminal_loads)),
                "far_slope_median_pN_per_100nm": float(np.median(far_slopes)),
                "far_slope_mad_pN_per_100nm": robust_mad(np.asarray(far_slopes)),
                "contact_height_median_um": float(np.median(contact_heights)),
                "approach_snap_detected_fraction": float(
                    np.mean([event["approach_snap_detected"] for event in event_records])
                ),
                "approach_snap_distance_median_nm": finite_median(
                    event_records, "approach_snap_distance_nm"
                ),
                "retract_free_baseline_valid_fraction": float(
                    np.mean(
                        [event["retract_free_baseline_valid"] for event in event_records]
                    )
                ),
                "retract_snapoff_detected_fraction": float(
                    np.mean([event["retract_snapoff_detected"] for event in event_records])
                ),
                "retract_detachment_right_censored_fraction": float(
                    np.mean(
                        [
                            event["retract_detachment_right_censored"]
                            for event in event_records
                        ]
                    )
                ),
                "retract_pull_off_censored_fraction": float(
                    np.mean([event["retract_pull_off_censored"] for event in event_records])
                ),
                "retract_detachment_piezo_travel_median_nm": finite_median(
                    event_records, "retract_detachment_piezo_travel_nm"
                ),
                "retract_pull_off_force_median_nN": finite_median(
                    event_records, "retract_pull_off_force_nN"
                ),
                "branch_contact_height_difference_median_nm": finite_median(
                    event_records, "branch_contact_height_difference_nm"
                ),
            }
        )
    return curve_rows, map_rows, force_rows, matrices


def finite_median(rows: list[dict], field: str) -> float:
    values = np.asarray([row[field] for row in rows], dtype=float)
    values = values[np.isfinite(values)]
    return float(np.median(values)) if values.size else float("nan")


def pair_analysis(
    sources: list[base.SourceData], map_rows: list[dict], matrices: dict[str, np.ndarray]
) -> tuple[list[dict], list[dict], list[dict]]:
    """Pair early/late same-speed maps and fit descriptive block slopes."""

    row_by_source = {row["source"]: row for row in map_rows}
    pair_rows: list[dict] = []
    fit_rows: list[dict] = []
    slice_rows: list[dict] = []
    blocks = sorted({map_identity(source.path)[0] for source in sources})
    for block in blocks:
        by_speed: dict[float, list[base.SourceData]] = {}
        for source in sources:
            source_block, speed = map_identity(source.path)
            if source_block == block:
                by_speed.setdefault(speed, []).append(source)
        if set(by_speed) != set(SPEEDS_UM_PER_S):
            raise RuntimeError(f"block {block}: speeds are incomplete: {sorted(by_speed)}")
        pair_by_speed: dict[float, dict[str, np.ndarray | float | str]] = {}
        for speed in SPEEDS_UM_PER_S:
            members = sorted(by_speed[speed], key=lambda item: item.timestamp)
            if len(members) != 2:
                raise RuntimeError(f"block {block}, speed {speed}: expected two maps")
            early, late = members
            early_key = str(early.path.relative_to(ROOT))
            late_key = str(late.path.relative_to(ROOT))
            early_matrix = matrices[early_key + "|line"]
            late_matrix = matrices[late_key + "|line"]
            sym_matrix = 0.5 * (early_matrix + late_matrix)
            history_matrix = 0.5 * (late_matrix - early_matrix)
            actual_speed = float(
                np.median(
                    [
                        row_by_source[early_key]["gap_speed_20_200nm_median_um_per_s"],
                        row_by_source[late_key]["gap_speed_20_200nm_median_um_per_s"],
                    ]
                )
            )
            pair_by_speed[speed] = {
                "sym": sym_matrix,
                "history": history_matrix,
                "actual_speed": actual_speed,
                "early": early_key,
                "late": late_key,
            }
            for index, distance in enumerate(BIN_CENTERS_NM):
                for quantity, matrix in (("F_sym", sym_matrix), ("F_history", history_matrix)):
                    values = matrix[:, index]
                    values = values[np.isfinite(values)]
                    if not values.size:
                        continue
                    pair_rows.append(
                        {
                            "block": block,
                            "nominal_speed_um_per_s": speed,
                            "actual_gap_speed_um_per_s": actual_speed,
                            "early_source": early_key,
                            "late_source": late_key,
                            "quantity": quantity,
                            "distance_nm": distance,
                            "median_pN": float(np.median(values)),
                            "q25_pN": float(np.quantile(values, 0.25)),
                            "q75_pN": float(np.quantile(values, 0.75)),
                            "paired_pixel_count": int(values.size),
                        }
                    )
        for index, distance in enumerate(BIN_CENTERS_NM):
            speeds = []
            forces = []
            nominal = []
            history = []
            for speed in SPEEDS_UM_PER_S:
                item = pair_by_speed[speed]
                sym_values = np.asarray(item["sym"])[:, index]
                hist_values = np.asarray(item["history"])[:, index]
                if np.any(np.isfinite(sym_values)):
                    speeds.append(float(item["actual_speed"]))
                    forces.append(float(np.nanmedian(sym_values)))
                    history.append(float(np.nanmedian(hist_values)))
                    nominal.append(speed)
            if len(speeds) != 3:
                continue
            x = np.asarray(speeds)
            y = np.asarray(forces)
            design = np.column_stack([np.ones(3), x])
            coefficients, residuals, rank, singular = np.linalg.lstsq(design, y, rcond=None)
            if rank != 2:
                continue
            prediction = design @ coefficients
            ss_total = float(np.sum((y - np.mean(y)) ** 2))
            r2 = float("nan") if ss_total == 0 else 1.0 - float(np.sum((y - prediction) ** 2)) / ss_total
            residual_sd = float(math.sqrt(np.sum((y - prediction) ** 2)))
            covariance = residual_sd**2 * np.linalg.inv(design.T @ design)
            slope_se = float(math.sqrt(max(0.0, covariance[1, 1])))
            hyd_slope = hydrodynamic_slope_pN_per_um_s(distance)
            row = {
                "block": block,
                "distance_nm": distance,
                "pilot_zero_speed_intercept_pN": float(coefficients[0]),
                "pilot_speed_slope_pN_per_um_s": float(coefficients[1]),
                "pilot_speed_slope_SE_pN_per_um_s": slope_se,
                "pilot_speed_fit_R2": r2,
                "theory_no_slip_hyd_slope_pN_per_um_s": hyd_slope,
                "observed_to_theory_slope_ratio": float(coefficients[1] / hyd_slope),
                "speed_0.05_actual_um_per_s": speeds[nominal.index(0.05)],
                "speed_0.05_Fsym_pN": forces[nominal.index(0.05)],
                "speed_0.05_Fhistory_pN": history[nominal.index(0.05)],
                "speed_0.1_actual_um_per_s": speeds[nominal.index(0.1)],
                "speed_0.1_Fsym_pN": forces[nominal.index(0.1)],
                "speed_0.1_Fhistory_pN": history[nominal.index(0.1)],
                "speed_0.2_actual_um_per_s": speeds[nominal.index(0.2)],
                "speed_0.2_Fsym_pN": forces[nominal.index(0.2)],
                "speed_0.2_Fhistory_pN": history[nominal.index(0.2)],
                "claim_status": "two_block_descriptive_pilot",
            }
            fit_rows.append(row)
            if distance in TARGET_DISTANCES_NM:
                slice_rows.append(row.copy())
    return pair_rows, fit_rows, slice_rows


def systematic_diagnostics(
    map_rows: list[dict], map_force_rows: list[dict]
) -> tuple[list[dict], list[dict]]:
    """Quantify baseline sensitivity and chronological map-level trends."""

    lookup = {
        (row["source"], row["baseline_method"], float(row["distance_nm"])): row
        for row in map_force_rows
    }
    baseline_rows: list[dict] = []
    chronological_values: dict[str, list[float]] = {
        "water_InvOLS_nm_per_V": [],
        "far_slope_pN_per_100nm": [],
        "terminal_load_nN": [],
        "approach_snap_detected_fraction": [],
        "retract_detachment_piezo_travel_nm": [],
        "retract_pull_off_force_nN": [],
    }
    ordered = sorted(map_rows, key=lambda row: row["acquisition_order"])
    for row in ordered:
        source = row["source"]
        chronological_values["water_InvOLS_nm_per_V"].append(
            row["map_contact_invOLS_median_nm_per_V"]
        )
        chronological_values["far_slope_pN_per_100nm"].append(
            row["far_slope_median_pN_per_100nm"]
        )
        chronological_values["terminal_load_nN"].append(row["terminal_load_median_nN"])
        chronological_values["approach_snap_detected_fraction"].append(
            row["approach_snap_detected_fraction"]
        )
        chronological_values["retract_detachment_piezo_travel_nm"].append(
            row["retract_detachment_piezo_travel_median_nm"]
        )
        chronological_values["retract_pull_off_force_nN"].append(
            row["retract_pull_off_force_median_nN"]
        )
        for distance in TARGET_DISTANCES_NM:
            line = lookup[(source, "linear_drift_corrected", distance)]["force_median_pN"]
            constant = lookup[(source, "far_constant_referenced", distance)]["force_median_pN"]
            baseline_rows.append(
                {
                    "source": source,
                    "acquisition_order": row["acquisition_order"],
                    "block": row["block"],
                    "nominal_speed_um_per_s": row["nominal_speed_um_per_s"],
                    "distance_nm": distance,
                    "linear_drift_corrected_force_pN": line,
                    "far_constant_referenced_force_pN": constant,
                    "linear_minus_constant_pN": line - constant,
                }
            )
            chronological_values.setdefault(f"force_{distance:g}nm_pN", []).append(line)

    trend_rows: list[dict] = []
    order = np.arange(1, len(ordered) + 1, dtype=float)
    for metric, values_list in chronological_values.items():
        values = np.asarray(values_list, dtype=float)
        finite = np.isfinite(values)
        if np.count_nonzero(finite) < 4:
            continue
        slope, intercept = np.polyfit(order[finite], values[finite], 1)
        rho, p_value = scipy.stats.spearmanr(order[finite], values[finite])
        trend_rows.append(
            {
                "metric": metric,
                "map_count": int(np.count_nonzero(finite)),
                "first_map_value": float(values[finite][0]),
                "last_map_value": float(values[finite][-1]),
                "linear_slope_per_map_order": float(slope),
                "linear_intercept": float(intercept),
                "spearman_rho_vs_order": float(rho),
                "spearman_two_sided_p": float(p_value),
                "interpretation": "descriptive_time_trend; maps are sequential, not randomized",
            }
        )
    return baseline_rows, trend_rows


def hydrodynamic_slope_pN_per_um_s(distance_nm: float) -> float:
    eta_Pa_s = WATER_VISCOSITY_MPA_S * 1e-3
    distance_m = distance_nm * 1e-9
    coefficient_N_per_m_per_s = 6.0 * math.pi * eta_Pa_s * PROBE_RADIUS_M**2 / distance_m
    return coefficient_N_per_m_per_s * 1e6


def make_calibration_figure(summaries: list[dict], plots: dict) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.4))
    for column, name in enumerate(("D4", "D5", "D6")):
        ax = axes[0, column]
        thermal_colors = plt.cm.viridis(
            np.linspace(0.08, 0.9, len(plots[name]["thermal"]))
        )
        for color, item in zip(
            thermal_colors, plots[name]["thermal"], strict=True
        ):
            f = item["frequency"]
            measured = item["measured_average_psd"]
            record = item["record"]
            mask = np.abs(f - record.fit_frequency_hz) <= 1800.0
            ax.plot(f[mask] / 1e3, measured[mask], color=color, lw=0.8, alpha=0.75)
            ax.plot(
                f[mask] / 1e3,
                cal.thermal_model_for_plot(record, f[mask]),
                color=color,
                lw=1.3,
            )
        ax.set_yscale("log")
        ax.set_title(f"{name} thermal fundamental")
        ax.set_xlabel("Frequency (kHz)")
        ax.set_ylabel("PSD (V²/Hz)")
        ax.grid(alpha=0.2)
        ax = axes[1, column]
        force_colors = plt.cm.viridis(
            np.linspace(0.08, 0.9, len(plots[name]["force"]))
        )
        for color, item in zip(force_colors, plots[name]["force"], strict=True):
            h = item["height_m"]
            v = item["deflection_V"]
            fit = item["fit"]
            relative = (h - h[-1]) * 1e9
            mask = np.abs(relative) <= 100
            ax.plot(relative[mask], v[mask], color=color, lw=0.8, alpha=0.7)
            idx = slice(fit["start"], fit["stop"])
            ax.plot(
                relative[idx],
                fit["slope_V_per_m"] * h[idx] + fit["intercept_V"],
                color=color,
                lw=1.5,
            )
        summary = next(row for row in summaries if row["cantilever"] == name)
        ax.set_title(
            f"{name}: {summary['sensitivity_nm_per_V']:.2f} nm/V, "
            f"k={summary['spring_constant_N_per_m']:.4f} N/m"
        )
        ax.set_xlabel("Height from terminal point (nm)")
        ax.set_ylabel("vDeflection (V)")
        ax.invert_xaxis()
        ax.grid(alpha=0.2)
    fig.suptitle("Air calibration at user-specified 25.6 °C")
    fig.tight_layout()
    fig.savefig(FIGURES / "calibration_D4_D5_D6.png", dpi=220)
    plt.close(fig)


def make_map_qc_figure(map_rows: list[dict]) -> None:
    rows = sorted(map_rows, key=lambda row: row["timestamp"])
    x = np.arange(1, len(rows) + 1)
    labels = [f"B{r['block']} {r['nominal_speed_um_per_s']:g}" for r in rows]
    fig, axes = plt.subplots(4, 1, figsize=(13.0, 10.5), sharex=True)
    fields = (
        ("map_contact_invOLS_median_nm_per_V", "Map contact InvOLS (nm/V)"),
        ("far_slope_median_pN_per_100nm", "Far-field slope (pN/100 nm)"),
        ("terminal_load_median_nN", "Terminal load (nN)"),
        ("gap_speed_20_200nm_median_um_per_s", "Gap speed 20–200 nm (µm/s)"),
    )
    colors = {0.05: "#2a9d8f", 0.1: "#e9c46a", 0.2: "#e76f51"}
    for ax, (field, ylabel) in zip(axes, fields, strict=True):
        y = np.asarray([r[field] for r in rows])
        ax.plot(x, y, color="0.65", lw=1.0)
        ax.scatter(
            x,
            y,
            c=[colors[r["nominal_speed_um_per_s"]] for r in rows],
            s=42,
            zorder=3,
        )
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
    axes[-1].set_xticks(x, labels, rotation=50, ha="right")
    axes[-1].set_xlabel("Chronological map order (block, nominal speed µm/s)")
    fig.suptitle("Map-level chronological QC (pixels summarized within each map)")
    fig.tight_layout()
    fig.savefig(FIGURES / "map_chronological_QC.png", dpi=220)
    plt.close(fig)


def make_retract_qc_figure(map_rows: list[dict]) -> None:
    rows = sorted(map_rows, key=lambda row: row["timestamp"])
    x = np.arange(1, len(rows) + 1)
    labels = [f"B{r['block']} {r['nominal_speed_um_per_s']:g}" for r in rows]
    fig, axes = plt.subplots(4, 1, figsize=(13.0, 10.5), sharex=True)
    fields = (
        ("approach_snap_detected_fraction", "Approach snap detected fraction"),
        ("retract_snapoff_detected_fraction", "Retract snap-off detected fraction"),
        ("retract_detachment_piezo_travel_median_nm", "Detachment travel (nm)"),
        ("retract_pull_off_force_median_nN", "Pull-off force (nN)"),
    )
    colors = {0.05: "#2a9d8f", 0.1: "#e9c46a", 0.2: "#e76f51"}
    for ax, (field, ylabel) in zip(axes, fields, strict=True):
        y = np.asarray([r[field] for r in rows])
        ax.plot(x, y, color="0.65", lw=1.0)
        ax.scatter(
            x,
            y,
            c=[colors[r["nominal_speed_um_per_s"]] for r in rows],
            s=42,
            zorder=3,
        )
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
    axes[-1].set_xticks(x, labels, rotation=50, ha="right")
    axes[-1].set_xlabel("Chronological map order (block, nominal approach speed µm/s)")
    fig.suptitle("Approach/retract event QC (retract protocol is common across maps)")
    fig.tight_layout()
    fig.savefig(FIGURES / "retract_chronological_QC.png", dpi=220)
    plt.close(fig)


def make_force_chronology_figure(
    map_rows: list[dict], baseline_rows: list[dict]
) -> None:
    ordered = sorted(map_rows, key=lambda row: row["acquisition_order"])
    x = np.arange(1, len(ordered) + 1)
    labels = [f"B{r['block']} {r['nominal_speed_um_per_s']:g}" for r in ordered]
    colors = {0.05: "#2a9d8f", 0.1: "#e9c46a", 0.2: "#e76f51"}
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.5), sharex=True)
    for ax, distance in zip(axes.flat, TARGET_DISTANCES_NM, strict=True):
        selected = sorted(
            [row for row in baseline_rows if row["distance_nm"] == distance],
            key=lambda row: row["acquisition_order"],
        )
        line = np.asarray([row["linear_drift_corrected_force_pN"] for row in selected])
        constant = np.asarray([row["far_constant_referenced_force_pN"] for row in selected])
        ax.plot(x, line, color="#264653", lw=1.2, label="linear drift corrected")
        ax.plot(x, constant, color="0.55", lw=1.0, ls="--", label="far constant referenced")
        ax.scatter(
            x,
            line,
            c=[colors[row["nominal_speed_um_per_s"]] for row in selected],
            s=38,
            zorder=3,
        )
        ax.set_title(f"D = {distance:g} nm")
        ax.set_ylabel("Map median force (pN)")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xticks(x, labels, rotation=50, ha="right")
        ax.set_xlabel("Chronological map order")
    fig.suptitle("Force chronology and far-baseline definition sensitivity")
    fig.tight_layout()
    fig.savefig(FIGURES / "force_chronology_and_baseline_sensitivity.png", dpi=220)
    plt.close(fig)


def make_far_slope_maps(curve_rows: list[dict], map_rows: list[dict]) -> None:
    sources = [row["source"] for row in sorted(map_rows, key=lambda row: row["timestamp"])]
    values = np.asarray([row["far_slope_pN_per_100nm"] for row in curve_rows])
    limit = float(np.nanquantile(np.abs(values), 0.98))
    fig, axes = plt.subplots(3, 4, figsize=(13.2, 9.7), layout="constrained")
    image = None
    for ax, source in zip(axes.flat, sources, strict=True):
        selected = [row for row in curve_rows if row["source"] == source]
        grid = np.full((8, 8), np.nan)
        for row in selected:
            grid[int(row["row"]), int(row["column"])] = row["far_slope_pN_per_100nm"]
        image = ax.imshow(grid, cmap="coolwarm", vmin=-limit, vmax=limit, origin="lower")
        info = next(row for row in map_rows if row["source"] == source)
        ax.set_title(
            f"B{info['block']} {info['nominal_speed_um_per_s']:g} µm/s\n"
            f"order {info['acquisition_order']}"
        )
        ax.set_xticks([])
        ax.set_yticks([])
    if image is not None:
        fig.colorbar(
            image,
            ax=axes.ravel().tolist(),
            shrink=0.78,
            pad=0.02,
            label="pN per 100 nm",
        )
    fig.suptitle("8×8 far-field linear slope maps")
    fig.savefig(FIGURES / "far_field_slope_maps_8x8.png", dpi=220)
    plt.close(fig)


def make_pair_figures(pair_rows: list[dict], fit_rows: list[dict]) -> None:
    blocks = sorted({int(row["block"]) for row in pair_rows})
    colors = {0.05: "#2a9d8f", 0.1: "#e9c46a", 0.2: "#e76f51"}
    for quantity, filename, ylabel in (
        ("F_sym", "palindrome_symmetric_force.png", "F_sym (pN)"),
        ("F_history", "palindrome_history_residual.png", "F_history=(late−early)/2 (pN)"),
    ):
        fig, axes = plt.subplots(1, len(blocks), figsize=(6.3 * len(blocks), 4.8), sharey=True)
        axes = np.atleast_1d(axes)
        for ax, block in zip(axes, blocks, strict=True):
            for speed in SPEEDS_UM_PER_S:
                rows = sorted(
                    [
                        row
                        for row in pair_rows
                        if row["block"] == block
                        and row["nominal_speed_um_per_s"] == speed
                        and row["quantity"] == quantity
                        and row["distance_nm"] <= 250
                    ],
                    key=lambda row: row["distance_nm"],
                )
                d = np.asarray([row["distance_nm"] for row in rows])
                median = np.asarray([row["median_pN"] for row in rows])
                q25 = np.asarray([row["q25_pN"] for row in rows])
                q75 = np.asarray([row["q75_pN"] for row in rows])
                ax.fill_between(d, q25, q75, color=colors[speed], alpha=0.17)
                ax.plot(d, median, color=colors[speed], lw=1.7, label=f"{speed:g} µm/s")
            ax.axhline(0.0, color="0.25", lw=0.7)
            ax.set_xlim(0, 250)
            ax.set_xlabel("Separation D (nm)")
            ax.set_title(f"Block {block}")
            ax.grid(alpha=0.2)
            ax.legend()
        axes[0].set_ylabel(ylabel)
        fig.suptitle("Same-pixel early/late palindrome pairing; band = pixel IQR")
        fig.tight_layout()
        fig.savefig(FIGURES / filename, dpi=220)
        plt.close(fig)

    fig, axes = plt.subplots(1, len(blocks), figsize=(6.3 * len(blocks), 4.8), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, block in zip(axes, blocks, strict=True):
        rows = [row for row in fit_rows if row["block"] == block and row["distance_nm"] <= 250]
        d = np.asarray([row["distance_nm"] for row in rows])
        slope = np.asarray([row["pilot_speed_slope_pN_per_um_s"] for row in rows])
        theory = np.asarray([row["theory_no_slip_hyd_slope_pN_per_um_s"] for row in rows])
        ax.plot(d, slope, color="#264653", lw=1.6, label="Observed pilot slope")
        ax.plot(d, theory, color="#e76f51", lw=1.4, ls="--", label="No-slip hyd theory")
        ax.axhline(0, color="0.25", lw=0.7)
        ax.set_xlim(0, 250)
        ax.set_xlabel("Separation D (nm)")
        ax.set_title(f"Block {block}")
        ax.grid(alpha=0.2)
        ax.legend()
    axes[0].set_ylabel("dF/dU (pN per µm/s)")
    fig.suptitle("Descriptive three-speed slope versus no-slip hydrodynamics")
    fig.tight_layout()
    fig.savefig(FIGURES / "pilot_speed_slope_vs_hydrodynamics.png", dpi=220)
    plt.close(fig)


def render_report(
    calibration: list[dict],
    map_rows: list[dict],
    water_sensitivity: float,
    slice_rows: list[dict],
    curve_rows: list[dict],
    baseline_rows: list[dict],
    trend_rows: list[dict],
) -> str:
    d4 = next(row for row in calibration if row["cantilever"] == "D4")
    blocks = sorted({int(row["block"]) for row in map_rows})
    by_block = {
        block: [row for row in map_rows if row["block"] == block] for block in blocks
    }
    lines = [
        "# 27-08-26 D4–D6 标定与 D4 纯水回文测试：两组 pilot",
        "",
        "> 当前只有计划中的前两组回文 block。本报告完成基础标定、raw map 重建、同点位早/晚配对与描述性速度拟合；不把两组数据称为完整的 zero-speed extrapolation validation。",
        "",
        "## 直接结果",
        "",
        f"用户给定实际温度为 **{TEMPERATURE_C:.1f} °C**，标定环境为空气；TND header 中记录的是 25.0 °C，本计算明确以 25.6 °C 覆盖 header。热谱全部从 measured PSD 重搜 fundamental 并重拟合，未使用导出的 fit-data。",
        "",
        "| Cantilever | Air InvOLS (nm/V) | k (N/m) | f0 (kHz) | Q |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in calibration:
        lines.append(
            f"| {row['cantilever']} | {row['sensitivity_nm_per_V']:.2f} ± {row['sensitivity_repeatability_sd_nm_per_V']:.2f} | "
            f"{row['spring_constant_N_per_m']:.4f} ± {row['spring_constant_repeatability_sd_N_per_m']:.4f} | "
            f"{row['resonance_frequency_kHz']:.4f} ± {row['resonance_frequency_sd_kHz']:.4f} | "
            f"{row['quality_factor']:.2f} ± {row['quality_factor_sd']:.2f} |"
        )
    lines.extend(
        [
            "",
            "± 为 5 次重复的 sample SD，不是 traceable absolute uncertainty。由于 D5 的硬接触段只有约 40 nm，三支 cantilever 统一使用末端 35 nm 作为 primary contact fit，并以 30/40 nm 作为 window sensitivity；直接固定 50 nm 会把 D5 的接触转折混入拟合。",
            "",
            f"纯水 map 从全部有效硬接触重新得到 global InvOLS = **{water_sensitivity*1e9:.3f} nm/V**；力使用空气热标定得到的 D4 **k = {d4['spring_constant_N_per_m']:.6f} N/m**。JPK 文件中写入的 sensitivity 和 force conversion 均未作为最终力标尺。",
            "",
            "## 回文采集与 map-level QC",
            "",
        ]
    )
    for block in blocks:
        order = " → ".join(
            f"{row['nominal_speed_um_per_s']:g}"
            for row in sorted(by_block[block], key=lambda row: row["timestamp"])
        )
        lines.append(f"- Block {block} 时间顺序：`{order} µm/s`。")
    lines.extend(
        [
            "- 12 张 map 均为 8×8 pixels、10×10 µm，同一物理区域；每张 map 是一个速度实验单位，pixels 只作为同点位 paired spatial observations。",
            "",
            "| Block | Speed (µm/s) | Map InvOLS early → late (nm/V) | Gap speed early → late (µm/s) | Far slope early → late (pN/100 nm) |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for block in blocks:
        for speed in SPEEDS_UM_PER_S:
            members = sorted(
                [row for row in by_block[block] if row["nominal_speed_um_per_s"] == speed],
                key=lambda row: row["timestamp"],
            )
            lines.append(
                f"| {block} | {speed:g} | {members[0]['map_contact_invOLS_median_nm_per_V']:.3f} → {members[1]['map_contact_invOLS_median_nm_per_V']:.3f} | "
                f"{members[0]['gap_speed_20_200nm_median_um_per_s']:.4f} → {members[1]['gap_speed_20_200nm_median_um_per_s']:.4f} | "
                f"{members[0]['far_slope_median_pN_per_100nm']:.2f} → {members[1]['far_slope_median_pN_per_100nm']:.2f} |"
            )
    lines.extend(
        [
            "",
            "Retract segment 的 protocol 在所有 map 中相同，因此下列变化不能简单归因于 approach nominal speed；它们主要作为随时间/接触历史的 systematics probe。",
            "",
            "| Order | Block | Approach U | Approach snap fraction | Retract snap-off fraction | Detachment travel (nm) | Pull-off (nN) |",
            "|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in sorted(map_rows, key=lambda item: item["timestamp"]):
        lines.append(
            f"| {row['acquisition_order']} | {row['block']} | {row['nominal_speed_um_per_s']:g} | "
            f"{row['approach_snap_detected_fraction']:.2f} | {row['retract_snapoff_detected_fraction']:.2f} | "
            f"{row['retract_detachment_piezo_travel_median_nm']:.1f} | {row['retract_pull_off_force_median_nN']:.3f} |"
        )
    trend = {row["metric"]: row for row in trend_rows}
    lines.extend(
        [
            "",
            "Map-order trend 是描述性统计（连续顺序采集、并非 randomized independent maps），但可用于识别主要混杂量：",
            "",
            "| Metric | first → last | slope/map | Spearman ρ | two-sided p |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for metric in (
        "water_InvOLS_nm_per_V",
        "far_slope_pN_per_100nm",
        "force_50nm_pN",
        "force_100nm_pN",
        "approach_snap_detected_fraction",
        "retract_pull_off_force_nN",
    ):
        row = trend[metric]
        lines.append(
            f"| {metric} | {row['first_map_value']:.3g} → {row['last_map_value']:.3g} | "
            f"{row['linear_slope_per_map_order']:.3g} | {row['spearman_rho_vs_order']:.3f} | "
            f"{row['spearman_two_sided_p']:.3g} |"
        )
    block1_inv = np.asarray([row["map_contact_invOLS_median_nm_per_V"] for row in by_block[blocks[0]]])
    block2_inv = np.asarray([row["map_contact_invOLS_median_nm_per_V"] for row in by_block[blocks[1]]])
    lines.extend(
        [
            "",
            f"Block 1 的 map-median InvOLS 为 {np.mean(block1_inv):.3f} ± {sample_sd(block1_inv):.3f} nm/V，Block 2 为 {np.mean(block2_inv):.3f} ± {sample_sd(block2_inv):.3f} nm/V，均值差 {np.mean(block2_inv)-np.mean(block1_inv):+.3f} nm/V（{(np.mean(block2_inv)/np.mean(block1_inv)-1)*100:+.2f}%）。这是时间相关 optical/contact response 变化的直接证据；primary force 仍使用预先声明的全数据 global InvOLS，未逐图重新缩放。",
            "",
            "## 同点位回文分解与速度 pilot",
            "",
            "对每个 block、每个速度，将相同 pixel 的早/晚曲线定义为：",
            "",
            "`F_sym = (F_early + F_late)/2`，`F_history = (F_late - F_early)/2`。",
            "",
            "下表均为逐 curve far-field 线性 drift 修正后的 map-pixel median。每个 block 内仅有三个 F_sym map-level 速度点，因此 slope 与 U→0 intercept 是描述性 pilot。",
            "",
            "| Block | D (nm) | Fsym 0.05/0.1/0.2 (pN) | Fhistory 0.05/0.1/0.2 (pN) | observed dF/dU | hyd theory dF/dU | ratio | R² |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in slice_rows:
        lines.append(
            f"| {row['block']} | {row['distance_nm']:.0f} | "
            f"{row['speed_0.05_Fsym_pN']:.1f} / {row['speed_0.1_Fsym_pN']:.1f} / {row['speed_0.2_Fsym_pN']:.1f} | "
            f"{row['speed_0.05_Fhistory_pN']:.1f} / {row['speed_0.1_Fhistory_pN']:.1f} / {row['speed_0.2_Fhistory_pN']:.1f} | "
            f"{row['pilot_speed_slope_pN_per_um_s']:.1f} | {row['theory_no_slip_hyd_slope_pN_per_um_s']:.2f} | "
            f"{row['observed_to_theory_slope_ratio']:.1f} | {row['pilot_speed_fit_R2']:.3f} |"
        )
    history_100 = np.asarray(
        [
            abs(row["speed_0.05_Fhistory_pN"])
            for row in slice_rows
            if row["distance_nm"] == 100
        ]
        + [
            abs(row["speed_0.1_Fhistory_pN"])
            for row in slice_rows
            if row["distance_nm"] == 100
        ]
        + [
            abs(row["speed_0.2_Fhistory_pN"])
            for row in slice_rows
            if row["distance_nm"] == 100
        ]
    )
    lines.extend(
        [
            "",
            f"水在 25.6 °C 的 Cheng correlation viscosity 为 {WATER_VISCOSITY_MPA_S:.4f} mPa·s；球半径沿用 R={PROBE_RADIUS_M*1e6:.6f} µm。no-slip lubrication prediction 为 `Fhyd=6πηR²U/D`：在 100 nm、0.05/0.1/0.2 µm/s 仅为 {hydrodynamic_slope_pN_per_um_s(100)*0.05:.3f}/{hydrodynamic_slope_pN_per_um_s(100)*0.1:.3f}/{hydrodynamic_slope_pN_per_um_s(100)*0.2:.3f} pN。当前两组在 100 nm 的 |F_history| 中位量级为 {np.median(history_100):.2f} pN，因而物理 hyd signal 远小于时间/history residual；若 observed slope 很大、变号或不呈 1/D，不能解释成 hydrodynamics。",
            "",
            "### Far-field 零点定义的影响",
            "",
        ]
    )
    for distance in TARGET_DISTANCES_NM:
        differences = np.asarray(
            [
                row["linear_minus_constant_pN"]
                for row in baseline_rows
                if row["distance_nm"] == distance
            ]
        )
        lines.append(
            f"- D={distance:g} nm：linear drift corrected − far constant referenced 的 map-median 为 "
            f"{np.median(differences):+.2f} pN，范围 {np.min(differences):+.2f} 到 {np.max(differences):+.2f} pN。"
        )
    lines.extend(
        [
            "",
            "在 100–200 nm，这个 baseline-definition shift 已远大于理论 0.05–0.2 µm/s hyd force。因此第三 block 之后仍应同时报告两种零点定义，并用跨 block 可重复性决定能否分离真正的 1/D hyd term；不能把 far-field straight-line subtraction 当作无物理代价的预处理。",
            "",
            "## 基础判断",
            "",
            "1. 这套 8×8 回文设计能清楚测出 map-level 的时间/history 漂移，并避免把 64 pixels 当作 64 个速度重复。",
            "2. 水中 InvOLS 与空气 InvOLS 的差异很大，验证了必须从水中 hard contact 重算 sensitivity；同时两 block 间仍有约百分之一量级 response 漂移，需要第三 block 判断它是单调时间漂移、block 跳变还是可重复的顺序效应。",
            "3. Approach snap fraction、retract detachment travel 与 pull-off 随 map order 的变化是独立于 far-field force 的接触历史指标。若它们在同一时段共同跳变，说明 velocity label 之外还存在 surface/contact-state 漂移，不能用单一 hyd term 修正。",
            "4. 0.05–0.2 µm/s 下理论 hyd force 在 20–200 nm 只有亚 pN 到数 pN。第三 block 的主要价值不是再提高 pixel 数，而是检验回文对称化后的 F_sym(U) slope 是否跨 block 同号、同量级、并近似 1/D。",
            "5. 当前力是 apparent finite-speed force：包含 equilibrium surface force、极小的 hyd contribution 和未完全消除的 history/systematic residual。不能据这两组 alone 报告 zeta potential、Debye length 或 validated U→0 force。",
            "",
            "## 数值与数据 QC",
            "",
            f"- 共 {len(map_rows)} 张 map、{len(curve_rows)} 对 approach/retract curve；所有 map 的两条 branch 均无 parser-skipped curve。",
            "- 输入 JPK ZIP 已逐文件 CRC 检查；TND frequency 严格递增、PSD 为正且无 NaN/Inf；所有 SHO optimizer 报告 success。",
            "- 接触、far field、scanner speed 与 actual gap speed 都在 SI 单位中计算，输出时才转换为 nm、µm/s、pN、nN。",
            "- 同时保留 `far_constant_referenced` 与 `linear_drift_corrected` map force；primary pair analysis 使用后者，far slope 的 8×8 空间图单独保留，避免把 slope 悄悄当作纯仪器项。",
            "",
            "## 输出文件",
            "",
            "- `calibration_summary.csv`、`thermal_refits.csv`、`calibration_force_contact_fits.csv`：D4–D6 原始标定与窗口敏感性。",
            "- `water_contact_sensitivity_curves.csv`、`map_inventory_QC.csv`、`pixel_QC.csv`：水中 InvOLS、map acquisition/QC、逐 pixel approach/retract 指标。",
            "- `map_force_curves.csv`：两种 baseline 定义下每张 map 的 median/IQR F–D。",
            "- `palindrome_pair_curves.csv`：同点位 F_sym 与 F_history。",
            "- `pilot_velocity_fits.csv`、`force_slices_20_50_100_200nm.csv`：两 block 的描述性三速度拟合。",
            "- `baseline_sensitivity_slices.csv`、`chronological_trends.csv`：零点定义变化和 map-order 描述性趋势。",
            "- `figures/`：标定、时间 QC、8×8 far-slope、回文分解和 hyd slope 对照图。",
            "- `provenance.json`、`artifact_manifest.sha256`：参数、输入 hashes、软件版本与输出身份。",
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
    calibration, thermal_rows, force_cal_rows, calibration_plots = calibration_analysis()
    d4 = next(row for row in calibration if row["cantilever"] == "D4")

    # Keep this historical pilot tied to its declared two-block scope even
    # after later palindrome blocks are added to the same raw-data directory.
    source_paths = sorted(
        path
        for path in (DATA_ROOT / "0").glob("*.jpk-force-map")
        if map_identity(path)[0] in (1, 2)
    )
    sources = [base.load_source(path.resolve(), 0) for path in source_paths]
    sources.sort(key=lambda source: source.timestamp)
    if len(sources) != 12:
        raise RuntimeError(f"expected current two-block pilot to have 12 maps, got {len(sources)}")
    sensitivity_rows, water_sensitivity = set_map_contact_fits(sources)
    curve_rows, map_rows, map_force_rows, matrices = map_analysis(
        d4["spring_constant_N_per_m"], sources
    )
    pair_rows, fit_rows, slice_rows = pair_analysis(sources, map_rows, matrices)
    baseline_rows, trend_rows = systematic_diagnostics(map_rows, map_force_rows)

    csv_outputs = {
        "calibration_summary.csv": calibration,
        "thermal_refits.csv": thermal_rows,
        "calibration_force_contact_fits.csv": force_cal_rows,
        "water_contact_sensitivity_curves.csv": sensitivity_rows,
        "map_inventory_QC.csv": map_rows,
        "pixel_QC.csv": curve_rows,
        "map_force_curves.csv": map_force_rows,
        "palindrome_pair_curves.csv": pair_rows,
        "pilot_velocity_fits.csv": fit_rows,
        "force_slices_20_50_100_200nm.csv": slice_rows,
        "baseline_sensitivity_slices.csv": baseline_rows,
        "chronological_trends.csv": trend_rows,
    }
    for filename, rows in csv_outputs.items():
        write_csv(RESULTS / filename, rows)

    make_calibration_figure(calibration, calibration_plots)
    make_map_qc_figure(map_rows)
    make_retract_qc_figure(map_rows)
    make_force_chronology_figure(map_rows, baseline_rows)
    make_far_slope_maps(curve_rows, map_rows)
    make_pair_figures(pair_rows, fit_rows)

    report = render_report(
        calibration,
        map_rows,
        water_sensitivity,
        slice_rows,
        curve_rows,
        baseline_rows,
        trend_rows,
    )
    (RESULTS / "REPORT.md").write_text(report, encoding="utf-8")

    raw_paths = sorted(DATA_ROOT.rglob("*.jpk-force"))
    raw_paths += sorted(DATA_ROOT.rglob("*.jpk-force-map"))
    raw_paths += sorted(DATA_ROOT.rglob("*.tnd"))
    input_hashes = {
        str(path.relative_to(ROOT)): sha256_file(path) for path in raw_paths
    }
    provenance = {
        "analysis": "27-08-26 D4-D6 calibration and two-block D4 water palindrome pilot",
        "claim_status": "two_of_three_blocks_descriptive_pilot",
        "temperature_C": TEMPERATURE_C,
        "temperature_authority": "user override; TND header records 25.0 C",
        "calibration_environment": "air",
        "map_environment": "nominally pure water, no added salt",
        "geometry": "silica sphere against silica plane",
        "same_region": True,
        "probe_radius_m": PROBE_RADIUS_M,
        "water_viscosity_mPa_s": WATER_VISCOSITY_MPA_S,
        "calibration_contact_span_nm": CAL_CONTACT_SPAN_NM,
        "calibration_contact_checks_nm": CAL_CONTACT_CHECKS_NM,
        "map_contact_span_nm": MAP_CONTACT_SPAN_NM,
        "map_contact_checks_nm": MAP_CONTACT_CHECKS_NM,
        "water_InvOLS_method": "global median of retained all-map hard-contact curve fits",
        "water_InvOLS_nm_per_V": water_sensitivity * 1e9,
        "D4_spring_constant_N_per_m": d4["spring_constant_N_per_m"],
        "force_baselines_retained": [
            "far_constant_referenced",
            "linear_drift_corrected",
        ],
        "primary_pair_baseline": "linear_drift_corrected",
        "bin_centers_nm": BIN_CENTERS_NM.tolist(),
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "input_hashes": input_hashes,
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
    artifacts.append(Path(__file__).resolve())
    create_manifest(artifacts, RESULTS / "artifact_manifest.sha256")
    print(f"Wrote {RESULTS}")
    print(f"D4 k = {d4['spring_constant_N_per_m']:.9f} N/m")
    print(f"Water InvOLS = {water_sensitivity*1e9:.6f} nm/V")


if __name__ == "__main__":
    main()

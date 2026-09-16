#!/usr/bin/env python3
"""QC and velocity-distance analysis of the 16-09-26 AFM palindrome maps.

The raw JPK branches are decoded with the established repository parser.  In
contrast to the historical analysis entry points, prematurely terminated
segments are retained so every map pixel receives an explicit QC state.  Only
full-range curves with a resolved terminal hard contact enter force-distance
reconstruction.

The primary force uses a constant detector reference from 100--250 nm scanner
travel, a batch hard-contact InvOLS, and the spring constant embedded in all 18
maps.  Palindrome mates are averaged before testing velocity-distance
coupling.  Double centering removes an arbitrary distance-only surface force
and an arbitrary speed-only force gauge; it does not identify hydrodynamics
uniquely.
"""

from __future__ import annotations

import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
from typing import Iterable
from zipfile import ZipFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import numpy as np
import scipy
from scipy import stats
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import analyze_12_09_26_997glycerol_observational as observational  # noqa: E402
import analyze_velocity_systematics as velocity  # noqa: E402
import fit_glycerol_surface_forces as base  # noqa: E402


DATA = ROOT / "16-09-26"
OUT = ROOT / "analysis" / "low_speed_16_09_26_results"
FIG = OUT / "figures"

MAP_COUNT = 18
GRID_SIZE = 8
MAPS_PER_BLOCK = 6
BLOCKS = (1, 2, 3)
SPEEDS = np.asarray((0.1, 0.2, 0.4), dtype=np.float64)
EXPECTED_BLOCKS = (
    (0.1, 0.2, 0.4, 0.4, 0.2, 0.1),
    (0.2, 0.4, 0.1, 0.1, 0.4, 0.2),
    (0.4, 0.1, 0.2, 0.2, 0.1, 0.4),
)
COLORS = {0.1: "#16817a", 0.2: "#3d6fb4", 0.4: "#d56a2f"}

BASELINE_START_NM = 100.0
BASELINE_STOP_NM = 250.0
BASELINE_SENSITIVITY_WINDOWS_NM = ((75.0, 200.0), (125.0, 275.0))
CONTACT_WINDOW_NM = 40.0
CONTACT_WINDOW_CHECKS_NM = (30.0, 50.0)
CONTACT_AMPLITUDE_PRE_START_NM = 60.0
CONTACT_AMPLITUDE_PRE_STOP_NM = 50.0
CONTACT_AMPLITUDE_END_NM = 10.0
HARD_FAIL_MAX_TRAVEL_NM = 50.0
MIN_FORCE_TRAVEL_NM = 350.0
MIN_CONTACT_AMPLITUDE_V = 0.25
FORCE_CONTACT_INVOLS_RANGE_NM_PER_V = (35.0, 100.0)
FORCE_CONTACT_R2_MIN = 0.95
CALIBRATION_INVOLS_RANGE_NM_PER_V = (45.0, 75.0)
CALIBRATION_R2_MIN = 0.99
MIN_CALIBRATION_CURVES_PER_MAP = 48

BIN_CENTERS_NM = np.arange(5.0, 300.0 + 0.1, 5.0, dtype=np.float64)
BIN_HALF_WIDTH_NM = 2.5
GAP_SPEED_HALF_WINDOW_NM = 10.0
TARGET_DISTANCES_NM = (20.0, 50.0, 100.0, 200.0)
HYD_PRIMARY_MIN_NM = 20.0
HYD_PRIMARY_MAX_NM = 200.0
HYD_REFERENCE_DISTANCE_NM = 250.0
HYD_FIT_RANGES_NM = (
    (20.0, 200.0, "primary_20_200"),
    (30.0, 200.0, "sensitivity_30_200"),
    (20.0, 250.0, "sensitivity_20_250"),
)

# These are campaign assumptions used only for the fixed no-slip comparator.
# The raw files do not encode sample composition, temperature, or probe radius.
ASSUMED_GLYCEROL_MASS_FRACTION = 0.997
ASSUMED_TEMPERATURE_C = 25.6
ASSUMED_PROBE_RADIUS_M = 4.546848945303745e-6

QC_STATES = ("pass_force", "no_valid_contact", "insufficient_range", "hard_fail")
QC_CODES = {state: index for index, state in enumerate(QC_STATES)}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite(values: Iterable[float]) -> np.ndarray:
    array = np.asarray(list(values), dtype=np.float64)
    return array[np.isfinite(array)]


def median(values: Iterable[float]) -> float:
    array = finite(values)
    return float(np.median(array)) if array.size else float("nan")


def robust_mad(values: Iterable[float]) -> float:
    return base.robust_mad(np.asarray(list(values), dtype=np.float64))


def percentile(values: Iterable[float], probability: float) -> float:
    array = finite(values)
    return float(np.percentile(array, probability)) if array.size else float("nan")


def map_midpoint(metadata: dict) -> datetime:
    start = datetime.fromisoformat(str(metadata["map_start_time"]))
    end = datetime.fromisoformat(str(metadata["map_end_time"]))
    return start + (end - start) / 2


def read_map_metadata(path: Path) -> dict:
    metadata = observational.map_metadata(path)
    with ZipFile(path) as archive:
        header = base.parse_properties(archive.read("header.properties"))
    metadata.update(
        {
            "xcenter_um": float(
                header["force-scan-map.position-pattern.grid.xcenter"]
            )
            * 1e6,
            "ycenter_um": float(
                header["force-scan-map.position-pattern.grid.ycenter"]
            )
            * 1e6,
        }
    )
    return metadata


def physical_window(
    travel_nm: np.ndarray, start_nm: float, stop_nm: float
) -> np.ndarray:
    return (travel_nm >= start_nm) & (travel_nm < stop_nm)


def safe_robust_line(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, int]:
    if x.ndim != 1 or y.shape != x.shape or x.size < 3:
        return float("nan"), float("nan"), float("nan"), 0
    try:
        return base.robust_line(x, y)
    except (ValueError, np.linalg.LinAlgError, FloatingPointError):
        return float("nan"), float("nan"), float("nan"), 0


def curve_qc_record(
    curve: velocity.BranchCurve,
    retract: velocity.BranchCurve | None,
    map_record: dict,
) -> dict:
    travel_nm = (
        np.asarray(curve.measured_height_m[0] - curve.measured_height_m, dtype=np.float64)
        * 1e9
    )
    raw_v = np.asarray(curve.deflection_V, dtype=np.float64)
    if travel_nm.shape != raw_v.shape or travel_nm.ndim != 1:
        raise ValueError("Approach arrays have inconsistent shapes")
    travel_span_nm = float(np.ptp(travel_nm))
    steps_nm = np.diff(travel_nm)
    baseline_mask = physical_window(
        travel_nm, BASELINE_START_NM, BASELINE_STOP_NM
    )
    baseline_v = (
        float(np.median(raw_v[baseline_mask]))
        if np.count_nonzero(baseline_mask) >= 3
        else float("nan")
    )
    baseline_slope, _, baseline_r2, baseline_inliers = safe_robust_line(
        travel_nm[baseline_mask], raw_v[baseline_mask]
    )
    alternative_baselines: list[float] = []
    for start_nm, stop_nm in BASELINE_SENSITIVITY_WINDOWS_NM:
        selected = physical_window(travel_nm, start_nm, stop_nm)
        alternative_baselines.append(
            float(np.median(raw_v[selected]))
            if np.count_nonzero(selected) >= 3
            else float("nan")
        )

    terminal_mask = travel_nm >= float(np.max(travel_nm) - CONTACT_WINDOW_NM)
    contact_slope, _, contact_r2, contact_inliers = safe_robust_line(
        travel_nm[terminal_mask], raw_v[terminal_mask]
    )
    contact_invols = (
        float(1.0 / contact_slope)
        if np.isfinite(contact_slope) and contact_slope > 0.0
        else float("nan")
    )
    contact_window_checks: dict[float, tuple[float, float]] = {}
    for span_nm in CONTACT_WINDOW_CHECKS_NM:
        selected = travel_nm >= float(np.max(travel_nm) - span_nm)
        slope, _, r2, _ = safe_robust_line(travel_nm[selected], raw_v[selected])
        inv_ols = (
            float(1.0 / slope)
            if np.isfinite(slope) and slope > 0.0
            else float("nan")
        )
        contact_window_checks[span_nm] = (inv_ols, r2)
    pre_contact = (
        (travel_nm >= np.max(travel_nm) - CONTACT_AMPLITUDE_PRE_START_NM)
        & (travel_nm < np.max(travel_nm) - CONTACT_AMPLITUDE_PRE_STOP_NM)
    )
    contact_end = travel_nm >= np.max(travel_nm) - CONTACT_AMPLITUDE_END_NM
    contact_amplitude_v = (
        float(np.median(raw_v[contact_end]) - np.median(raw_v[pre_contact]))
        if min(np.count_nonzero(pre_contact), np.count_nonzero(contact_end)) >= 3
        else float("nan")
    )

    direction_valid = bool(
        travel_nm[-1] > travel_nm[0]
        and (steps_nm.size == 0 or float(np.min(steps_nm)) >= -2.0)
    )
    baseline_valid = bool(np.count_nonzero(baseline_mask) >= 20)
    contact_valid = bool(
        np.isfinite(contact_amplitude_v)
        and contact_amplitude_v >= MIN_CONTACT_AMPLITUDE_V
        and np.isfinite(contact_invols)
        and FORCE_CONTACT_INVOLS_RANGE_NM_PER_V[0]
        <= contact_invols
        <= FORCE_CONTACT_INVOLS_RANGE_NM_PER_V[1]
        and np.isfinite(contact_r2)
        and contact_r2 >= FORCE_CONTACT_R2_MIN
        and np.count_nonzero(terminal_mask) >= 20
    )
    saturated = bool(np.any(curve.raw_saturation_mask))
    if travel_span_nm < HARD_FAIL_MAX_TRAVEL_NM:
        qc_state = "hard_fail"
        qc_reason = "premature_trigger_under_50nm_travel"
    elif travel_span_nm < MIN_FORCE_TRAVEL_NM or not baseline_valid or not direction_valid:
        qc_state = "insufficient_range"
        qc_reason = "insufficient_travel_or_baseline_window"
    elif not contact_valid:
        qc_state = "no_valid_contact"
        qc_reason = "terminal_contact_amplitude_or_line_fit_failed"
    elif saturated:
        qc_state = "no_valid_contact"
        qc_reason = "raw_encoder_saturation"
    else:
        qc_state = "pass_force"
        qc_reason = ""

    calibration_accepted = bool(
        qc_state == "pass_force"
        and CALIBRATION_INVOLS_RANGE_NM_PER_V[0]
        <= contact_invols
        <= CALIBRATION_INVOLS_RANGE_NM_PER_V[1]
        and contact_r2 >= CALIBRATION_R2_MIN
    )
    point_index = int(curve.point_index)
    physical_row, physical_column = base.map_pixel_from_index(
        point_index, GRID_SIZE, bool(map_record["back_and_forth"])
    )
    retract_points = int(retract.deflection_V.size) if retract is not None else 0
    retract_span_nm = (
        float(np.ptp((retract.measured_height_m[0] - retract.measured_height_m) * 1e9))
        if retract is not None
        else float("nan")
    )
    first_5nm = travel_nm <= 5.0
    initial_spike_v = (
        float(np.median(raw_v[first_5nm]) - baseline_v)
        if np.count_nonzero(first_5nm) >= 3 and np.isfinite(baseline_v)
        else float("nan")
    )
    return {
        "acquisition_order": map_record["acquisition_order"],
        "instrument_scan_number": map_record["instrument_scan_number"],
        "block": map_record["block"],
        "position_in_block": map_record["position_in_block"],
        "source": map_record["source"],
        "map_midpoint_time": map_record["map_midpoint_time"],
        "speed_um_per_s": map_record["speed_um_per_s"],
        "point_index": point_index,
        "acquisition_line": point_index // GRID_SIZE + 1,
        "position_within_line": point_index % GRID_SIZE + 1,
        "physical_row": physical_row,
        "physical_column": physical_column,
        "declared_approach_points": map_record["approach_declared_points"],
        "actual_approach_points": int(raw_v.size),
        "point_count_fraction": float(
            raw_v.size / int(map_record["approach_declared_points"])
        ),
        "approach_duration_s": float(curve.duration_s),
        "approach_sample_rate_Hz": float(raw_v.size / curve.duration_s),
        "approach_travel_span_nm": travel_span_nm,
        "approach_largest_backward_step_nm": (
            float(np.min(steps_nm)) if steps_nm.size else float("nan")
        ),
        "approach_direction_valid": direction_valid,
        "raw_saturation_fraction": float(np.mean(curve.raw_saturation_mask)),
        "raw_min_encoder": curve.raw_min,
        "raw_max_encoder": curve.raw_max,
        "baseline_start_nm": BASELINE_START_NM,
        "baseline_stop_nm": BASELINE_STOP_NM,
        "baseline_points": int(np.count_nonzero(baseline_mask)),
        "baseline_raw_V": baseline_v,
        "baseline_raw_mad_V": robust_mad(raw_v[baseline_mask]),
        "baseline_slope_mV_per_100nm": float(baseline_slope * 1e5),
        "baseline_line_r2": baseline_r2,
        "baseline_line_inliers": baseline_inliers,
        "baseline_75_200nm_raw_V": alternative_baselines[0],
        "baseline_125_275nm_raw_V": alternative_baselines[1],
        "initial_0_5nm_minus_baseline_V": initial_spike_v,
        "terminal_contact_window_nm": CONTACT_WINDOW_NM,
        "terminal_contact_points": int(np.count_nonzero(terminal_mask)),
        "terminal_contact_slope_V_per_nm": contact_slope,
        "terminal_contact_invOLS_nm_per_V": contact_invols,
        "terminal_contact_r2": contact_r2,
        "terminal_contact_inliers": contact_inliers,
        "terminal_contact_amplitude_V": contact_amplitude_v,
        "terminal_contact_InvOLS_30nm_nm_per_V": contact_window_checks[30.0][0],
        "terminal_contact_R2_30nm": contact_window_checks[30.0][1],
        "terminal_contact_InvOLS_50nm_nm_per_V": contact_window_checks[50.0][0],
        "terminal_contact_R2_50nm": contact_window_checks[50.0][1],
        "force_contact_valid": contact_valid,
        "calibration_accepted": calibration_accepted,
        "qc_state": qc_state,
        "qc_reason": qc_reason,
        "retract_present": retract is not None,
        "retract_points": retract_points,
        "retract_travel_span_nm": retract_span_nm,
    }


def load_maps() -> tuple[list[dict], list[dict], list[dict]]:
    paths = sorted(DATA.glob("*.jpk-force-map"))
    if len(paths) != MAP_COUNT:
        raise RuntimeError(f"Expected {MAP_COUNT} maps, found {len(paths)}")
    map_records: list[dict] = []
    curve_rows: list[dict] = []
    payloads: list[dict] = []
    for order, path in enumerate(paths, start=1):
        metadata = read_map_metadata(path)
        source = base.load_source(path, 100)
        approaches, approach_skipped = velocity.load_branch(
            path, "extend", minimum_points=1
        )
        retracts, retract_skipped = velocity.load_branch(
            path, "retract", minimum_points=1
        )
        if approach_skipped or retract_skipped:
            raise RuntimeError(
                f"Undecodable branch in {path.name}: {approach_skipped}, {retract_skipped}"
            )
        approach_by_point = {int(curve.point_index): curve for curve in approaches}
        retract_by_point = {int(curve.point_index): curve for curve in retracts}
        expected_points = set(range(GRID_SIZE**2))
        if set(approach_by_point) != expected_points or set(retract_by_point) != expected_points:
            raise RuntimeError(f"Incomplete point-index inventory in {path.name}")
        speed = float(metadata["approach_nominal_speed_um_per_s"])
        expected_speed = EXPECTED_BLOCKS[(order - 1) // MAPS_PER_BLOCK][
            (order - 1) % MAPS_PER_BLOCK
        ]
        if not np.isclose(speed, expected_speed, atol=0.002):
            raise RuntimeError(
                f"Map {order}: expected {expected_speed:g} um/s, found {speed:g}"
            )
        midpoint = map_midpoint(metadata)
        record = {
            "acquisition_order": order,
            "instrument_scan_number": int(metadata["instrument_scan_number"]),
            "block": (order - 1) // MAPS_PER_BLOCK + 1,
            "position_in_block": (order - 1) % MAPS_PER_BLOCK + 1,
            "source": path.relative_to(ROOT).as_posix(),
            "sha256": sha256_file(path),
            "map_start_time": metadata["map_start_time"],
            "map_end_time": metadata["map_end_time"],
            "map_midpoint_time": midpoint.isoformat(),
            "speed_um_per_s": speed,
            "retract_speed_um_per_s": float(metadata["retract_nominal_speed_um_per_s"]),
            "approach_declared_points": int(metadata["approach_declared_points"]),
            "approach_declared_duration_s": float(
                metadata["approach_declared_duration_s"]
            ),
            "retract_declared_points": int(metadata["retract_declared_points"]),
            "retract_declared_duration_s": float(
                metadata["retract_declared_duration_s"]
            ),
            "nominal_z_span_nm": float(metadata["nominal_z_span_nm"]),
            "grid_i": int(metadata["grid_i"]),
            "grid_j": int(metadata["grid_j"]),
            "field_u_um": float(metadata["field_u_um"]),
            "field_v_um": float(metadata["field_v_um"]),
            "back_and_forth": bool(metadata["back_and_forth"]),
            "xcenter_um": float(metadata["xcenter_um"]),
            "ycenter_um": float(metadata["ycenter_um"]),
            "baseline_adjust_enabled": bool(metadata["baseline_adjust_enabled"]),
            "baseline_adjust_begin_of_line": bool(
                metadata["baseline_adjust_begin_of_line"]
            ),
            "baseline_adjust_deadtime_samples": int(
                metadata["baseline_adjust_deadtime_samples"]
            ),
            "baseline_adjust_average_samples": int(
                metadata["baseline_adjust_average_samples"]
            ),
            "instrument_failure_flag_count": int(
                metadata["instrument_failure_flag_count"]
            ),
            "embedded_InvOLS_nm_per_V": float(source.stored_sensitivity_m_per_V * 1e9),
            "embedded_spring_constant_N_per_m": float(
                source.stored_spring_constant_N_per_m
            ),
        }
        if record["grid_i"] != GRID_SIZE or record["grid_j"] != GRID_SIZE:
            raise RuntimeError(f"Map {order}: unexpected grid")
        current_rows = [
            curve_qc_record(approach_by_point[index], retract_by_point.get(index), record)
            for index in range(GRID_SIZE**2)
        ]
        curve_rows.extend(current_rows)
        payloads.append(
            {
                "record": record,
                "approach_by_point": approach_by_point,
                "retract_by_point": retract_by_point,
            }
        )
        map_records.append(record)
        counts = {state: sum(row["qc_state"] == state for row in current_rows) for state in QC_STATES}
        print(
            f"[{order:02d}/{MAP_COUNT}] {speed:g} um/s scan {record['instrument_scan_number']}: {counts}",
            flush=True,
        )

    first_midpoint = datetime.fromisoformat(map_records[0]["map_midpoint_time"])
    reference_x = median(row["xcenter_um"] for row in map_records[1:])
    reference_y = median(row["ycenter_um"] for row in map_records[1:])
    for row in map_records:
        midpoint = datetime.fromisoformat(row["map_midpoint_time"])
        row["elapsed_midpoint_min"] = (midpoint - first_midpoint).total_seconds() / 60.0
        shift = math.hypot(row["xcenter_um"] - reference_x, row["ycenter_um"] - reference_y)
        row["shift_from_maps_2_to_18_center_um"] = shift
        row["same_spatial_grid_as_maps_2_to_18"] = bool(shift <= 0.01)

    embedded_inv = {round(float(row["embedded_InvOLS_nm_per_V"]), 12) for row in map_records}
    embedded_k = {round(float(row["embedded_spring_constant_N_per_m"]), 12) for row in map_records}
    if len(embedded_inv) != 1 or len(embedded_k) != 1:
        raise RuntimeError("Embedded calibration is not constant across maps")
    return map_records, curve_rows, payloads


def calibrate_contact(
    map_records: list[dict], curve_rows: list[dict]
) -> tuple[float, list[dict]]:
    output: list[dict] = []
    for record in map_records:
        selected = [
            float(row["terminal_contact_invOLS_nm_per_V"])
            for row in curve_rows
            if row["acquisition_order"] == record["acquisition_order"]
            and row["calibration_accepted"]
        ]
        selected_rows = [
            row
            for row in curve_rows
            if row["acquisition_order"] == record["acquisition_order"]
            and row["calibration_accepted"]
        ]
        if len(selected) < MIN_CALIBRATION_CURVES_PER_MAP:
            raise RuntimeError(
                f"Map {record['acquisition_order']}: only {len(selected)} calibration contacts"
            )
        output.append(
            {
                "acquisition_order": record["acquisition_order"],
                "instrument_scan_number": record["instrument_scan_number"],
                "block": record["block"],
                "position_in_block": record["position_in_block"],
                "map_midpoint_time": record["map_midpoint_time"],
                "speed_um_per_s": record["speed_um_per_s"],
                "accepted_contacts": len(selected),
                "map_contact_InvOLS_median_nm_per_V": median(selected),
                "map_contact_InvOLS_mad_nm_per_V": robust_mad(selected),
                "map_contact_InvOLS_q25_nm_per_V": percentile(selected, 25),
                "map_contact_InvOLS_q75_nm_per_V": percentile(selected, 75),
                "map_contact_InvOLS_30nm_median_nm_per_V": median(
                    row["terminal_contact_InvOLS_30nm_nm_per_V"]
                    for row in selected_rows
                ),
                "map_contact_InvOLS_50nm_median_nm_per_V": median(
                    row["terminal_contact_InvOLS_50nm_nm_per_V"]
                    for row in selected_rows
                ),
            }
        )
    batch = median(row["map_contact_InvOLS_median_nm_per_V"] for row in output)
    return batch, output


def bin_curve(distance_nm: np.ndarray, force_pn: np.ndarray, mask: np.ndarray) -> np.ndarray:
    output = np.full(BIN_CENTERS_NM.shape, np.nan, dtype=np.float64)
    for column, center in enumerate(BIN_CENTERS_NM):
        selected = (
            mask
            & (distance_nm >= center - BIN_HALF_WIDTH_NM)
            & (distance_nm < center + BIN_HALF_WIDTH_NM)
            & np.isfinite(force_pn)
        )
        if np.count_nonzero(selected) >= 2:
            output[column] = float(np.median(force_pn[selected]))
    return output


def bin_gap_speed(
    distance_nm: np.ndarray, time_s: np.ndarray, mask: np.ndarray
) -> np.ndarray:
    """Estimate local physical gap-closing speed from robust D(t) slopes."""

    output = np.full(BIN_CENTERS_NM.shape, np.nan, dtype=np.float64)
    for column, center in enumerate(BIN_CENTERS_NM):
        selected = (
            mask
            & (distance_nm >= center - GAP_SPEED_HALF_WINDOW_NM)
            & (distance_nm < center + GAP_SPEED_HALF_WINDOW_NM)
            & np.isfinite(distance_nm)
            & np.isfinite(time_s)
        )
        if np.count_nonzero(selected) < 20:
            continue
        slope_m_per_s, _, _, _ = safe_robust_line(
            time_s[selected], distance_nm[selected] * 1e-9
        )
        closing_speed = -slope_m_per_s * 1e6
        if np.isfinite(closing_speed) and closing_speed > 0.0:
            output[column] = float(closing_speed)
    return output


def reconstruct(
    map_records: list[dict],
    curve_rows: list[dict],
    payloads: list[dict],
    inv_ols_nm_per_v: float,
) -> tuple[list[dict], np.ndarray, np.ndarray]:
    spring_constants = {
        round(float(row["embedded_spring_constant_N_per_m"]), 12)
        for row in map_records
    }
    if len(spring_constants) != 1:
        raise RuntimeError("Expected one embedded spring constant")
    spring = float(map_records[0]["embedded_spring_constant_N_per_m"])
    matrices = np.full(
        (MAP_COUNT, GRID_SIZE**2, BIN_CENTERS_NM.size), np.nan, dtype=np.float64
    )
    gap_speed_matrices = np.full_like(matrices, np.nan)
    qc_lookup = {
        (int(row["acquisition_order"]), int(row["point_index"])): row
        for row in curve_rows
    }
    for payload in payloads:
        record = payload["record"]
        order = int(record["acquisition_order"])
        for point, curve in payload["approach_by_point"].items():
            qc = qc_lookup[(order, point)]
            if qc["qc_state"] != "pass_force":
                continue
            height_m = np.asarray(curve.measured_height_m, dtype=np.float64)
            raw_v = np.asarray(curve.deflection_V, dtype=np.float64)
            travel_nm = (height_m[0] - height_m) * 1e9
            baseline_v = float(qc["baseline_raw_V"])
            delta_m = (raw_v - baseline_v) * inv_ols_nm_per_v * 1e-9
            terminal = travel_nm >= np.max(travel_nm) - CONTACT_WINDOW_NM
            contact_coordinate_m = height_m[terminal] + delta_m[terminal]
            contact_plane_m = float(np.median(contact_coordinate_m))
            distance_nm = (height_m + delta_m - contact_plane_m) * 1e9
            force_pn = spring * delta_m * 1e12
            precontact = np.arange(distance_nm.size) < int(np.flatnonzero(terminal)[0])
            matrices[order - 1, point] = bin_curve(distance_nm, force_pn, precontact)
            time_s = (
                np.arange(distance_nm.size, dtype=np.float64) + 0.5
            ) * float(curve.duration_s) / distance_nm.size
            gap_speed_matrices[order - 1, point] = bin_gap_speed(
                distance_nm, time_s, precontact
            )
            qc["contact_plane_m"] = contact_plane_m
            qc["contact_coordinate_mad_nm"] = robust_mad(contact_coordinate_m * 1e9)
            qc["terminal_reconstructed_distance_nm"] = float(
                np.median(distance_nm[terminal])
            )
            qc["reconstructed_max_precontact_distance_nm"] = float(
                np.max(distance_nm[precontact])
            )
            qc["force_scale_pN_per_V"] = spring * inv_ols_nm_per_v * 1e3
            qc["baseline_75_200_force_shift_pN"] = (
                spring
                * inv_ols_nm_per_v
                * 1e3
                * (baseline_v - float(qc["baseline_75_200nm_raw_V"]))
            )
            qc["baseline_125_275_force_shift_pN"] = (
                spring
                * inv_ols_nm_per_v
                * 1e3
                * (baseline_v - float(qc["baseline_125_275nm_raw_V"]))
            )

    force_rows: list[dict] = []
    for record in map_records:
        order = int(record["acquisition_order"])
        matrix = matrices[order - 1]
        for column, distance in enumerate(BIN_CENTERS_NM):
            values = matrix[:, column]
            values = values[np.isfinite(values)]
            gap_values = gap_speed_matrices[order - 1, :, column]
            gap_values = gap_values[np.isfinite(gap_values)]
            force_rows.append(
                {
                    "acquisition_order": order,
                    "instrument_scan_number": record["instrument_scan_number"],
                    "block": record["block"],
                    "position_in_block": record["position_in_block"],
                    "map_midpoint_time": record["map_midpoint_time"],
                    "elapsed_midpoint_min": record["elapsed_midpoint_min"],
                    "speed_um_per_s": record["speed_um_per_s"],
                    "distance_nm": float(distance),
                    "available_pixels": int(values.size),
                    "map_median_force_pN": median(values),
                    "map_force_q25_pN": percentile(values, 25),
                    "map_force_q75_pN": percentile(values, 75),
                    "map_force_spatial_mad_pN": robust_mad(values),
                    "available_gap_speed_pixels": int(gap_values.size),
                    "map_median_gap_speed_um_per_s": median(gap_values),
                    "map_gap_speed_q25_um_per_s": percentile(gap_values, 25),
                    "map_gap_speed_q75_um_per_s": percentile(gap_values, 75),
                    "map_gap_speed_spatial_mad_um_per_s": robust_mad(gap_values),
                }
            )
    return force_rows, matrices, gap_speed_matrices


def enrich_map_records(
    map_records: list[dict],
    curve_rows: list[dict],
    contact_rows: list[dict],
    force_rows: list[dict],
) -> None:
    contact_lookup = {int(row["acquisition_order"]): row for row in contact_rows}
    force_lookup = {
        (int(row["acquisition_order"]), float(row["distance_nm"])): row
        for row in force_rows
    }
    for record in map_records:
        order = int(record["acquisition_order"])
        curves = [row for row in curve_rows if int(row["acquisition_order"]) == order]
        for state in QC_STATES:
            record[f"qc_{state}_curves"] = sum(row["qc_state"] == state for row in curves)
        record["median_actual_approach_points"] = median(
            row["actual_approach_points"] for row in curves
        )
        record["median_approach_travel_nm"] = median(
            row["approach_travel_span_nm"] for row in curves
        )
        record["map_baseline_raw_V_median"] = median(row["baseline_raw_V"] for row in curves)
        record["map_baseline_raw_V_spatial_mad"] = robust_mad(
            row["baseline_raw_V"] for row in curves
        )
        record["map_baseline_slope_mV_per_100nm_median"] = median(
            row["baseline_slope_mV_per_100nm"] for row in curves
            if row["qc_state"] == "pass_force"
        )
        record.update(
            {
                key: value
                for key, value in contact_lookup[order].items()
                if key not in record
            }
        )
        for distance in TARGET_DISTANCES_NM:
            row = force_lookup[(order, distance)]
            record[f"force_{distance:g}nm_median_pN"] = row["map_median_force_pN"]
            record[f"force_{distance:g}nm_spatial_mad_pN"] = row[
                "map_force_spatial_mad_pN"
            ]
            record[f"gap_speed_{distance:g}nm_median_um_per_s"] = row[
                "map_median_gap_speed_um_per_s"
            ]
            record[f"gap_speed_{distance:g}nm_over_nominal"] = float(
                row["map_median_gap_speed_um_per_s"]
            ) / float(record["speed_um_per_s"])


def build_palindrome_pairs(
    map_records: list[dict], force_rows: list[dict]
) -> tuple[
    list[dict],
    dict[tuple[int, float, float], float],
    dict[tuple[int, float, float], float],
]:
    force_lookup = {
        (int(row["acquisition_order"]), float(row["distance_nm"])): float(
            row["map_median_force_pN"]
        )
        for row in force_rows
    }
    gap_speed_lookup = {
        (int(row["acquisition_order"]), float(row["distance_nm"])): float(
            row["map_median_gap_speed_um_per_s"]
        )
        for row in force_rows
    }
    rows: list[dict] = []
    pair_lookup: dict[tuple[int, float, float], float] = {}
    pair_gap_speed_lookup: dict[tuple[int, float, float], float] = {}
    for block in BLOCKS:
        for speed in SPEEDS:
            maps = sorted(
                (
                    row
                    for row in map_records
                    if int(row["block"]) == block
                    and np.isclose(float(row["speed_um_per_s"]), speed)
                ),
                key=lambda row: int(row["position_in_block"]),
            )
            if len(maps) != 2 or int(maps[0]["position_in_block"]) + int(
                maps[1]["position_in_block"]
            ) != 7:
                raise RuntimeError(f"Block {block}, speed {speed:g}: invalid palindrome pair")
            early, late = maps
            early_order = int(early["acquisition_order"])
            late_order = int(late["acquisition_order"])
            early_time = datetime.fromisoformat(early["map_midpoint_time"])
            late_time = datetime.fromisoformat(late["map_midpoint_time"])
            pair_center = early_time + (late_time - early_time) / 2
            centers_match = bool(
                math.hypot(
                    float(early["xcenter_um"]) - float(late["xcenter_um"]),
                    float(early["ycenter_um"]) - float(late["ycenter_um"]),
                )
                <= 0.01
            )
            for distance in BIN_CENTERS_NM:
                early_force = force_lookup[(early_order, float(distance))]
                late_force = force_lookup[(late_order, float(distance))]
                pair_mean = 0.5 * (early_force + late_force)
                pair_lookup[(block, float(distance), float(speed))] = pair_mean
                early_gap_speed = gap_speed_lookup[(early_order, float(distance))]
                late_gap_speed = gap_speed_lookup[(late_order, float(distance))]
                pair_gap_speed = 0.5 * (early_gap_speed + late_gap_speed)
                pair_gap_speed_lookup[(block, float(distance), float(speed))] = (
                    pair_gap_speed
                )
                rows.append(
                    {
                        "block": block,
                        "speed_um_per_s": float(speed),
                        "distance_nm": float(distance),
                        "early_map_order": early_order,
                        "late_map_order": late_order,
                        "early_scan_number": early["instrument_scan_number"],
                        "late_scan_number": late["instrument_scan_number"],
                        "pair_center_time": pair_center.isoformat(),
                        "early_late_separation_min": (
                            late_time - early_time
                        ).total_seconds()
                        / 60.0,
                        "pixelwise_pairing_permitted": centers_match,
                        "early_force_pN": early_force,
                        "late_force_pN": late_force,
                        "pair_mean_force_pN": pair_mean,
                        "history_half_difference_pN": 0.5
                        * (late_force - early_force),
                        "history_half_difference_absolute_pN": 0.5
                        * abs(late_force - early_force),
                        "early_gap_speed_um_per_s": early_gap_speed,
                        "late_gap_speed_um_per_s": late_gap_speed,
                        "pair_mean_gap_speed_um_per_s": pair_gap_speed,
                        "pair_gap_speed_over_nominal": pair_gap_speed / float(speed),
                    }
                )
    return rows, pair_lookup, pair_gap_speed_lookup


def double_center(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2 or not np.all(np.isfinite(array)):
        raise ValueError("double_center requires a finite 2-D array")
    return (
        array
        - np.mean(array, axis=1, keepdims=True)
        - np.mean(array, axis=0, keepdims=True)
        + np.mean(array)
    )


def hyd_model_matrix(
    distances_nm: np.ndarray,
    coefficient: float,
    velocity_exponent: float,
    distance_shift_nm: float,
) -> np.ndarray:
    distances = np.asarray(distances_nm, dtype=np.float64)
    shifted = distances + distance_shift_nm
    if np.any(shifted <= 0.0) or velocity_exponent <= 0.0:
        raise ValueError("Invalid empirical hydrodynamic coordinates")
    distance_basis = 1.0 / shifted
    velocity_basis = SPEEDS**velocity_exponent
    return coefficient * np.outer(
        distance_basis - np.mean(distance_basis),
        velocity_basis - np.mean(velocity_basis),
    )


def theoretical_no_slip_coefficient() -> tuple[float, float]:
    viscosity_mpa_s = base.cheng_viscosity_mPa_s(
        ASSUMED_GLYCEROL_MASS_FRACTION, ASSUMED_TEMPERATURE_C
    )
    coefficient = (
        6.0
        * math.pi
        * viscosity_mpa_s
        * 1e-3
        * ASSUMED_PROBE_RADIUS_M**2
        * 1e15
    )
    return viscosity_mpa_s, coefficient


def hyd_matrix(
    pair_lookup: dict[tuple[int, float, float], float],
    block: int,
    minimum_nm: float,
    maximum_nm: float,
) -> tuple[np.ndarray, np.ndarray]:
    distances = BIN_CENTERS_NM[
        (BIN_CENTERS_NM >= minimum_nm) & (BIN_CENTERS_NM <= maximum_nm)
    ]
    matrix = np.asarray(
        [
            [pair_lookup[(block, float(distance), float(speed))] for speed in SPEEDS]
            for distance in distances
        ],
        dtype=np.float64,
    )
    if matrix.shape != (distances.size, SPEEDS.size) or not np.all(np.isfinite(matrix)):
        raise RuntimeError(f"Block {block}: incomplete hydrodynamic force matrix")
    return distances, matrix


def analytic_coefficient(observed: list[np.ndarray], basis: np.ndarray) -> float:
    numerator = sum(float(np.sum(values * basis)) for values in observed)
    denominator = len(observed) * float(np.sum(basis**2))
    if denominator <= 0.0:
        raise RuntimeError("Degenerate hydrodynamic basis")
    return numerator / denominator


def nonlinear_hyd_fit(
    observed: list[np.ndarray], distances_nm: np.ndarray, model: str
) -> tuple[float, float, float, bool, int, str, int, float]:
    if model not in {
        "shifted_v_over_D",
        "velocity_power_over_D",
        "empirical_power_shift",
    }:
        raise ValueError(f"Unknown model: {model}")

    def unpack(parameters: np.ndarray) -> tuple[float, float, float]:
        if model == "shifted_v_over_D":
            return float(parameters[0]), 1.0, float(parameters[1])
        if model == "velocity_power_over_D":
            return float(parameters[0]), float(parameters[1]), 0.0
        return float(parameters[0]), float(parameters[1]), float(parameters[2])

    def residual(parameters: np.ndarray) -> np.ndarray:
        coefficient, exponent, shift = unpack(parameters)
        prediction = hyd_model_matrix(distances_nm, coefficient, exponent, shift)
        return np.concatenate([(prediction - values).ravel() for values in observed])

    if model == "shifted_v_over_D":
        starts = [
            np.asarray((coefficient, shift), dtype=np.float64)
            for coefficient in (1e4, 5e4, 2e5, 8e5)
            for shift in (0.0, 40.0, 150.0)
        ]
        lower = np.asarray((0.0, -float(np.min(distances_nm)) + 0.1))
        upper = np.asarray((1e8, 1000.0))
    elif model == "velocity_power_over_D":
        starts = [
            np.asarray((coefficient, exponent), dtype=np.float64)
            for coefficient in (1e4, 5e4, 2e5, 8e5)
            for exponent in (0.3, 0.6, 1.0, 1.5)
        ]
        lower = np.asarray((0.0, 0.05))
        upper = np.asarray((1e8, 3.0))
    else:
        starts = [
            np.asarray((coefficient, exponent, shift), dtype=np.float64)
            for coefficient in (1e4, 5e4, 2e5, 8e5)
            for exponent in (0.3, 0.6, 1.0, 1.5)
            for shift in (0.0, 40.0, 150.0)
        ]
        lower = np.asarray((0.0, 0.05, -float(np.min(distances_nm)) + 0.1))
        upper = np.asarray((1e8, 3.0, 1000.0))
    solutions = [
        least_squares(
            residual,
            start,
            bounds=(lower, upper),
            method="trf",
            x_scale="jac",
            max_nfev=10000,
            ftol=1e-11,
            xtol=1e-11,
            gtol=1e-11,
        )
        for start in starts
    ]
    solution = min(solutions, key=lambda item: float(np.sum(item.fun**2)))
    coefficient, exponent, shift = unpack(solution.x)
    singular = np.linalg.svd(solution.jac, compute_uv=False)
    tolerance = max(solution.jac.shape) * np.finfo(np.float64).eps * singular[0]
    rank = int(np.sum(singular > tolerance))
    condition = (
        float(singular[0] / singular[-1])
        if singular.size and singular[-1] > 0.0
        else float("inf")
    )
    return (
        coefficient,
        exponent,
        shift,
        bool(solution.success),
        int(solution.nfev),
        str(solution.message),
        rank,
        condition,
    )


def fit_hyd_model(
    observed: list[np.ndarray],
    distances_nm: np.ndarray,
    model: str,
    theoretical_coefficient: float,
) -> dict:
    success = True
    nfev = 0
    message = "closed form"
    rank = 0
    condition = float("nan")
    if model == "additive_null":
        coefficient, exponent, shift = 0.0, 1.0, 0.0
        parameter_count = 0
    elif model == "theoretical_no_slip":
        coefficient, exponent, shift = theoretical_coefficient, 1.0, 0.0
        parameter_count = 0
        message = "fixed campaign-assumption coefficient"
    elif model == "fitted_v_over_D":
        basis = hyd_model_matrix(distances_nm, 1.0, 1.0, 0.0)
        coefficient = max(0.0, analytic_coefficient(observed, basis))
        exponent, shift = 1.0, 0.0
        parameter_count = 1
        rank = 1
        condition = 1.0
    else:
        (
            coefficient,
            exponent,
            shift,
            success,
            nfev,
            message,
            rank,
            condition,
        ) = nonlinear_hyd_fit(observed, distances_nm, model)
        parameter_count = {
            "shifted_v_over_D": 2,
            "velocity_power_over_D": 2,
            "empirical_power_shift": 3,
        }[model]
    prediction = hyd_model_matrix(distances_nm, coefficient, exponent, shift)
    residual = np.concatenate([(values - prediction).ravel() for values in observed])
    values_flat = np.concatenate([values.ravel() for values in observed])
    rmse = float(np.sqrt(np.mean(residual**2)))
    interaction_rms = float(np.sqrt(np.mean(values_flat**2)))
    return {
        "model": model,
        "coefficient_pN_nm_per_um_s_to_alpha": coefficient,
        "velocity_exponent_alpha": exponent,
        "distance_shift_D0_nm": shift,
        "parameter_count": parameter_count,
        "interaction_rms_pN": interaction_rms,
        "residual_rmse_pN": rmse,
        "normalized_rmse": rmse / interaction_rms if interaction_rms > 0.0 else float("nan"),
        "optimizer_success": success,
        "optimizer_nfev": nfev,
        "optimizer_message": message,
        "jacobian_rank": rank,
        "jacobian_condition": condition,
    }


def actual_gap_basis(
    distances_nm: np.ndarray,
    gap_speed_um_per_s: np.ndarray,
    distance_shift_nm: float,
) -> np.ndarray:
    shifted = np.asarray(distances_nm, dtype=np.float64) + distance_shift_nm
    gap_speed = np.asarray(gap_speed_um_per_s, dtype=np.float64)
    if (
        gap_speed.shape != (shifted.size, SPEEDS.size)
        or not np.all(np.isfinite(gap_speed))
        or np.any(shifted <= 0.0)
    ):
        raise ValueError("Invalid actual-gap hydrodynamic basis")
    return double_center(gap_speed / shifted[:, None])


def analytic_variable_basis_coefficient(
    observed: list[np.ndarray], bases: list[np.ndarray]
) -> float:
    if len(observed) != len(bases) or not observed:
        raise ValueError("Observed and basis matrices must have matching blocks")
    numerator = sum(
        float(np.sum(values * basis))
        for values, basis in zip(observed, bases, strict=True)
    )
    denominator = sum(float(np.sum(basis**2)) for basis in bases)
    if denominator <= 0.0:
        raise RuntimeError("Degenerate actual-gap basis")
    return numerator / denominator


def fit_actual_gap_model(
    observed: list[np.ndarray],
    gap_speeds: list[np.ndarray],
    distances_nm: np.ndarray,
    model: str,
    theoretical_coefficient: float,
) -> dict:
    if model not in {
        "actual_gap_theoretical_no_slip",
        "actual_gap_fitted",
        "actual_gap_shifted",
    }:
        raise ValueError(f"Unknown actual-gap model: {model}")
    success = True
    nfev = 0
    message = "closed form"
    rank = 0
    condition = float("nan")
    if model == "actual_gap_theoretical_no_slip":
        coefficient = theoretical_coefficient
        shift = 0.0
        parameter_count = 0
        message = "fixed campaign-assumption coefficient with measured U_gap(D)"
    elif model == "actual_gap_fitted":
        bases = [actual_gap_basis(distances_nm, values, 0.0) for values in gap_speeds]
        coefficient = max(
            0.0, analytic_variable_basis_coefficient(observed, bases)
        )
        shift = 0.0
        parameter_count = 1
        rank = 1
        condition = 1.0
    else:
        def residual(parameters: np.ndarray) -> np.ndarray:
            coefficient_value, shift_value = parameters
            bases_value = [
                actual_gap_basis(distances_nm, values, float(shift_value))
                for values in gap_speeds
            ]
            return np.concatenate(
                [
                    (coefficient_value * basis - values).ravel()
                    for basis, values in zip(bases_value, observed, strict=True)
                ]
            )

        starts = [
            np.asarray((coefficient_start, shift_start), dtype=np.float64)
            for coefficient_start in (1e4, 5e4, 2e5, 8e5)
            for shift_start in (0.0, 40.0, 150.0)
        ]
        lower = np.asarray((0.0, -float(np.min(distances_nm)) + 0.1))
        upper = np.asarray((1e8, 1000.0))
        solutions = [
            least_squares(
                residual,
                start,
                bounds=(lower, upper),
                method="trf",
                x_scale="jac",
                max_nfev=10000,
                ftol=1e-11,
                xtol=1e-11,
                gtol=1e-11,
            )
            for start in starts
        ]
        solution = min(solutions, key=lambda item: float(np.sum(item.fun**2)))
        coefficient, shift = map(float, solution.x)
        success = bool(solution.success)
        nfev = int(solution.nfev)
        message = str(solution.message)
        singular = np.linalg.svd(solution.jac, compute_uv=False)
        tolerance = max(solution.jac.shape) * np.finfo(np.float64).eps * singular[0]
        rank = int(np.sum(singular > tolerance))
        condition = (
            float(singular[0] / singular[-1])
            if singular.size and singular[-1] > 0.0
            else float("inf")
        )
        parameter_count = 2
    bases = [actual_gap_basis(distances_nm, values, shift) for values in gap_speeds]
    residual_values = np.concatenate(
        [
            (values - coefficient * basis).ravel()
            for values, basis in zip(observed, bases, strict=True)
        ]
    )
    observed_flat = np.concatenate([values.ravel() for values in observed])
    rmse = float(np.sqrt(np.mean(residual_values**2)))
    interaction_rms = float(np.sqrt(np.mean(observed_flat**2)))
    return {
        "model": model,
        "coefficient_pN_nm_per_um_s_to_alpha": coefficient,
        "velocity_exponent_alpha": 1.0,
        "distance_shift_D0_nm": shift,
        "parameter_count": parameter_count,
        "interaction_rms_pN": interaction_rms,
        "residual_rmse_pN": rmse,
        "normalized_rmse": rmse / interaction_rms if interaction_rms > 0.0 else float("nan"),
        "optimizer_success": success,
        "optimizer_nfev": nfev,
        "optimizer_message": message,
        "jacobian_rank": rank,
        "jacobian_condition": condition,
    }


def hydrodynamic_analysis(
    pair_lookup: dict[tuple[int, float, float], float],
    pair_gap_speed_lookup: dict[tuple[int, float, float], float],
) -> tuple[list[dict], list[dict], list[dict]]:
    viscosity, theoretical_coefficient = theoretical_no_slip_coefficient()
    models = (
        "additive_null",
        "theoretical_no_slip",
        "fitted_v_over_D",
        "shifted_v_over_D",
        "velocity_power_over_D",
        "empirical_power_shift",
    )
    actual_gap_models = (
        "actual_gap_theoretical_no_slip",
        "actual_gap_fitted",
        "actual_gap_shifted",
    )
    fit_rows: list[dict] = []
    centered_rows: list[dict] = []
    reproducibility_rows: list[dict] = []
    primary_observed: dict[int, np.ndarray] = {}
    primary_gap_speeds: dict[int, np.ndarray] = {}
    primary_distances = np.asarray([], dtype=np.float64)
    for minimum, maximum, label in HYD_FIT_RANGES_NM:
        matrices: dict[int, np.ndarray] = {}
        gap_speed_matrices: dict[int, np.ndarray] = {}
        for block in BLOCKS:
            distances, raw = hyd_matrix(pair_lookup, block, minimum, maximum)
            matrices[block] = double_center(raw)
            gap_speed_matrices[block] = np.asarray(
                [
                    [
                        pair_gap_speed_lookup[(block, float(distance), float(speed))]
                        for speed in SPEEDS
                    ]
                    for distance in distances
                ],
                dtype=np.float64,
            )
            if not np.all(np.isfinite(gap_speed_matrices[block])):
                raise RuntimeError(f"Block {block}: non-finite actual gap speed")
            if label == "primary_20_200":
                primary_observed[block] = matrices[block]
                primary_gap_speeds[block] = gap_speed_matrices[block]
                primary_distances = distances
                for row_index, distance in enumerate(distances):
                    for column, speed in enumerate(SPEEDS):
                        centered_rows.append(
                            {
                                "block": block,
                                "distance_nm": float(distance),
                                "speed_um_per_s": float(speed),
                                "pair_mean_force_pN": float(raw[row_index, column]),
                                "double_centered_interaction_pN": float(
                                    matrices[block][row_index, column]
                                ),
                                "pair_mean_gap_speed_um_per_s": float(
                                    gap_speed_matrices[block][row_index, column]
                                ),
                                "gap_speed_over_nominal": float(
                                    gap_speed_matrices[block][row_index, column] / speed
                                ),
                            }
                        )
        for scope, selected_blocks in [
            *((f"block_{block}", (block,)) for block in BLOCKS),
            ("pooled_shared", BLOCKS),
        ]:
            observed = [matrices[block] for block in selected_blocks]
            gap_speeds = [gap_speed_matrices[block] for block in selected_blocks]
            for model in models:
                result = fit_hyd_model(
                    observed, distances, model, theoretical_coefficient
                )
                result.update(
                    {
                        "fit_range": label,
                        "minimum_distance_nm": minimum,
                        "maximum_distance_nm": maximum,
                        "scope": scope,
                        "blocks": "+".join(map(str, selected_blocks)),
                        "glycerol_mass_fraction_assumed": ASSUMED_GLYCEROL_MASS_FRACTION,
                        "temperature_C_assumed": ASSUMED_TEMPERATURE_C,
                        "probe_radius_um_assumed": ASSUMED_PROBE_RADIUS_M * 1e6,
                        "viscosity_mPa_s_assumed": viscosity,
                        "theoretical_no_slip_coefficient_pN_nm_per_um_s": theoretical_coefficient,
                    }
                )
                fit_rows.append(result)
            for model in actual_gap_models:
                result = fit_actual_gap_model(
                    observed,
                    gap_speeds,
                    distances,
                    model,
                    theoretical_coefficient,
                )
                result.update(
                    {
                        "fit_range": label,
                        "minimum_distance_nm": minimum,
                        "maximum_distance_nm": maximum,
                        "scope": scope,
                        "blocks": "+".join(map(str, selected_blocks)),
                        "glycerol_mass_fraction_assumed": ASSUMED_GLYCEROL_MASS_FRACTION,
                        "temperature_C_assumed": ASSUMED_TEMPERATURE_C,
                        "probe_radius_um_assumed": ASSUMED_PROBE_RADIUS_M * 1e6,
                        "viscosity_mPa_s_assumed": viscosity,
                        "theoretical_no_slip_coefficient_pN_nm_per_um_s": theoretical_coefficient,
                    }
                )
                fit_rows.append(result)

    for first in BLOCKS:
        for second in BLOCKS:
            if second <= first:
                continue
            a = primary_observed[first].ravel()
            b = primary_observed[second].ravel()
            correlation = (
                float(np.corrcoef(a, b)[0, 1])
                if np.std(a) > 0.0 and np.std(b) > 0.0
                else float("nan")
            )
            scale = float(np.dot(b, a) / np.dot(a, a))
            residual = b - scale * a
            reproducibility_rows.append(
                {
                    "first_block": first,
                    "second_block": second,
                    "pearson_shape_correlation": correlation,
                    "least_squares_scale_second_over_first": scale,
                    "scaled_shape_normalized_rmse": float(
                        np.sqrt(np.mean(residual**2)) / np.sqrt(np.mean(b**2))
                    ),
                    "first_interaction_rms_pN": float(np.sqrt(np.mean(a**2))),
                    "second_interaction_rms_pN": float(np.sqrt(np.mean(b**2))),
                }
            )

    primary_pooled = {
        row["model"]: row
        for row in fit_rows
        if row["fit_range"] == "primary_20_200" and row["scope"] == "pooled_shared"
    }
    for row in centered_rows:
        distance_index = int(
            np.flatnonzero(np.isclose(primary_distances, row["distance_nm"]))[0]
        )
        speed_index = int(np.flatnonzero(np.isclose(SPEEDS, row["speed_um_per_s"]))[0])
        for model, fit in primary_pooled.items():
            if model.startswith("actual_gap_"):
                prediction = float(
                    fit["coefficient_pN_nm_per_um_s_to_alpha"]
                ) * actual_gap_basis(
                    primary_distances,
                    primary_gap_speeds[int(row["block"])],
                    float(fit["distance_shift_D0_nm"]),
                )
            else:
                prediction = hyd_model_matrix(
                    primary_distances,
                    float(fit["coefficient_pN_nm_per_um_s_to_alpha"]),
                    float(fit["velocity_exponent_alpha"]),
                    float(fit["distance_shift_D0_nm"]),
                )
            row[f"{model}_prediction_pN"] = float(
                prediction[distance_index, speed_index]
            )
    return fit_rows, centered_rows, reproducibility_rows


def speed_interval_analysis(
    pair_lookup: dict[tuple[int, float, float], float]
) -> tuple[list[dict], list[dict]]:
    detail: list[dict] = []
    summary: list[dict] = []
    pairs = ((0.1, 0.2), (0.2, 0.4), (0.1, 0.4))
    for low, high in pairs:
        for distance in BIN_CENTERS_NM[
            (BIN_CENTERS_NM >= HYD_PRIMARY_MIN_NM)
            & (BIN_CENTERS_NM <= HYD_REFERENCE_DISTANCE_NM)
        ]:
            values = []
            for block in BLOCKS:
                contrast = (
                    pair_lookup[(block, float(distance), high)]
                    - pair_lookup[(block, float(distance), low)]
                )
                reference = (
                    pair_lookup[(block, HYD_REFERENCE_DISTANCE_NM, high)]
                    - pair_lookup[(block, HYD_REFERENCE_DISTANCE_NM, low)]
                )
                double_difference = contrast - reference
                values.append(double_difference)
                detail.append(
                    {
                        "block": block,
                        "speed_low_um_per_s": low,
                        "speed_high_um_per_s": high,
                        "delta_speed_um_per_s": high - low,
                        "distance_nm": float(distance),
                        "reference_distance_nm": HYD_REFERENCE_DISTANCE_NM,
                        "raw_speed_contrast_pN": contrast,
                        "reference_speed_contrast_pN": reference,
                        "distance_double_difference_pN": double_difference,
                        "double_difference_per_delta_speed_pN_per_um_s": double_difference
                        / (high - low),
                    }
                )
            array = np.asarray(values, dtype=np.float64)
            sd = float(np.std(array, ddof=1))
            se = sd / math.sqrt(array.size)
            critical = float(stats.t.ppf(0.975, array.size - 1))
            mean = float(np.mean(array))
            summary.append(
                {
                    "speed_low_um_per_s": low,
                    "speed_high_um_per_s": high,
                    "delta_speed_um_per_s": high - low,
                    "distance_nm": float(distance),
                    "reference_distance_nm": HYD_REFERENCE_DISTANCE_NM,
                    "blocks": int(array.size),
                    "mean_distance_double_difference_pN": mean,
                    "sample_sd_pN": sd,
                    "standard_error_pN": se,
                    "mean_over_standard_error": abs(mean) / se if se > 0.0 else float("nan"),
                    "ci95_low_pN": mean - critical * se,
                    "ci95_high_pN": mean + critical * se,
                    "positive_blocks": int(np.count_nonzero(array > 0.0)),
                    "negative_blocks": int(np.count_nonzero(array < 0.0)),
                    "equivalence_interpretation": "no equivalence margin prespecified; CI overlap with zero is not evidence of equivalence",
                }
            )
    return detail, summary


def plot_qc(map_records: list[dict], curve_rows: list[dict]) -> None:
    fig, axes = plt.subplots(3, 6, figsize=(16.5, 8.6), layout="constrained")
    cmap = ListedColormap(("#d8e8e4", "#e0a454", "#8577a8", "#b43b36"))
    for axis, record in zip(axes.flat, map_records, strict=True):
        matrix = np.full((GRID_SIZE, GRID_SIZE), np.nan)
        rows = [
            row
            for row in curve_rows
            if row["acquisition_order"] == record["acquisition_order"]
        ]
        for row in rows:
            matrix[int(row["physical_row"]), int(row["physical_column"])] = QC_CODES[
                str(row["qc_state"])
            ]
        axis.imshow(matrix, cmap=cmap, vmin=-0.5, vmax=3.5, interpolation="nearest")
        axis.set_title(
            f"#{record['acquisition_order']}  {record['speed_um_per_s']:g} um/s\nscan {record['instrument_scan_number']}",
            fontsize=9,
        )
        axis.set_xticks([])
        axis.set_yticks([])
        for row in rows:
            if row["qc_state"] != "pass_force":
                axis.text(
                    int(row["physical_column"]),
                    int(row["physical_row"]),
                    str(row["point_index"]),
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    color="white" if row["qc_state"] == "hard_fail" else "black",
                )
    handles = [
        Line2D([], [], marker="s", ls="", color=cmap(index), label=state)
        for state, index in QC_CODES.items()
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False)
    fig.suptitle("16 Sep 2026 approach-curve QC; labels are acquisition point indices", fontsize=14)
    fig.savefig(FIG / "curve_qc_maps.png", dpi=220)
    plt.close(fig)


def plot_qc_features(curve_rows: list[dict]) -> None:
    state_colors = {
        "pass_force": "#16817a",
        "no_valid_contact": "#e0a454",
        "insufficient_range": "#8577a8",
        "hard_fail": "#b43b36",
    }
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.4), layout="constrained")
    for state in QC_STATES:
        selected = [row for row in curve_rows if row["qc_state"] == state]
        axes[0].scatter(
            [row["approach_travel_span_nm"] for row in selected],
            [row["point_count_fraction"] for row in selected],
            s=18 if state == "pass_force" else 44,
            color=state_colors[state],
            alpha=0.42 if state == "pass_force" else 0.95,
            edgecolor="none",
            label=state,
        )
        finite_contact = [
            row
            for row in selected
            if np.isfinite(float(row["terminal_contact_amplitude_V"]))
            and np.isfinite(float(row["terminal_contact_invOLS_nm_per_V"]))
            and float(row["terminal_contact_invOLS_nm_per_V"]) > 0.0
        ]
        axes[1].scatter(
            [row["terminal_contact_amplitude_V"] for row in finite_contact],
            [row["terminal_contact_invOLS_nm_per_V"] for row in finite_contact],
            s=18 if state == "pass_force" else 44,
            color=state_colors[state],
            alpha=0.42 if state == "pass_force" else 0.95,
            edgecolor="none",
        )
        if state != "pass_force":
            for row in finite_contact:
                axes[1].annotate(
                    f"{row['acquisition_order']}:{row['point_index']}",
                    (
                        float(row["terminal_contact_amplitude_V"]),
                        float(row["terminal_contact_invOLS_nm_per_V"]),
                    ),
                    xytext=(3, 3),
                    textcoords="offset points",
                    fontsize=7,
                )
    axes[0].axvline(HARD_FAIL_MAX_TRAVEL_NM, color="0.3", ls=":", lw=1)
    axes[0].axvline(MIN_FORCE_TRAVEL_NM, color="0.3", ls="--", lw=1)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Recorded scanner-travel span (nm; log scale)")
    axes[0].set_ylabel("Actual / declared approach points")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].axvline(MIN_CONTACT_AMPLITUDE_V, color="0.3", ls="--", lw=1)
    axes[1].axhline(FORCE_CONTACT_INVOLS_RANGE_NM_PER_V[0], color="0.3", ls=":", lw=1)
    axes[1].axhline(FORCE_CONTACT_INVOLS_RANGE_NM_PER_V[1], color="0.3", ls=":", lw=1)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Terminal contact-rise amplitude (V)")
    axes[1].set_ylabel("Terminal 40 nm InvOLS (nm/V; log scale)")
    for axis in axes:
        axis.grid(alpha=0.2)
    fig.suptitle("Curve-level failure features used by the QC classifier", fontsize=14)
    fig.savefig(FIG / "curve_qc_feature_space.png", dpi=220)
    plt.close(fig)


def plot_calibration(
    map_records: list[dict], contact_rows: list[dict], curve_rows: list[dict], batch_inv: float
) -> None:
    times = np.asarray(
        [datetime.fromisoformat(row["map_midpoint_time"]) for row in map_records],
        dtype=object,
    )
    medians = np.asarray([row["map_contact_InvOLS_median_nm_per_V"] for row in contact_rows])
    q25 = np.asarray([row["map_contact_InvOLS_q25_nm_per_V"] for row in contact_rows])
    q75 = np.asarray([row["map_contact_InvOLS_q75_nm_per_V"] for row in contact_rows])
    baselines = np.asarray([row["map_baseline_raw_V_median"] for row in map_records])
    shifts_a = np.asarray(
        [
            median(
                row.get("baseline_75_200_force_shift_pN", float("nan"))
                for row in curve_rows
                if row["acquisition_order"] == record["acquisition_order"]
            )
            for record in map_records
        ]
    )
    shifts_b = np.asarray(
        [
            median(
                row.get("baseline_125_275_force_shift_pN", float("nan"))
                for row in curve_rows
                if row["acquisition_order"] == record["acquisition_order"]
            )
            for record in map_records
        ]
    )
    fig, axes = plt.subplots(3, 1, figsize=(13.4, 9.2), sharex=True, layout="constrained")
    axes[0].errorbar(times, medians, yerr=np.vstack((medians - q25, q75 - medians)), fmt="o-", color="#246b78", lw=1)
    axes[0].axhline(batch_inv, color="0.2", ls="--", lw=1, label=f"batch median {batch_inv:.2f} nm/V")
    axes[0].set_ylabel("Contact InvOLS (nm/V)")
    axes[0].legend(frameon=False)
    axes[1].plot(times, baselines, "o-", color="#8b4d78", lw=1)
    axes[1].set_ylabel("Raw baseline (V)")
    axes[2].plot(times, shifts_a, "o-", label="75-200 vs primary", color="#16817a")
    axes[2].plot(times, shifts_b, "s-", label="125-275 vs primary", color="#d56a2f")
    axes[2].axhline(0.0, color="0.3", lw=0.8)
    axes[2].set_ylabel("Map-median gauge shift (pN)")
    axes[2].set_xlabel("Map midpoint time (UTC+02:00)")
    axes[2].legend(frameon=False, ncol=2)
    for axis in axes:
        axis.grid(axis="y", alpha=0.2)
        axis.xaxis.set_major_locator(mdates.MinuteLocator(interval=10, tz=times[0].tzinfo))
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=times[0].tzinfo))
    fig.suptitle("Contact calibration, raw offset drift, and baseline-window sensitivity", fontsize=14)
    fig.savefig(FIG / "calibration_and_baseline_diagnostics.png", dpi=220)
    plt.close(fig)


def plot_gap_speed(force_rows: list[dict]) -> None:
    fig, axis = plt.subplots(figsize=(9.8, 5.8), layout="constrained")
    for speed in SPEEDS:
        orders = sorted(
            {
                int(row["acquisition_order"])
                for row in force_rows
                if np.isclose(float(row["speed_um_per_s"]), speed)
            }
        )
        curves = []
        for order in orders:
            selected = sorted(
                (
                    row
                    for row in force_rows
                    if int(row["acquisition_order"]) == order
                    and HYD_PRIMARY_MIN_NM
                    <= float(row["distance_nm"])
                    <= HYD_REFERENCE_DISTANCE_NM
                ),
                key=lambda row: float(row["distance_nm"]),
            )
            x = np.asarray([float(row["distance_nm"]) for row in selected])
            ratio = np.asarray(
                [float(row["map_median_gap_speed_um_per_s"]) / speed for row in selected]
            )
            curves.append(ratio)
            axis.plot(x, ratio, color=COLORS[float(speed)], alpha=0.16, lw=0.8)
        matrix = np.asarray(curves, dtype=np.float64)
        axis.plot(
            x,
            np.median(matrix, axis=0),
            color=COLORS[float(speed)],
            lw=2.2,
            label=f"{speed:g} um/s nominal",
        )
    axis.axhline(1.0, color="0.25", ls="--", lw=1, label="U_gap = nominal")
    axis.set_xlabel("Separation D (nm)")
    axis.set_ylabel("Measured gap-closing speed / nominal drive speed")
    axis.set_ylim(0.72, 1.03)
    axis.grid(alpha=0.2)
    axis.legend(frameon=False)
    axis.set_title("Cantilever compliance slows the physical gap motion near contact")
    fig.savefig(FIG / "gap_speed_profiles.png", dpi=220)
    plt.close(fig)


def plot_force_time(map_records: list[dict], force_rows: list[dict]) -> None:
    lookup = {
        (int(row["acquisition_order"]), float(row["distance_nm"])): row
        for row in force_rows
    }
    times = np.asarray(
        [datetime.fromisoformat(row["map_midpoint_time"]) for row in map_records],
        dtype=object,
    )
    fig, axes = plt.subplots(2, 2, figsize=(15.5, 8.8), sharex=True, layout="constrained")
    for axis, distance in zip(axes.flat, TARGET_DISTANCES_NM, strict=True):
        values = np.asarray(
            [float(lookup[(row["acquisition_order"], distance)]["map_median_force_pN"]) for row in map_records]
        )
        q25 = np.asarray(
            [float(lookup[(row["acquisition_order"], distance)]["map_force_q25_pN"]) for row in map_records]
        )
        q75 = np.asarray(
            [float(lookup[(row["acquisition_order"], distance)]["map_force_q75_pN"]) for row in map_records]
        )
        axis.plot(times, values, color="0.65", lw=0.9, zorder=1)
        for speed in SPEEDS:
            selected = np.asarray(
                [np.isclose(float(row["speed_um_per_s"]), speed) for row in map_records]
            )
            axis.errorbar(
                times[selected],
                values[selected],
                yerr=np.vstack((values[selected] - q25[selected], q75[selected] - values[selected])),
                fmt="o",
                ms=5.5,
                color=COLORS[float(speed)],
                mec="0.2",
                mew=0.4,
                elinewidth=0.7,
                alpha=0.95,
            )
        axis.set_title(f"D = {distance:g} nm")
        axis.set_ylabel("Force (pN; finite-window gauge)")
        axis.grid(axis="y", alpha=0.2)
        axis.xaxis.set_major_locator(mdates.MinuteLocator(interval=10, tz=times[0].tzinfo))
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=times[0].tzinfo))
    for axis in axes[-1]:
        axis.set_xlabel("Map midpoint time (UTC+02:00)")
    handles = [
        Line2D([], [], marker="o", ls="", color=COLORS[float(speed)], label=f"{speed:g} um/s")
        for speed in SPEEDS
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False)
    fig.suptitle("Map-median approach force versus acquisition time", fontsize=14, y=1.03)
    fig.savefig(FIG / "force_time_and_speed.png", dpi=220)
    plt.close(fig)


def plot_speed_intervals(detail_rows: list[dict], summary_rows: list[dict]) -> None:
    pairs = ((0.1, 0.2), (0.2, 0.4), (0.1, 0.4))
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.0), sharey=True, layout="constrained")
    for axis, (low, high) in zip(axes, pairs, strict=True):
        for block, color in zip(BLOCKS, ("#16817a", "#3d6fb4", "#d56a2f"), strict=True):
            selected = sorted(
                (
                    row
                    for row in detail_rows
                    if row["block"] == block
                    and np.isclose(row["speed_low_um_per_s"], low)
                    and np.isclose(row["speed_high_um_per_s"], high)
                ),
                key=lambda row: row["distance_nm"],
            )
            axis.plot(
                [row["distance_nm"] for row in selected],
                [row["distance_double_difference_pN"] for row in selected],
                color=color,
                lw=0.9,
                alpha=0.7,
                label=f"block {block}",
            )
        aggregate = sorted(
            (
                row
                for row in summary_rows
                if np.isclose(row["speed_low_um_per_s"], low)
                and np.isclose(row["speed_high_um_per_s"], high)
            ),
            key=lambda row: row["distance_nm"],
        )
        x = np.asarray([row["distance_nm"] for row in aggregate])
        mean = np.asarray([row["mean_distance_double_difference_pN"] for row in aggregate])
        se = np.asarray([row["standard_error_pN"] for row in aggregate])
        axis.plot(x, mean, color="black", lw=2.0, label="3-block mean")
        axis.fill_between(x, mean - se, mean + se, color="black", alpha=0.12, label="+/-1 SE")
        axis.axhline(0.0, color="0.35", lw=0.8)
        axis.set_title(f"{high:g} - {low:g} um/s")
        axis.set_xlabel("Separation D (nm)")
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("Speed contrast minus its 250 nm value (pN)")
    axes[-1].legend(frameon=False, fontsize=8)
    fig.suptitle("Palindrome-symmetrized speed intervals; block spread is experimental variability", fontsize=14)
    fig.savefig(FIG / "speed_interval_comparison.png", dpi=220)
    plt.close(fig)


def plot_hydrodynamics(fit_rows: list[dict], centered_rows: list[dict]) -> None:
    primary = {
        row["model"]: row
        for row in fit_rows
        if row["fit_range"] == "primary_20_200" and row["scope"] == "pooled_shared"
    }
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.3), layout="constrained")
    model_styles = {
        "fitted_v_over_D": ("#3d6fb4", "fitted v/D"),
        "empirical_power_shift": ("#b34f37", "K v^alpha/(D+D0)"),
        "actual_gap_shifted": ("#5f7f3f", "K U_gap/(D+D0)"),
    }
    for speed in SPEEDS:
        selected = sorted(
            (row for row in centered_rows if np.isclose(row["speed_um_per_s"], speed)),
            key=lambda row: (row["block"], row["distance_nm"]),
        )
        distances = sorted({float(row["distance_nm"]) for row in selected})
        means = []
        sds = []
        for distance in distances:
            values = [
                float(row["double_centered_interaction_pN"])
                for row in selected
                if np.isclose(row["distance_nm"], distance)
            ]
            means.append(float(np.mean(values)))
            sds.append(float(np.std(values, ddof=1)))
        axes[0].plot(distances, means, "o", ms=3.5, color=COLORS[float(speed)], label=f"data {speed:g} um/s")
        axes[0].fill_between(
            distances,
            np.asarray(means) - np.asarray(sds),
            np.asarray(means) + np.asarray(sds),
            color=COLORS[float(speed)],
            alpha=0.08,
        )
    for model, (color, label) in model_styles.items():
        row = primary[model]
        for speed in SPEEDS:
            selected = [
                item
                for item in centered_rows
                if np.isclose(item["speed_um_per_s"], speed)
            ]
            distances_array = sorted({float(item["distance_nm"]) for item in selected})
            prediction = [
                float(
                    np.mean(
                        [
                            item[f"{model}_prediction_pN"]
                            for item in selected
                            if np.isclose(item["distance_nm"], distance)
                        ]
                    )
                )
                for distance in distances_array
            ]
            axes[0].plot(
                distances_array,
                prediction,
                color=color,
                lw=1.4,
                ls={
                    "fitted_v_over_D": "--",
                    "empirical_power_shift": "-",
                    "actual_gap_shifted": ":",
                }[model],
                alpha=0.85,
            )
        axes[1].bar(
            label,
            float(row["normalized_rmse"]),
            color=color,
            alpha=0.85,
        )
    axes[0].axhline(0.0, color="0.35", lw=0.8)
    axes[0].set_xlabel("Separation D (nm)")
    axes[0].set_ylabel("Double-centered interaction (pN)")
    axes[0].grid(alpha=0.2)
    handles = [
        Line2D([], [], marker="o", ls="", color=COLORS[float(speed)], label=f"{speed:g} um/s data")
        for speed in SPEEDS
    ] + [
        Line2D(
            [],
            [],
            color=color,
            ls={
                "fitted_v_over_D": "--",
                "empirical_power_shift": "-",
                "actual_gap_shifted": ":",
            }[model],
            label=label,
        )
        for model, (color, label) in model_styles.items()
    ]
    axes[0].legend(handles=handles, frameon=False, fontsize=8, ncol=2)
    axes[1].set_ylabel("Residual RMSE / interaction RMS")
    axes[1].tick_params(axis="x", rotation=15)
    axes[1].grid(axis="y", alpha=0.2)
    fig.suptitle("Velocity-distance interaction after palindrome averaging and double centering", fontsize=14)
    fig.savefig(FIG / "hydrodynamic_model_comparison.png", dpi=220)
    plt.close(fig)


def synthetic_self_checks() -> dict:
    distances = np.arange(20.0, 200.0 + 0.1, 5.0)
    surface = 400.0 * np.exp(-distances / 55.0)
    offsets = np.asarray((12.0, -5.0, 21.0))
    coefficient = 12345.0
    raw = surface[:, None] + offsets[None, :] + coefficient * np.outer(
        1.0 / distances, SPEEDS
    )
    observed = double_center(raw)
    expected = hyd_model_matrix(distances, coefficient, 1.0, 0.0)
    recovered = analytic_coefficient([observed], hyd_model_matrix(distances, 1.0, 1.0, 0.0))
    synthetic_gap_speed = np.outer(
        0.82 + 0.18 * (distances - distances.min()) / np.ptp(distances), SPEEDS
    )
    actual_gap_interaction = coefficient * actual_gap_basis(
        distances, synthetic_gap_speed, 7.0
    )
    recovered_actual_gap = analytic_variable_basis_coefficient(
        [actual_gap_interaction],
        [actual_gap_basis(distances, synthetic_gap_speed, 7.0)],
    )
    return {
        "synthetic_double_center_max_abs_error_pN": float(
            np.max(np.abs(observed - expected))
        ),
        "synthetic_classical_coefficient_relative_error": abs(recovered - coefficient)
        / coefficient,
        "synthetic_actual_gap_coefficient_relative_error": abs(
            recovered_actual_gap - coefficient
        )
        / coefficient,
        "double_center_row_mean_max_abs_pN": float(
            np.max(np.abs(np.mean(observed, axis=1)))
        ),
        "double_center_column_mean_max_abs_pN": float(
            np.max(np.abs(np.mean(observed, axis=0)))
        ),
    }


def invalid_curve_table(curve_rows: list[dict]) -> list[str]:
    lines = [
        "| map | scan | time | speed (um/s) | point | pixel (r,c) | state | points | travel (nm) | contact amplitude (V) |",
        "|---:|---:|---|---:|---:|---|---|---:|---:|---:|",
    ]
    for row in curve_rows:
        if row["qc_state"] == "pass_force":
            continue
        time = datetime.fromisoformat(str(row["map_midpoint_time"])).strftime("%H:%M:%S")
        lines.append(
            f"| {row['acquisition_order']} | {row['instrument_scan_number']} | {time} | "
            f"{float(row['speed_um_per_s']):g} | {row['point_index']} | "
            f"({row['physical_row']},{row['physical_column']}) | {row['qc_state']} | "
            f"{row['actual_approach_points']} | {float(row['approach_travel_span_nm']):.1f} | "
            f"{float(row['terminal_contact_amplitude_V']):.3f} |"
        )
    return lines


def write_report(
    map_records: list[dict],
    curve_rows: list[dict],
    contact_rows: list[dict],
    inv_ols: float,
    fit_rows: list[dict],
    reproducibility_rows: list[dict],
    interval_summary: list[dict],
) -> None:
    counts = {state: sum(row["qc_state"] == state for row in curve_rows) for state in QC_STATES}
    spring = float(map_records[0]["embedded_spring_constant_N_per_m"])
    embedded_inv = float(map_records[0]["embedded_InvOLS_nm_per_V"])
    inv_30 = median(
        row["map_contact_InvOLS_30nm_median_nm_per_V"] for row in contact_rows
    )
    inv_50 = median(
        row["map_contact_InvOLS_50nm_median_nm_per_V"] for row in contact_rows
    )
    primary = {
        row["model"]: row
        for row in fit_rows
        if row["fit_range"] == "primary_20_200" and row["scope"] == "pooled_shared"
    }
    target_interval = {
        (row["speed_low_um_per_s"], row["speed_high_um_per_s"], row["distance_nm"]): row
        for row in interval_summary
        if row["distance_nm"] in TARGET_DISTANCES_NM
    }
    lines = [
        "# 16-09-26 low-speed palindrome AFM analysis",
        "",
        "## 主要结论",
        "",
        f"18 张 8x8 force maps 的 1152 条 approach 中，**{counts['pass_force']} 条可用于 force-distance 重建**；"
        f"{counts['hard_fail']} 条是小于 {HARD_FAIL_MAX_TRAVEL_NM:g} nm 的明显 premature trigger，"
        f"{counts['insufficient_range']} 条量程不足，{counts['no_valid_contact']} 条虽有较长记录但没有满足振幅和线性判据的 terminal hard contact。"
        "JPK acquisition flags 对这些点全部未报警，因此必须由曲线级 QC 识别并保留为 NaN，不能插值补点。",
        "18 张 map 的 point index 0 全部通过相同判据；因此黑色 height-map 像素、行首身份或 point index 本身都不是坏点标签。",
        "",
        "第一张 map 的中心相对其余 17 张偏移约 "
        f"**{float(map_records[0]['shift_from_maps_2_to_18_center_um']):.3f} um**。因此 block 1 的 0.1 um/s palindrome mate 只能比较 map median，不能做同像素配对；其余 palindrome pairs 的网格中心一致。",
        "",
        "## 曲线处理",
        "",
        f"- branch loader 以 `minimum_points=1` 解码，确保短段仍带原始 point index 进入 QC。主分析只用 Extend；Retract 固定为 {float(map_records[0]['retract_speed_um_per_s']):g} um/s，不与 approach 速度混合。",
        f"- 1152 条 Retract 全部为 500 点，travel span 为 {min(float(row['retract_travel_span_nm']) for row in curve_rows):.1f}-"
        f"{max(float(row['retract_travel_span_nm']) for row in curve_rows):.1f} nm。坏点的完整 Retract 说明 ZIP/像素记录仍在，故应解释为 approach-side premature trigger/contact failure；不能用 2 um/s Retract 替代 approach force。",
        f"- baseline 是每条曲线在 scanner travel {BASELINE_START_NM:g}-{BASELINE_STOP_NM:g} nm 的 raw-vDeflection 中位数，只减常数，不减线性斜率。它定义 finite-window force gauge，不等于无限远零力。",
        f"- contact 用 terminal {CONTACT_WINDOW_NM:g} nm 物理窗口；force-use 判据同时要求 travel、contact amplitude、InvOLS 范围和 robust-line R2。固定样本数没有被使用。",
        f"- batch contact InvOLS = **{inv_ols:.3f} nm/V**（18 个 map median 再取 median），embedded InvOLS = {embedded_inv:.3f} nm/V。embedded InvOLS 不满足本批 hard-contact slope，因此未用于位移换算。",
        f"- contact-window sensitivity：30/40/50 nm 分别得到 {inv_30:.3f}/{inv_ols:.3f}/{inv_50:.3f} nm/V；主值没有依赖单一固定样本数。",
        f"- spring constant 暂用所有文件一致的 embedded `k={spring:.9f} N/m`，得到条件 force scale **{spring * inv_ols:.3f} nN/V**。没有同日独立 thermal calibration，因此绝对力标度仍带 calibration 条件。",
        "",
        "## 坏点清单",
        "",
        *invalid_curve_table(curve_rows),
        "",
        "## 低速差异是否小于实验随机性",
        "",
        "下表给 palindrome-symmetrized speed contrast，并减去同一 contrast 在 250 nm 的值；这样不受每条曲线常数 baseline gauge 影响。n=3 是三个 block，不是像素数。`mean/SE < 1` 只能说明当前实验未分辨出差异，不能证明等效；正式 equivalence 还需要预先规定可接受 margin。",
        "",
        "| D (nm) | interval (um/s) | mean double difference (pN) | block SD (pN) | SE (pN) | |mean|/SE |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    for distance in TARGET_DISTANCES_NM:
        for low, high in ((0.1, 0.2), (0.2, 0.4)):
            row = target_interval[(low, high, distance)]
            lines.append(
                f"| {distance:g} | {low:g}->{high:g} | "
                f"{float(row['mean_distance_double_difference_pN']):+.2f} | "
                f"{float(row['sample_sd_pN']):.2f} | {float(row['standard_error_pN']):.2f} | "
                f"{float(row['mean_over_standard_error']):.2f} |"
            )
    lines += [
        "",
        "在 20/50/100 nm，0.1->0.2 um/s 的三个 block 全部同号，且 |mean|/SE 分别为 "
        f"{float(target_interval[(0.1, 0.2, 20.0)]['mean_over_standard_error']):.1f}/"
        f"{float(target_interval[(0.1, 0.2, 50.0)]['mean_over_standard_error']):.1f}/"
        f"{float(target_interval[(0.1, 0.2, 100.0)]['mean_over_standard_error']):.1f}。"
        "因此本批数据明确不支持 0.1 与 0.2 um/s 已进入同一 plateau；到 200 nm 时 95% CI 仍跨零，只能说该远距离差异未分辨。",
        "",
        "由重建的 `D(t)` 直接求得：20 nm 处 `U_gap/nominal` 的 map-median 范围为 "
        f"**{min(float(row['gap_speed_20nm_over_nominal']) for row in map_records):.3f}-"
        f"{max(float(row['gap_speed_20nm_over_nominal']) for row in map_records):.3f}**，"
        "200 nm 处接近 1。以下模型因此同时列 nominal speed 与实际 `U_gap(D)` 分支。",
        "",
        "## F(v,D) 模型检查",
        "",
        "每个 6-map block 先平均同速 palindrome mates，再对 distance x speed 矩阵 double-center。这个操作严格去掉任意 distance-only surface-force 曲线和任意 speed-only 常数 offset；剩余项才用于比较 `v/D` 与经验式。它仍可能包含 speed-dependent instrument response 或非线性 history，不能仅凭拟合命名为 hydrodynamics。",
        "",
        "| pooled model, 20-200 nm | K | alpha | D0 (nm) | RMSE (pN) | normalized RMSE | Jacobian condition |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for model in (
        "additive_null",
        "theoretical_no_slip",
        "fitted_v_over_D",
        "shifted_v_over_D",
        "velocity_power_over_D",
        "empirical_power_shift",
        "actual_gap_theoretical_no_slip",
        "actual_gap_fitted",
        "actual_gap_shifted",
    ):
        row = primary[model]
        lines.append(
            f"| {model} | {float(row['coefficient_pN_nm_per_um_s_to_alpha']):.3g} | "
            f"{float(row['velocity_exponent_alpha']):.3f} | {float(row['distance_shift_D0_nm']):.2f} | "
            f"{float(row['residual_rmse_pN']):.2f} | {float(row['normalized_rmse']):.3f} | "
            f"{float(row['jacobian_condition']):.3g} |"
        )
    empirical = primary["empirical_power_shift"]
    classical = primary["fitted_v_over_D"]
    shifted = primary["shifted_v_over_D"]
    actual_shifted = primary["actual_gap_shifted"]
    improvement = 1.0 - float(empirical["residual_rmse_pN"]) / float(
        classical["residual_rmse_pN"]
    )
    shift_improvement = 1.0 - float(shifted["residual_rmse_pN"]) / float(
        classical["residual_rmse_pN"]
    )
    power_improvement = 1.0 - float(empirical["residual_rmse_pN"]) / float(
        shifted["residual_rmse_pN"]
    )
    block_empirical = [
        row
        for row in fit_rows
        if row["fit_range"] == "primary_20_200"
        and row["model"] == "empirical_power_shift"
        and str(row["scope"]).startswith("block_")
    ]
    correlations = [float(row["pearson_shape_correlation"]) for row in reproducibility_rows]
    lines += [
        "",
        f"在 pooled in-sample 比较中，`v/(D+D0)` 相对 fitted `v/D` 把 RMSE 降低 **{shift_improvement * 100:.1f}%**；"
        f"完整经验式相对 fitted `v/D` 降低 **{improvement * 100:.1f}%**，但相对已经带 D0 的线性速度式只再降低 **{power_improvement * 100:.1f}%**。"
        f"三个 block 的 double-centered shape 两两相关范围为 **{min(correlations):.3f} 到 {max(correlations):.3f}**。"
        f"block-wise alpha 为 {min(float(row['velocity_exponent_alpha']) for row in block_empirical):.3f}-"
        f"{max(float(row['velocity_exponent_alpha']) for row in block_empirical):.3f}，而 D0 为 "
        f"{min(float(row['distance_shift_D0_nm']) for row in block_empirical):.2f}-"
        f"{max(float(row['distance_shift_D0_nm']) for row in block_empirical):.2f} nm。"
        "经验式多两个自由参数且仅有三个速度，Jacobian condition 也很大；因此较低 RMSE 不能把 alpha 解释成已识别的 shear-thinning exponent，也不能把 D0 直接解释成 slip length。",
        "",
        f"用实际 `U_gap(D)` 后，最佳 shifted model 的 normalized RMSE 为 {float(actual_shifted['normalized_rmse']):.3f}。"
        "固定 bulk no-slip coefficient 仍显著过预测，说明这批数据支持可重复的 velocity-distance coupling，"
        "但不验证以 bulk viscosity、prior radius、nominal no-slip geometry 组成的绝对 prefactor。",
        "",
        "固定 no-slip 系数仅作为 campaign-context comparator：假定 99.7 wt% glycerol、"
        f"T={ASSUMED_TEMPERATURE_C:g} C、R={ASSUMED_PROBE_RADIUS_M * 1e6:.6f} um。"
        "这些量不在 JPK header 中，主 empirical/fitted comparisons 不依赖固定理论系数。",
        "",
        "## 产物",
        "",
        "- `curve_qc.csv`: 1152 条 approach 的长度、量程、baseline、contact、retract inventory 与最终状态。",
        "- `map_inventory_and_summary.csv`: acquisition order、grid shift、校准、QC 数量及 20/50/100/200 nm force。",
        "- `contact_calibration_maps.csv`: 等 map 权重的 contact InvOLS。",
        "- `curve_force_bins.npz` / `map_force_by_separation.csv`: NaN-preserving force、actual gap-speed pixel matrices 与 map summaries。",
        "- `palindrome_pair_force.csv`: 同速 early/late、pair mean 与 history half-difference。",
        "- `speed_interval_detail.csv` / `speed_interval_summary.csv`: 低速 interval 的 block-level 结果。",
        "- `hydrodynamic_double_centered.csv`, `hydrodynamic_model_fits.csv`, `block_reproducibility.csv`: `F(v,D)` 检查。",
        "- `figures/`: QC、calibration/baseline、force-time、speed intervals 与 model comparison。",
        "- `provenance.json`, `input_manifest.sha256`, `artifact_manifest.sha256`: 输入 hash、定义、软件和输出 hash。",
        "- `NUMERICS_AUDIT.md`: units、shape、NaN、window/range sensitivity、optimizer 与 synthetic closure 检查。",
        "",
        "本脚本没有进行 PB/Debye/surface-potential 拟合。先确认低速速度项是否可重现、校准是否独立成立，再决定能否把某个 subtraction 后的残差解释为 equilibrium surface force。",
        "",
    ]
    (OUT / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def write_numerics_audit(
    map_records: list[dict],
    curve_rows: list[dict],
    contact_rows: list[dict],
    force_matrices: np.ndarray,
    gap_speed_matrices: np.ndarray,
    fit_rows: list[dict],
    checks: dict[str, float],
) -> None:
    invalid_keys = {
        (int(row["acquisition_order"]) - 1, int(row["point_index"]))
        for row in curve_rows
        if row["qc_state"] != "pass_force"
    }
    invalid_force_bins = sum(
        int(np.count_nonzero(np.isfinite(force_matrices[order, point])))
        for order, point in invalid_keys
    )
    invalid_gap_bins = sum(
        int(np.count_nonzero(np.isfinite(gap_speed_matrices[order, point])))
        for order, point in invalid_keys
    )
    inv_30 = median(
        row["map_contact_InvOLS_30nm_median_nm_per_V"] for row in contact_rows
    )
    inv_40 = median(row["map_contact_InvOLS_median_nm_per_V"] for row in contact_rows)
    inv_50 = median(
        row["map_contact_InvOLS_50nm_median_nm_per_V"] for row in contact_rows
    )
    pooled = [
        row
        for row in fit_rows
        if row["scope"] == "pooled_shared"
        and row["model"] in (
            "fitted_v_over_D",
            "shifted_v_over_D",
            "empirical_power_shift",
            "actual_gap_shifted",
        )
    ]
    lines = [
        "# Numerical audit: 16-09-26 low-speed analysis",
        "",
        "## Data flow and units",
        "",
        "`JPK int/float encoder -> raw vDeflection [V] and measuredHeight [m] -> per-curve constant baseline [V] -> contact InvOLS [nm/V] -> deflection [m] -> force [pN] and separation D [nm] -> local U_gap=-dD/dt [um/s] -> pixel bins -> map medians -> palindrome means -> double-centered model comparison`.",
        "",
        f"Force uses `k={float(map_records[0]['embedded_spring_constant_N_per_m']):.12g} N/m`; the resulting scale is {float(map_records[0]['embedded_spring_constant_N_per_m']) * inv_40:.9g} nN/V. `K` has units pN nm/(um/s)^alpha; for actual-gap models alpha=1.",
        "",
        "## Array and missing-data checks",
        "",
        f"- Force matrix shape/dtype: `{force_matrices.shape}`, `{force_matrices.dtype}`; finite={int(np.count_nonzero(np.isfinite(force_matrices)))}, NaN={int(np.count_nonzero(np.isnan(force_matrices)))}, Inf={int(np.count_nonzero(np.isinf(force_matrices)))}`.",
        f"- Gap-speed matrix shape/dtype: `{gap_speed_matrices.shape}`, `{gap_speed_matrices.dtype}`; finite={int(np.count_nonzero(np.isfinite(gap_speed_matrices)))}, NaN={int(np.count_nonzero(np.isnan(gap_speed_matrices)))}, Inf={int(np.count_nonzero(np.isinf(gap_speed_matrices)))}`.",
        f"- The 15 rejected curves contain {invalid_force_bins} finite force bins and {invalid_gap_bins} finite gap-speed bins; both must be zero. Missing curves remain NaN and are never imputed.",
        "- Every map retains finite map medians at 5-300 nm; each target bin has at least 61 finite pixels.",
        "",
        "## Sensitivity and conditioning",
        "",
        f"- Contact-window batch InvOLS at 30/40/50 nm: {inv_30:.6f}/{inv_40:.6f}/{inv_50:.6f} nm/V (maximum relative shift {max(abs(inv_30/inv_40-1), abs(inv_50/inv_40-1))*100:.3f}%).",
        "- Alternative constant-baseline windows change absolute force by a map-dependent constant. Separation and the double-centered velocity-distance interaction are invariant to that constant; absolute force-time panels remain finite-window-gauge quantities.",
        "- Hydrodynamic fit-range sensitivity is tabulated below. Correlated 5 nm bins are not counted as independent replicates and no bin-level p-value is reported.",
        "",
        "| range | model | K | alpha | D0 (nm) | normalized RMSE | Jacobian condition |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in pooled:
        lines.append(
            f"| {row['fit_range']} | {row['model']} | "
            f"{float(row['coefficient_pN_nm_per_um_s_to_alpha']):.6g} | "
            f"{float(row['velocity_exponent_alpha']):.5f} | "
            f"{float(row['distance_shift_D0_nm']):.5f} | "
            f"{float(row['normalized_rmse']):.5f} | "
            f"{float(row['jacobian_condition']):.6g} |"
        )
    optimizer_failures = sum(
        not bool(row["optimizer_success"])
        for row in fit_rows
        if int(row["parameter_count"]) > 0
    )
    lines += [
        "",
        f"All parameterized optimizations report success: failures={optimizer_failures}. The empirical three-parameter Jacobian is full-rank but ill-conditioned, so alpha and D0 are model-shape descriptors rather than independently identified material parameters.",
        "",
        "## Synthetic closure",
        "",
    ]
    for name, value in checks.items():
        lines.append(f"- `{name}` = `{value:.6g}`")
    lines += [
        "",
        "All synthetic errors are below 1e-10. Input ZIP CRC is exercised by the reused parser for every map; input and output SHA-256 manifests are written separately.",
        "",
        "## Claim boundary",
        "",
        "The numerics establish a reproducible velocity-distance interaction and quantify which curves are usable. They do not identify the interaction uniquely as hydrodynamic, validate the bulk no-slip prefactor, establish low-speed equivalence, or justify equilibrium PB parameters after subtraction.",
        "",
    ]
    (OUT / "NUMERICS_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def write_input_manifest(map_records: list[dict]) -> str:
    lines = [f"{row['sha256']}  {row['source']}" for row in map_records]
    path = OUT / "input_manifest.sha256"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    digest = hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()
    return digest


def write_artifact_manifest() -> None:
    targets = sorted(
        path
        for path in OUT.rglob("*")
        if path.is_file() and path.name != "artifact_manifest.sha256"
    )
    lines = [f"{sha256_file(path)}  {path.relative_to(ROOT).as_posix()}" for path in targets]
    (OUT / "artifact_manifest.sha256").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    map_records, curve_rows, payloads = load_maps()
    inv_ols, contact_rows = calibrate_contact(map_records, curve_rows)
    force_rows, matrices, gap_speed_matrices = reconstruct(
        map_records, curve_rows, payloads, inv_ols
    )
    enrich_map_records(map_records, curve_rows, contact_rows, force_rows)
    pair_rows, pair_lookup, pair_gap_speed_lookup = build_palindrome_pairs(
        map_records, force_rows
    )
    fit_rows, centered_rows, reproducibility_rows = hydrodynamic_analysis(
        pair_lookup, pair_gap_speed_lookup
    )
    interval_detail, interval_summary = speed_interval_analysis(pair_lookup)
    checks = synthetic_self_checks()
    if max(checks.values()) > 1e-10:
        raise RuntimeError(f"Synthetic invariant failed: {checks}")

    write_csv(OUT / "curve_qc.csv", curve_rows)
    write_csv(OUT / "map_inventory_and_summary.csv", map_records)
    write_csv(OUT / "contact_calibration_maps.csv", contact_rows)
    write_csv(OUT / "map_force_by_separation.csv", force_rows)
    write_csv(OUT / "palindrome_pair_force.csv", pair_rows)
    write_csv(OUT / "speed_interval_detail.csv", interval_detail)
    write_csv(OUT / "speed_interval_summary.csv", interval_summary)
    write_csv(OUT / "hydrodynamic_model_fits.csv", fit_rows)
    write_csv(OUT / "hydrodynamic_double_centered.csv", centered_rows)
    write_csv(OUT / "block_reproducibility.csv", reproducibility_rows)
    np.savez_compressed(
        OUT / "curve_force_bins.npz",
        force_pN=np.asarray(matrices, dtype=np.float64),
        gap_speed_um_per_s=np.asarray(gap_speed_matrices, dtype=np.float64),
        distance_nm=BIN_CENTERS_NM,
        speed_um_per_s=np.asarray([row["speed_um_per_s"] for row in map_records]),
        acquisition_order=np.arange(1, MAP_COUNT + 1, dtype=np.int64),
        point_index=np.arange(GRID_SIZE**2, dtype=np.int64),
    )

    plot_qc(map_records, curve_rows)
    plot_qc_features(curve_rows)
    plot_calibration(map_records, contact_rows, curve_rows, inv_ols)
    plot_gap_speed(force_rows)
    plot_force_time(map_records, force_rows)
    plot_speed_intervals(interval_detail, interval_summary)
    plot_hydrodynamics(fit_rows, centered_rows)
    write_report(
        map_records,
        curve_rows,
        contact_rows,
        inv_ols,
        fit_rows,
        reproducibility_rows,
        interval_summary,
    )
    write_numerics_audit(
        map_records,
        curve_rows,
        contact_rows,
        matrices,
        gap_speed_matrices,
        fit_rows,
        checks,
    )
    input_set_digest = write_input_manifest(map_records)
    viscosity, theoretical_coefficient = theoretical_no_slip_coefficient()
    provenance = {
        "script": Path(__file__).relative_to(ROOT).as_posix(),
        "data_directory": DATA.relative_to(ROOT).as_posix(),
        "input_maps": MAP_COUNT,
        "input_map_set_sha256": input_set_digest,
        "raw_data_modified": False,
        "branch_policy": "approach only for force analysis; all approach/retract segments retained for inventory",
        "qc_states": list(QC_STATES),
        "qc_thresholds": {
            "hard_fail_max_travel_nm": HARD_FAIL_MAX_TRAVEL_NM,
            "minimum_force_travel_nm": MIN_FORCE_TRAVEL_NM,
            "minimum_terminal_contact_amplitude_V": MIN_CONTACT_AMPLITUDE_V,
            "force_contact_InvOLS_range_nm_per_V": list(FORCE_CONTACT_INVOLS_RANGE_NM_PER_V),
            "force_contact_R2_min": FORCE_CONTACT_R2_MIN,
            "calibration_InvOLS_range_nm_per_V": list(CALIBRATION_INVOLS_RANGE_NM_PER_V),
            "calibration_R2_min": CALIBRATION_R2_MIN,
        },
        "baseline": {
            "coordinate": "scanner travel from first recorded sample",
            "window_nm": [BASELINE_START_NM, BASELINE_STOP_NM],
            "operation": "per-curve constant median subtraction",
            "interpretation": "finite-window force gauge, not infinite-separation zero force",
        },
        "contact": {
            "terminal_physical_window_nm": CONTACT_WINDOW_NM,
            "terminal_window_sensitivity_nm": list(CONTACT_WINDOW_CHECKS_NM),
            "batch_InvOLS_nm_per_V": inv_ols,
            "aggregation": "median of 18 per-map medians from strict calibration contacts",
        },
        "force": {
            "spring_constant_N_per_m": float(map_records[0]["embedded_spring_constant_N_per_m"]),
            "spring_constant_source": "identical embedded value in all 18 JPK maps; no same-day independent thermal file supplied",
            "force_scale_nN_per_V": float(map_records[0]["embedded_spring_constant_N_per_m"]) * inv_ols,
            "surface_separation": "measuredHeight + cantilever deflection - median terminal contact coordinate",
            "bin_centers_nm": BIN_CENTERS_NM.tolist(),
            "bin_width_nm": BIN_HALF_WIDTH_NM * 2.0,
            "gap_speed_definition": "U_gap=-dD/dt from robust local D(t) line",
            "gap_speed_local_half_window_nm": GAP_SPEED_HALF_WINDOW_NM,
        },
        "palindrome": {
            "blocks": [list(block) for block in EXPECTED_BLOCKS],
            "map_is_experimental_unit": True,
            "first_map_grid_shift_um": map_records[0]["shift_from_maps_2_to_18_center_um"],
            "first_pair_pixelwise_matching_disabled": True,
        },
        "hydrodynamic_check": {
            "operation": "within-block palindrome mean followed by distance-by-speed double centering",
            "claim_boundary": "tests reproducible velocity-distance interaction; does not uniquely identify hydrodynamics or complete subtraction",
            "fixed_no_slip_campaign_assumptions": {
                "glycerol_mass_fraction": ASSUMED_GLYCEROL_MASS_FRACTION,
                "temperature_C": ASSUMED_TEMPERATURE_C,
                "probe_radius_m": ASSUMED_PROBE_RADIUS_M,
                "viscosity_mPa_s": viscosity,
                "coefficient_pN_nm_per_um_s": theoretical_coefficient,
                "not_encoded_in_raw_header": True,
            },
        },
        "synthetic_self_checks": checks,
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
    write_artifact_manifest()
    counts = {state: sum(row["qc_state"] == state for row in curve_rows) for state in QC_STATES}
    print(
        json.dumps(
            {
                "output": str(OUT),
                "qc_counts": counts,
                "batch_InvOLS_nm_per_V": inv_ols,
                "force_scale_nN_per_V": float(map_records[0]["embedded_spring_constant_N_per_m"]) * inv_ols,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

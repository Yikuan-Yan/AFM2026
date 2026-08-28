#!/usr/bin/env python3
"""Paired-map velocity diagnostics and evidence-bounded joint AFM fitting.

This script starts from the raw 20-08-26 JPK archives.  It reuses the audited
raw channel decoder and hard-contact calibration, but it does not use JPK's
embedded force or sensitivity.  The primary calculations are:

1. same-pixel differences between sequential maps;
2. separation-resolved tests of the sphere-plane lubrication signature;
3. approach/retract odd-even diagnostics only on verified free retract data;
4. joint equal-potential silica sphere-plane PB fits with map-level optical
   scale and linear-drift nuisance parameters.

The map, not the pixel, is the velocity experimental unit.  Pixel pairing and
spatial block bootstrap quantify within-area heterogeneity; they do not create
independent velocity replication.  Ten wt% glycerol is retained as QC only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.optimize import least_squares
from scipy.signal import savgol_filter
from scipy.stats import spearmanr

import analyze_velocity_systematics as events
import fit_glycerol_surface_forces as base


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT / "20-08-26"
DEFAULT_RESULTS = ROOT / "analysis" / "velocity_joint_fit_results"

PRIMARY_CONCENTRATIONS = (0, 20, 30, 40)
QC_EXCLUDED_CONCENTRATION = 10
RANDOM_SEED = 20260821

BIN_CENTERS_NM = np.arange(20.0, 905.0, 5.0, dtype=np.float64)
BIN_HALF_WIDTH_NM = 2.5
JOINT_DMIN_NM = 25.0
JOINT_DMAX_NM = 250.0
PAIR_DMIN_NM = 25.0
PAIR_DMAX_NM = 250.0
SELECTED_MAP_DISTANCES_NM = (30.0, 50.0, 100.0, 200.0)
SPATIAL_BLOCK_SIZE = 4
BOOTSTRAP_REPLICATES = 1000
MIN_PAIR_PIXELS = 12
MIN_MAP_BIN_PIXELS = 12
MIN_BRANCH_SUPPORT_FRACTION = 0.50


@dataclass
class PointCurve:
    source: str
    concentration_wt_percent: int
    timestamp: str
    point_index: int
    acquisition_order: int
    row: int
    column: int
    sensitivity_m_per_V: float
    approach_speed_um_per_s: float
    approach_contact_response_ratio: float
    approach_contact_invOLS_nm_per_V: float
    approach_contact_height_um: float
    approach_terminal_load_nN: float
    approach_far_slope_pN_per_100nm: float
    approach_snap_detected: bool
    approach_snap_distance_nm: float
    approach_barrier_force_nN: float
    retract_free_baseline_valid: bool
    retract_snapoff_detected: bool
    retract_pull_off_censored: bool
    retract_detachment_right_censored: bool
    retract_detachment_piezo_travel_nm: float
    retract_pull_off_force_nN: float
    force_raw_reference_pN: np.ndarray
    force_linear_corrected_pN: np.ndarray
    gap_closing_speed_um_per_s: np.ndarray
    retract_free_force_corrected_pN: np.ndarray
    retract_opening_speed_um_per_s: np.ndarray


@dataclass
class MapData:
    source: str
    concentration_wt_percent: int
    timestamp: str
    speed_um_per_s: float
    speed_label_um_per_s: float
    load_regime: str
    primary_included: bool
    sensitivity_nm_per_V: float
    contact_response_center: float
    contact_response_mad: float
    contact_invOLS_center_nm_per_V: float
    terminal_load_center_nN: float
    far_slope_center_pN_per_100nm: float
    points: dict[tuple[int, int], PointCurve]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite(values: Iterable[float]) -> np.ndarray:
    array = np.asarray(list(values), dtype=np.float64)
    return array[np.isfinite(array)]


def safe_median(values: Iterable[float]) -> float:
    array = finite(values)
    return float(np.median(array)) if array.size else float("nan")


def speed_label(speed_um_per_s: float) -> float:
    labels = base.SPEED_LABELS_UM_PER_S
    return float(labels[np.argmin(np.abs(labels - float(speed_um_per_s)))])


def hydrodynamic_force_pN(
    distance_nm: np.ndarray | float,
    speed_um_per_s: np.ndarray | float,
    viscosity_mPa_s: float,
) -> np.ndarray:
    """Leading no-slip sphere-plane lubrication force magnitude.

    ``F = 6*pi*eta*R^2*U/D``.  Positive ``U`` denotes approach/closing and
    returns a positive resisting force.  Inputs are converted to SI before the
    returned force is converted to pN.
    """

    distance_m = np.asarray(distance_nm, dtype=np.float64) * 1e-9
    speed_m_per_s = np.asarray(speed_um_per_s, dtype=np.float64) * 1e-6
    viscosity_Pa_s = float(viscosity_mPa_s) * 1e-3
    if (
        np.any(~np.isfinite(distance_m))
        or np.any(distance_m <= 0.0)
        or np.any(~np.isfinite(speed_m_per_s))
        or not np.isfinite(viscosity_Pa_s)
        or viscosity_Pa_s <= 0.0
    ):
        raise ValueError("hydrodynamic inputs must be finite with D and eta positive")
    force_N = (
        6.0
        * np.pi
        * viscosity_Pa_s
        * base.PROBE_RADIUS_M**2
        * speed_m_per_s
        / distance_m
    )
    return force_N * 1e12


def bin_median(
    distance_nm: np.ndarray,
    values: np.ndarray,
    usable: np.ndarray,
) -> np.ndarray:
    distance = np.asarray(distance_nm, dtype=np.float64)
    data = np.asarray(values, dtype=np.float64)
    mask = np.asarray(usable, dtype=bool)
    if distance.shape != data.shape or distance.shape != mask.shape:
        raise ValueError("binning arrays must have identical shapes")
    output = np.full(BIN_CENTERS_NM.shape, np.nan, dtype=np.float64)
    valid = (
        mask
        & np.isfinite(distance)
        & np.isfinite(data)
        & (distance >= BIN_CENTERS_NM[0] - BIN_HALF_WIDTH_NM)
        & (distance < BIN_CENTERS_NM[-1] + BIN_HALF_WIDTH_NM)
    )
    if not np.any(valid):
        return output
    indices = np.floor(
        (distance[valid] - (BIN_CENTERS_NM[0] - BIN_HALF_WIDTH_NM)) / 5.0
    ).astype(np.int64)
    selected = data[valid]
    for index in np.unique(indices):
        if 0 <= index < output.size:
            output[index] = float(np.median(selected[indices == index]))
    return output


def smoothed_velocity_um_per_s(
    coordinate_nm: np.ndarray,
    duration_s: float,
    closing: bool,
) -> np.ndarray:
    values = np.asarray(coordinate_nm, dtype=np.float64)
    if values.ndim != 1 or values.size < 11 or duration_s <= 0.0:
        raise ValueError("velocity calculation needs a sufficiently long trace")
    delta_t = float(duration_s) / values.size
    maximum = min(51, values.size if values.size % 2 else values.size - 1)
    window = max(11, maximum)
    if window % 2 == 0:
        window -= 1
    derivative_nm_per_s = savgol_filter(
        values,
        window_length=window,
        polyorder=2,
        deriv=1,
        delta=delta_t,
        mode="interp",
    )
    velocity = (-1.0 if closing else 1.0) * derivative_nm_per_s / 1e3
    return np.asarray(velocity, dtype=np.float64)


def _curve_key(curve: base.RawCurve | events.BranchCurve, order: int) -> int:
    return int(curve.point_index) if curve.point_index is not None else int(order)


def prepare_point_curve(
    source: base.SourceData,
    approach: base.RawCurve,
    raw_approach: events.BranchCurve | None,
    retract: events.BranchCurve | None,
    order: int,
    event_record: dict,
) -> PointCurve | None:
    if (
        source.source_type != "map"
        or approach.point_index is None
        or source.map_grid_i is None
        or approach.contact_fit is None
        or approach.far_field_fit is None
    ):
        return None
    sensitivity = float(source.sensitivity_used_m_per_V)
    try:
        app_far, app_delta, app_force_corrected = events.corrected_branch(
            approach, sensitivity, approach.far_field_fit
        )
        app_corrected_V = approach.deflection_V - base.baseline_voltage(
            approach.measured_height_m, approach.deflection_V, app_far
        )
        app_contact_height = events.contact_height_m(
            approach.measured_height_m,
            app_corrected_V,
            approach.contact_fit,
        )
    except (ValueError, FloatingPointError, np.linalg.LinAlgError):
        return None

    distance_nm = (
        approach.measured_height_m + app_delta - app_contact_height
    ) * 1e9
    try:
        closing_speed = smoothed_velocity_um_per_s(
            distance_nm, approach.duration_s, closing=True
        )
    except (ValueError, FloatingPointError):
        closing_speed = np.full(distance_nm.shape, np.nan)
    scanner_speed = events.estimate_speed_um_per_s(
        approach.measured_height_m, approach.duration_s
    )
    plausible_speed = (
        np.isfinite(closing_speed)
        & (closing_speed >= 0.20 * scanner_speed)
        & (closing_speed <= 3.0 * scanner_speed)
    )
    closing_speed = np.where(plausible_speed, closing_speed, scanner_speed)

    far_unsaturated = np.ones(approach.measured_height_m.shape, dtype=bool)
    if raw_approach is not None:
        far_unsaturated = ~raw_approach.raw_saturation_mask
    far_limit = min(app_far.n_points, approach.deflection_V.size)
    far_reference_values = approach.deflection_V[:far_limit][
        far_unsaturated[:far_limit]
    ]
    if far_reference_values.size < 20:
        return None
    voltage_reference = float(np.median(far_reference_values))
    force_raw_reference = (
        base.SPRING_CONSTANT_N_PER_M
        * sensitivity
        * (approach.deflection_V - voltage_reference)
        * 1e12
    )
    precontact = np.arange(distance_nm.size) < approach.contact_fit.start
    usable = precontact & far_unsaturated
    force_raw_binned = bin_median(distance_nm, force_raw_reference, usable)
    force_corrected_binned = bin_median(distance_nm, app_force_corrected, usable)
    speed_binned = bin_median(distance_nm, closing_speed, usable)

    retract_force_binned = np.full(BIN_CENTERS_NM.shape, np.nan)
    retract_speed_binned = np.full(BIN_CENTERS_NM.shape, np.nan)
    if (
        retract is not None
        and bool(event_record.get("retract_free_baseline_valid", False))
        and bool(event_record.get("retract_snapoff_detected", False))
        and not bool(event_record.get("retract_detachment_right_censored", False))
        and np.isfinite(
            float(event_record.get("retract_detachment_free_side_distance_nm", np.nan))
        )
    ):
        try:
            ret_far_fit = events.fit_far_field_free_prefix(retract)
            _, ret_delta, ret_force = events.corrected_branch(
                retract, sensitivity, ret_far_fit
            )
            ret_distance_nm = (
                retract.measured_height_m + ret_delta - app_contact_height
            ) * 1e9
            opening_speed = smoothed_velocity_um_per_s(
                ret_distance_nm, retract.duration_s, closing=False
            )
            ret_scanner_speed = events.estimate_speed_um_per_s(
                retract.measured_height_m, retract.duration_s
            )
            ret_plausible = (
                np.isfinite(opening_speed)
                & (opening_speed >= 0.20 * ret_scanner_speed)
                & (opening_speed <= 3.0 * ret_scanner_speed)
            )
            opening_speed = np.where(
                ret_plausible, opening_speed, ret_scanner_speed
            )
            free_threshold = float(
                event_record["retract_detachment_free_side_distance_nm"]
            )
            ret_usable = (
                ~retract.raw_saturation_mask
                & (ret_distance_nm >= free_threshold)
            )
            retract_force_binned = bin_median(
                ret_distance_nm, ret_force, ret_usable
            )
            retract_speed_binned = bin_median(
                ret_distance_nm, opening_speed, ret_usable
            )
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            pass

    row, column = base.map_pixel_from_index(
        approach.point_index,
        source.map_grid_i,
        source.map_back_and_forth,
    )
    return PointCurve(
        source=str(source.path.relative_to(ROOT)),
        concentration_wt_percent=source.concentration_wt_percent,
        timestamp=source.timestamp.isoformat(),
        point_index=int(approach.point_index),
        acquisition_order=int(order),
        row=int(row),
        column=int(column),
        sensitivity_m_per_V=sensitivity,
        approach_speed_um_per_s=float(scanner_speed),
        approach_contact_response_ratio=float(
            event_record.get("approach_contact_response_ratio", np.nan)
        ),
        approach_contact_invOLS_nm_per_V=float(
            event_record.get("approach_contact_invOLS_nm_per_V", np.nan)
        ),
        approach_contact_height_um=float(
            event_record.get("approach_contact_height_um", np.nan)
        ),
        approach_terminal_load_nN=float(
            event_record.get("approach_terminal_load_nN", np.nan)
        ),
        approach_far_slope_pN_per_100nm=float(
            event_record.get("approach_far_slope_pN_per_100nm", np.nan)
        ),
        approach_snap_detected=bool(
            event_record.get("approach_snap_detected", False)
        ),
        approach_snap_distance_nm=float(
            event_record.get("approach_snap_distance_nm", np.nan)
        ),
        approach_barrier_force_nN=float(
            event_record.get("approach_barrier_force_nN", np.nan)
        ),
        retract_free_baseline_valid=bool(
            event_record.get("retract_free_baseline_valid", False)
        ),
        retract_snapoff_detected=bool(
            event_record.get("retract_snapoff_detected", False)
        ),
        retract_pull_off_censored=bool(
            event_record.get("retract_pull_off_censored", False)
        ),
        retract_detachment_right_censored=bool(
            event_record.get("retract_detachment_right_censored", False)
        ),
        retract_detachment_piezo_travel_nm=float(
            event_record.get("retract_detachment_piezo_travel_nm", np.nan)
        ),
        retract_pull_off_force_nN=float(
            event_record.get("retract_pull_off_force_nN", np.nan)
        ),
        force_raw_reference_pN=force_raw_binned,
        force_linear_corrected_pN=force_corrected_binned,
        gap_closing_speed_um_per_s=speed_binned,
        retract_free_force_corrected_pN=retract_force_binned,
        retract_opening_speed_um_per_s=retract_speed_binned,
    )


def classify_load_regime(concentration: int, terminal_load_nN: float) -> str:
    if concentration == 30:
        return "high_load" if terminal_load_nN >= 12.0 else "low_load"
    return "main"


def load_raw_maps() -> tuple[list[base.SourceData], list[MapData], list[dict]]:
    sources: list[base.SourceData] = []
    for concentration in sorted(base.EPSILON_R):
        paths = sorted((DATA_ROOT / str(concentration)).glob("*.jpk-force"))
        paths += sorted((DATA_ROOT / str(concentration)).glob("*.jpk-force-map"))
        for path in paths:
            sources.append(base.load_source(path, concentration))
    sources.sort(key=lambda item: item.timestamp)
    base.calibrate_sensitivity(sources)

    map_data: list[MapData] = []
    event_rows: list[dict] = []
    for source in sources:
        if source.source_type != "map":
            continue
        raw_approaches, _ = events.load_branch(source.path, "extend")
        retracts, _ = events.load_branch(source.path, "retract")
        raw_by_key = {
            _curve_key(curve, order): curve
            for order, curve in enumerate(raw_approaches)
        }
        retract_by_key = {
            _curve_key(curve, order): curve
            for order, curve in enumerate(retracts)
        }
        point_curves: dict[tuple[int, int], PointCurve] = {}
        for order, approach in enumerate(source.curves):
            key = _curve_key(approach, order)
            retract = retract_by_key.get(key)
            try:
                record = events.analyze_pair(
                    approach, retract, source.sensitivity_used_m_per_V
                )
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                record = events._nan_record()
            point = prepare_point_curve(
                source,
                approach,
                raw_by_key.get(key),
                retract,
                order,
                record,
            )
            if point is None:
                continue
            point_curves[(point.row, point.column)] = point
            event_rows.append(
                {
                    "concentration_wt_percent": source.concentration_wt_percent,
                    "source": point.source,
                    "timestamp": point.timestamp,
                    "point_index": point.point_index,
                    "acquisition_order": point.acquisition_order,
                    "map_row": point.row,
                    "map_column": point.column,
                    "approach_speed_um_per_s": point.approach_speed_um_per_s,
                    "approach_contact_response_ratio": point.approach_contact_response_ratio,
                    "approach_contact_invOLS_nm_per_V": point.approach_contact_invOLS_nm_per_V,
                    "approach_contact_height_um": point.approach_contact_height_um,
                    "approach_terminal_load_nN": point.approach_terminal_load_nN,
                    "approach_far_slope_pN_per_100nm": point.approach_far_slope_pN_per_100nm,
                    "approach_snap_detected": point.approach_snap_detected,
                    "approach_snap_distance_nm": point.approach_snap_distance_nm
                    if point.approach_snap_detected
                    else np.nan,
                    "approach_barrier_force_nN": point.approach_barrier_force_nN
                    if point.approach_snap_detected
                    else np.nan,
                    "retract_free_baseline_valid": point.retract_free_baseline_valid,
                    "retract_snapoff_detected": point.retract_snapoff_detected,
                    "retract_pull_off_censored": point.retract_pull_off_censored,
                    "retract_detachment_right_censored": point.retract_detachment_right_censored,
                    "retract_detachment_piezo_travel_nm": point.retract_detachment_piezo_travel_nm
                    if point.retract_snapoff_detected
                    and not point.retract_detachment_right_censored
                    else np.nan,
                    "retract_pull_off_force_nN": point.retract_pull_off_force_nN
                    if not point.retract_pull_off_censored
                    else np.nan,
                }
            )
        if not point_curves:
            continue
        terminal = safe_median(
            point.approach_terminal_load_nN for point in point_curves.values()
        )
        map_data.append(
            MapData(
                source=str(source.path.relative_to(ROOT)),
                concentration_wt_percent=source.concentration_wt_percent,
                timestamp=source.timestamp.isoformat(),
                speed_um_per_s=safe_median(
                    point.approach_speed_um_per_s for point in point_curves.values()
                ),
                speed_label_um_per_s=speed_label(
                    safe_median(
                        point.approach_speed_um_per_s
                        for point in point_curves.values()
                    )
                ),
                load_regime=classify_load_regime(
                    source.concentration_wt_percent, terminal
                ),
                primary_included=source.concentration_wt_percent
                in PRIMARY_CONCENTRATIONS,
                sensitivity_nm_per_V=source.sensitivity_used_m_per_V * 1e9,
                contact_response_center=safe_median(
                    point.approach_contact_response_ratio
                    for point in point_curves.values()
                ),
                contact_response_mad=base.robust_mad(
                    finite(
                        point.approach_contact_response_ratio
                        for point in point_curves.values()
                    )
                ),
                contact_invOLS_center_nm_per_V=safe_median(
                    point.approach_contact_invOLS_nm_per_V
                    for point in point_curves.values()
                ),
                terminal_load_center_nN=terminal,
                far_slope_center_pN_per_100nm=safe_median(
                    point.approach_far_slope_pN_per_100nm
                    for point in point_curves.values()
                ),
                points=point_curves,
            )
        )
    map_data.sort(key=lambda item: item.timestamp)
    return sources, map_data, event_rows


PAIR_METRICS = {
    "contact_invOLS_nm_per_V": (
        "approach_contact_invOLS_nm_per_V",
        "nm/V",
    ),
    "contact_height_nm": ("approach_contact_height_um", "nm"),
    "terminal_load_nN": ("approach_terminal_load_nN", "nN"),
    "far_slope_pN_per_100nm": (
        "approach_far_slope_pN_per_100nm",
        "pN/100nm",
    ),
    "snap_distance_nm": ("approach_snap_distance_nm", "nm"),
    "barrier_force_nN": ("approach_barrier_force_nN", "nN"),
    "retract_detachment_piezo_travel_nm": (
        "retract_detachment_piezo_travel_nm",
        "nm",
    ),
    "retract_pull_off_force_nN": ("retract_pull_off_force_nN", "nN"),
}


def pair_id(first: MapData, second: MapData) -> str:
    return (
        f"c{first.concentration_wt_percent}_{first.load_regime}_"
        f"{first.timestamp[11:19].replace(':', '')}_{first.speed_label_um_per_s:g}_"
        f"to_{second.timestamp[11:19].replace(':', '')}_{second.speed_label_um_per_s:g}"
    )


def map_pairs(maps: list[MapData]) -> list[tuple[MapData, MapData]]:
    groups: dict[tuple[int, str], list[MapData]] = {}
    for item in maps:
        groups.setdefault(
            (item.concentration_wt_percent, item.load_regime), []
        ).append(item)
    output: list[tuple[MapData, MapData]] = []
    for members in groups.values():
        members.sort(key=lambda item: item.timestamp)
        output.extend(itertools.combinations(members, 2))
    return output


def metric_value(point: PointCurve, metric: str) -> float:
    attribute, _ = PAIR_METRICS[metric]
    value = float(getattr(point, attribute))
    if metric == "contact_height_nm":
        value *= 1e3
    if metric in ("snap_distance_nm", "barrier_force_nN"):
        return value if point.approach_snap_detected else float("nan")
    if metric == "retract_detachment_piezo_travel_nm":
        return (
            value
            if point.retract_snapoff_detected
            and not point.retract_detachment_right_censored
            else float("nan")
        )
    if metric == "retract_pull_off_force_nN":
        return value if not point.retract_pull_off_censored else float("nan")
    return value


def block_bootstrap_median_ci(
    rows: np.ndarray,
    columns: np.ndarray,
    values: np.ndarray,
    rng: np.random.Generator,
    replicates: int = BOOTSTRAP_REPLICATES,
) -> tuple[float, float, int]:
    row = np.asarray(rows, dtype=np.int64)
    column = np.asarray(columns, dtype=np.int64)
    data = np.asarray(values, dtype=np.float64)
    valid = np.isfinite(data)
    row = row[valid]
    column = column[valid]
    data = data[valid]
    if data.size < MIN_PAIR_PIXELS:
        return float("nan"), float("nan"), 0
    labels = (row // SPATIAL_BLOCK_SIZE) * 100 + column // SPATIAL_BLOCK_SIZE
    unique = np.unique(labels)
    blocks = [data[labels == label] for label in unique]
    if len(blocks) < 3:
        return float("nan"), float("nan"), len(blocks)
    boot = np.empty(replicates, dtype=np.float64)
    for index in range(replicates):
        selected = rng.integers(0, len(blocks), size=len(blocks))
        sample = np.concatenate([blocks[item] for item in selected])
        boot[index] = float(np.median(sample))
    low, high = np.quantile(boot, [0.025, 0.975])
    return float(low), float(high), len(blocks)


def robust_plane(
    rows: np.ndarray,
    columns: np.ndarray,
    values: np.ndarray,
) -> tuple[float, float, float, float, int]:
    row = np.asarray(rows, dtype=np.float64)
    column = np.asarray(columns, dtype=np.float64)
    data = np.asarray(values, dtype=np.float64)
    valid = np.isfinite(row + column + data)
    row = row[valid]
    column = column[valid]
    data = data[valid]
    if data.size < 12:
        return (float("nan"),) * 4 + (0,)
    design = np.column_stack(
        [
            np.ones(data.size),
            row - np.mean(row),
            column - np.mean(column),
        ]
    )
    keep = np.ones(data.size, dtype=bool)
    coefficients = np.full(3, np.nan)
    for _ in range(5):
        coefficients, _, rank, _ = np.linalg.lstsq(
            design[keep], data[keep], rcond=None
        )
        if rank != 3:
            return (float("nan"),) * 4 + (0,)
        residual = data - design @ coefficients
        sigma = base.robust_mad(residual[keep])
        if not np.isfinite(sigma) or sigma <= np.finfo(float).tiny:
            break
        new_keep = np.abs(residual - np.median(residual[keep])) <= 3.5 * sigma
        if np.count_nonzero(new_keep) < 12 or np.array_equal(new_keep, keep):
            break
        keep = new_keep
    prediction = design[keep] @ coefficients
    ss_res = float(np.sum((data[keep] - prediction) ** 2))
    ss_tot = float(np.sum((data[keep] - np.mean(data[keep])) ** 2))
    r2 = float("nan") if ss_tot <= 0.0 else 1.0 - ss_res / ss_tot
    return (
        float(coefficients[0]),
        float(coefficients[1]),
        float(coefficients[2]),
        r2,
        int(np.count_nonzero(keep)),
    )


def build_paired_metric_differences(
    maps: list[MapData],
    rng: np.random.Generator,
) -> tuple[list[dict], list[dict], list[dict]]:
    difference_rows: list[dict] = []
    summary_rows: list[dict] = []
    association_rows: list[dict] = []
    pairs = map_pairs(maps)
    for first, second in pairs:
        common = sorted(set(first.points) & set(second.points))
        identifier = pair_id(first, second)
        metric_differences: dict[str, dict[tuple[int, int], float]] = {}
        for metric, (_, unit) in PAIR_METRICS.items():
            current: dict[tuple[int, int], float] = {}
            for coordinate in common:
                first_value = metric_value(first.points[coordinate], metric)
                second_value = metric_value(second.points[coordinate], metric)
                difference = (
                    second_value - first_value
                    if np.isfinite(first_value) and np.isfinite(second_value)
                    else np.nan
                )
                current[coordinate] = difference
                if np.isfinite(difference):
                    difference_rows.append(
                        {
                            "pair_id": identifier,
                            "concentration_wt_percent": first.concentration_wt_percent,
                            "primary_included": first.primary_included,
                            "load_regime": first.load_regime,
                            "source_first": first.source,
                            "source_second": second.source,
                            "timestamp_first": first.timestamp,
                            "timestamp_second": second.timestamp,
                            "speed_first_um_per_s": first.speed_um_per_s,
                            "speed_second_um_per_s": second.speed_um_per_s,
                            "delta_speed_um_per_s": second.speed_um_per_s
                            - first.speed_um_per_s,
                            "same_speed_control": abs(
                                second.speed_um_per_s - first.speed_um_per_s
                            )
                            < 0.25,
                            "map_row": coordinate[0],
                            "map_column": coordinate[1],
                            "metric": metric,
                            "unit": unit,
                            "value_first": first_value,
                            "value_second": second_value,
                            "difference_second_minus_first": difference,
                        }
                    )
            metric_differences[metric] = current
            coordinates = [key for key, value in current.items() if np.isfinite(value)]
            values = np.asarray([current[key] for key in coordinates], dtype=float)
            row_values = np.asarray([key[0] for key in coordinates], dtype=int)
            column_values = np.asarray([key[1] for key in coordinates], dtype=int)
            if values.size:
                ci_low, ci_high, blocks = block_bootstrap_median_ci(
                    row_values, column_values, values, rng
                )
                _, row_gradient, column_gradient, plane_r2, plane_n = robust_plane(
                    row_values, column_values, values
                )
                acquisition = np.asarray(
                    [second.points[key].acquisition_order for key in coordinates],
                    dtype=float,
                )
                acquisition_correlation = (
                    float(np.corrcoef(acquisition, values)[0, 1])
                    if values.size >= 3
                    and np.std(acquisition) > 0.0
                    and np.std(values) > 0.0
                    else np.nan
                )
            else:
                ci_low = ci_high = np.nan
                blocks = 0
                row_gradient = column_gradient = plane_r2 = np.nan
                plane_n = 0
                acquisition_correlation = np.nan
            summary_rows.append(
                {
                    "pair_id": identifier,
                    "concentration_wt_percent": first.concentration_wt_percent,
                    "primary_included": first.primary_included,
                    "load_regime": first.load_regime,
                    "source_first": first.source,
                    "source_second": second.source,
                    "timestamp_first": first.timestamp,
                    "timestamp_second": second.timestamp,
                    "speed_first_um_per_s": first.speed_um_per_s,
                    "speed_second_um_per_s": second.speed_um_per_s,
                    "delta_speed_um_per_s": second.speed_um_per_s
                    - first.speed_um_per_s,
                    "same_speed_control": abs(
                        second.speed_um_per_s - first.speed_um_per_s
                    )
                    < 0.25,
                    "metric": metric,
                    "unit": unit,
                    "paired_pixel_count": int(values.size),
                    "median_difference_second_minus_first": safe_median(values),
                    "difference_mad": base.robust_mad(values),
                    "spatial_block_count": blocks,
                    "spatial_block_bootstrap_ci95_low": ci_low,
                    "spatial_block_bootstrap_ci95_high": ci_high,
                    "difference_row_gradient_per_pixel": row_gradient,
                    "difference_column_gradient_per_pixel": column_gradient,
                    "difference_plane_r2": plane_r2,
                    "difference_plane_inliers": plane_n,
                    "difference_correlation_with_acquisition_order": acquisition_correlation,
                    "inference_scope": "within_area_paired_not_independent_map_replication",
                }
            )

        snap = metric_differences["snap_distance_nm"]
        for predictor in (
            "contact_invOLS_nm_per_V",
            "contact_height_nm",
            "terminal_load_nN",
            "far_slope_pN_per_100nm",
        ):
            predictor_values = metric_differences[predictor]
            selected = [
                key
                for key in common
                if np.isfinite(snap.get(key, np.nan))
                and np.isfinite(predictor_values.get(key, np.nan))
            ]
            x = np.asarray([predictor_values[key] for key in selected], dtype=float)
            y = np.asarray([snap[key] for key in selected], dtype=float)
            if x.size >= 12 and np.std(x) > 0.0 and np.std(y) > 0.0:
                slope, intercept, _, inliers = base.robust_line(x, y)
                pearson = float(np.corrcoef(x, y)[0, 1])
                spear = float(spearmanr(x, y).statistic)
            else:
                slope = intercept = pearson = spear = np.nan
                inliers = 0
            association_rows.append(
                {
                    "pair_id": identifier,
                    "concentration_wt_percent": first.concentration_wt_percent,
                    "primary_included": first.primary_included,
                    "load_regime": first.load_regime,
                    "speed_first_um_per_s": first.speed_um_per_s,
                    "speed_second_um_per_s": second.speed_um_per_s,
                    "same_speed_control": abs(
                        second.speed_um_per_s - first.speed_um_per_s
                    )
                    < 0.25,
                    "response": "delta_snap_distance_nm",
                    "predictor": f"delta_{predictor}",
                    "paired_pixel_count": int(x.size),
                    "pearson_r": pearson,
                    "spearman_rho": spear,
                    "robust_slope": slope,
                    "robust_intercept": intercept,
                    "robust_inliers": inliers,
                }
            )
    return difference_rows, summary_rows, association_rows


def spatial_block_statistics(
    coordinates: list[tuple[int, int]],
    values: np.ndarray,
) -> tuple[float, float, float, int]:
    data = np.asarray(values, dtype=np.float64)
    if data.size != len(coordinates):
        raise ValueError("coordinate/value count mismatch")
    valid = np.isfinite(data)
    data = data[valid]
    selected_coordinates = [
        coordinate for coordinate, keep in zip(coordinates, valid) if keep
    ]
    if not data.size:
        return float("nan"), float("nan"), float("nan"), 0
    labels = np.asarray(
        [
            (row // SPATIAL_BLOCK_SIZE) * 100 + column // SPATIAL_BLOCK_SIZE
            for row, column in selected_coordinates
        ],
        dtype=np.int64,
    )
    block_medians = np.asarray(
        [float(np.median(data[labels == label])) for label in np.unique(labels)],
        dtype=np.float64,
    )
    block_mad = base.robust_mad(block_medians)
    block_se = (
        block_mad / math.sqrt(block_medians.size)
        if block_medians.size >= 2 and np.isfinite(block_mad)
        else float("nan")
    )
    return (
        float(np.median(data)),
        base.robust_mad(data),
        float(block_se),
        int(block_medians.size),
    )


def build_map_force_curves(maps: list[MapData]) -> list[dict]:
    rows: list[dict] = []
    for item in maps:
        coordinates = sorted(item.points)
        point_values = [item.points[key] for key in coordinates]
        for index, distance_nm in enumerate(BIN_CENTERS_NM):
            raw = np.asarray(
                [point.force_raw_reference_pN[index] for point in point_values],
                dtype=float,
            )
            corrected = np.asarray(
                [
                    point.force_linear_corrected_pN[index]
                    for point in point_values
                ],
                dtype=float,
            )
            gap_speed = np.asarray(
                [point.gap_closing_speed_um_per_s[index] for point in point_values],
                dtype=float,
            )
            raw_center, raw_mad, raw_block_se, raw_blocks = spatial_block_statistics(
                coordinates, raw
            )
            corrected_center, corrected_mad, _, _ = spatial_block_statistics(
                coordinates, corrected
            )
            speed_center, speed_mad, _, _ = spatial_block_statistics(
                coordinates, gap_speed
            )
            if not np.isfinite(raw_center):
                continue
            rows.append(
                {
                    "concentration_wt_percent": item.concentration_wt_percent,
                    "primary_included": item.primary_included,
                    "load_regime": item.load_regime,
                    "source": item.source,
                    "timestamp": item.timestamp,
                    "speed_um_per_s": item.speed_um_per_s,
                    "speed_label_um_per_s": item.speed_label_um_per_s,
                    "distance_nm": float(distance_nm),
                    "force_raw_far_constant_referenced_pN": raw_center,
                    "force_raw_spatial_mad_pN": raw_mad,
                    "force_raw_spatial_block_se_pN": raw_block_se,
                    "spatial_block_count": raw_blocks,
                    "force_linear_drift_corrected_pN": corrected_center,
                    "force_corrected_spatial_mad_pN": corrected_mad,
                    "gap_closing_speed_um_per_s": speed_center,
                    "gap_closing_speed_spatial_mad_um_per_s": speed_mad,
                    "valid_pixel_count": int(np.count_nonzero(np.isfinite(raw))),
                    "contact_response_center": item.contact_response_center,
                    "contact_response_mad": item.contact_response_mad,
                    "terminal_load_center_nN": item.terminal_load_center_nN,
                    "far_slope_center_pN_per_100nm": item.far_slope_center_pN_per_100nm,
                }
            )
    return rows


def build_paired_force_differences(
    maps: list[MapData],
    rng: np.random.Generator,
) -> list[dict]:
    rows: list[dict] = []
    for first, second in map_pairs(maps):
        common = sorted(set(first.points) & set(second.points))
        identifier = pair_id(first, second)
        viscosity = base.cheng_viscosity_mPa_s(
            first.concentration_wt_percent / 100.0,
            base.TEMPERATURE_C,
        )
        for index, distance_nm in enumerate(BIN_CENTERS_NM):
            coordinates: list[tuple[int, int]] = []
            raw_differences: list[float] = []
            corrected_differences: list[float] = []
            hydro_differences: list[float] = []
            gap_speed_differences: list[float] = []
            for coordinate in common:
                a = first.points[coordinate]
                b = second.points[coordinate]
                raw_a = a.force_raw_reference_pN[index]
                raw_b = b.force_raw_reference_pN[index]
                speed_a = a.gap_closing_speed_um_per_s[index]
                speed_b = b.gap_closing_speed_um_per_s[index]
                if not all(np.isfinite(value) for value in (raw_a, raw_b, speed_a, speed_b)):
                    continue
                coordinates.append(coordinate)
                raw_differences.append(float(raw_b - raw_a))
                corrected_a = a.force_linear_corrected_pN[index]
                corrected_b = b.force_linear_corrected_pN[index]
                corrected_differences.append(
                    float(corrected_b - corrected_a)
                    if np.isfinite(corrected_a) and np.isfinite(corrected_b)
                    else np.nan
                )
                hydro_a = float(
                    hydrodynamic_force_pN(distance_nm, speed_a, viscosity)
                )
                hydro_b = float(
                    hydrodynamic_force_pN(distance_nm, speed_b, viscosity)
                )
                hydro_differences.append(hydro_b - hydro_a)
                gap_speed_differences.append(speed_b - speed_a)
            raw_array = np.asarray(raw_differences, dtype=float)
            if raw_array.size < MIN_PAIR_PIXELS:
                continue
            row_array = np.asarray([key[0] for key in coordinates], dtype=int)
            column_array = np.asarray([key[1] for key in coordinates], dtype=int)
            ci_low, ci_high, block_count = block_bootstrap_median_ci(
                row_array,
                column_array,
                raw_array,
                rng,
                replicates=500,
            )
            _, _, block_se, _ = spatial_block_statistics(coordinates, raw_array)
            rows.append(
                {
                    "pair_id": identifier,
                    "concentration_wt_percent": first.concentration_wt_percent,
                    "primary_included": first.primary_included,
                    "load_regime": first.load_regime,
                    "source_first": first.source,
                    "source_second": second.source,
                    "timestamp_first": first.timestamp,
                    "timestamp_second": second.timestamp,
                    "speed_first_um_per_s": first.speed_um_per_s,
                    "speed_second_um_per_s": second.speed_um_per_s,
                    "delta_speed_um_per_s": second.speed_um_per_s
                    - first.speed_um_per_s,
                    "same_speed_control": abs(
                        second.speed_um_per_s - first.speed_um_per_s
                    )
                    < 0.25,
                    "distance_nm": float(distance_nm),
                    "paired_pixel_count": int(raw_array.size),
                    "spatial_block_count": block_count,
                    "median_raw_force_difference_pN": float(np.median(raw_array)),
                    "raw_force_difference_mad_pN": base.robust_mad(raw_array),
                    "raw_force_difference_block_se_pN": block_se,
                    "raw_force_difference_block_ci95_low_pN": ci_low,
                    "raw_force_difference_block_ci95_high_pN": ci_high,
                    "median_linear_corrected_force_difference_pN": safe_median(
                        corrected_differences
                    ),
                    "median_gap_speed_difference_um_per_s": safe_median(
                        gap_speed_differences
                    ),
                    "theoretical_no_slip_hydrodynamic_difference_pN": safe_median(
                        hydro_differences
                    ),
                    "viscosity_mPa_s": viscosity,
                    "inference_scope": "within_area_paired_not_independent_map_replication",
                }
            )
    return rows


def fit_pair_hydrodynamic_signature(rows: list[dict]) -> list[dict]:
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        grouped.setdefault(str(row["pair_id"]), []).append(row)
    output: list[dict] = []
    for identifier, records in grouped.items():
        selected = [
            row
            for row in records
            if PAIR_DMIN_NM <= float(row["distance_nm"]) <= PAIR_DMAX_NM
            and int(row["paired_pixel_count"]) >= MIN_PAIR_PIXELS
            and np.isfinite(float(row["median_raw_force_difference_pN"]))
            and np.isfinite(
                float(row["theoretical_no_slip_hydrodynamic_difference_pN"])
            )
        ]
        if len(selected) < 12:
            continue
        distance = np.asarray([float(row["distance_nm"]) for row in selected])
        force = np.asarray(
            [float(row["median_raw_force_difference_pN"]) for row in selected]
        )
        hydro = np.asarray(
            [
                float(row["theoretical_no_slip_hydrodynamic_difference_pN"])
                for row in selected
            ]
        )
        sigma = np.asarray(
            [
                max(
                    3.0,
                    float(row["raw_force_difference_block_se_pN"])
                    if np.isfinite(float(row["raw_force_difference_block_se_pN"]))
                    else 3.0,
                )
                for row in selected
            ]
        )
        centered = distance - np.mean(distance)
        null_design = np.column_stack([np.ones(distance.size), centered])
        weighted_null_design = null_design / sigma[:, None]
        weighted_force = force / sigma
        null_coefficients, _, null_rank, _ = np.linalg.lstsq(
            weighted_null_design, weighted_force, rcond=None
        )
        null_prediction = null_design @ null_coefficients
        null_residual = force - null_prediction
        null_weighted_rss = float(np.sum((null_residual / sigma) ** 2))
        same_speed = bool(selected[0]["same_speed_control"])
        hydro_scale = float(np.ptp(hydro))
        if same_speed or hydro_scale < 1e-6:
            chi = chi_se = chi_ci_low = chi_ci_high = np.nan
            prediction = null_prediction
            full_residual = null_residual
            full_rank = null_rank
            improvement = 0.0
            full_success = True
            full_weighted_rss = null_weighted_rss
            sign_agreement = np.nan
        else:
            full_design = np.column_stack(
                [np.ones(distance.size), centered, hydro]
            )
            weighted_full_design = full_design / sigma[:, None]
            coefficients, _, full_rank, _ = np.linalg.lstsq(
                weighted_full_design, weighted_force, rcond=None
            )
            chi = float(coefficients[2])
            prediction = full_design @ coefficients
            full_residual = force - prediction
            full_weighted_rss = float(np.sum((full_residual / sigma) ** 2))
            chi_se = np.nan
            if full_rank == 3 and distance.size > 3:
                try:
                    covariance = np.linalg.inv(
                        weighted_full_design.T @ weighted_full_design
                    )
                    covariance *= full_weighted_rss / (distance.size - 3)
                    chi_se = float(math.sqrt(max(0.0, covariance[2, 2])))
                except np.linalg.LinAlgError:
                    pass
            chi_ci_low = chi - 1.96 * chi_se if np.isfinite(chi_se) else np.nan
            chi_ci_high = chi + 1.96 * chi_se if np.isfinite(chi_se) else np.nan
            full_success = full_rank == 3
            improvement = (
                1.0 - full_weighted_rss / null_weighted_rss
                if null_weighted_rss > 0.0
                else np.nan
            )
            # Compare only the non-linear-in-D shape left after projecting both
            # the force and hydrodynamic basis away from the same weighted
            # intercept/linear-distance nuisance space.
            hydro_line, _, _, _ = np.linalg.lstsq(
                weighted_null_design, hydro / sigma, rcond=None
            )
            hydro_shape = hydro - null_design @ hydro_line
            shape_mask = np.abs(hydro_shape) > 0.25 * np.max(np.abs(hydro_shape))
            sign_agreement = (
                float(
                    np.mean(
                        np.sign(null_residual[shape_mask])
                        == np.sign(hydro_shape[shape_mask])
                    )
                )
                if np.any(shape_mask)
                else np.nan
            )
        force_scale = base.robust_mad(force)
        reference = selected[0]
        output.append(
            {
                "pair_id": identifier,
                "concentration_wt_percent": int(
                    reference["concentration_wt_percent"]
                ),
                "primary_included": bool(reference["primary_included"]),
                "load_regime": reference["load_regime"],
                "source_first": reference["source_first"],
                "source_second": reference["source_second"],
                "speed_first_um_per_s": float(reference["speed_first_um_per_s"]),
                "speed_second_um_per_s": float(reference["speed_second_um_per_s"]),
                "delta_speed_um_per_s": float(reference["delta_speed_um_per_s"]),
                "same_speed_control": same_speed,
                "distance_min_nm": float(np.min(distance)),
                "distance_max_nm": float(np.max(distance)),
                "distance_bin_count": int(distance.size),
                "chi_observed_over_no_slip": chi,
                "chi_local_se": chi_se,
                "chi_local_ci95_low": chi_ci_low,
                "chi_local_ci95_high": chi_ci_high,
                "hydrodynamic_term_weighted_rss_improvement_fraction": improvement,
                "detrended_shape_sign_agreement_with_hydrodynamic_difference": sign_agreement,
                "null_weighted_rss": null_weighted_rss,
                "full_weighted_rss": full_weighted_rss,
                "null_linear_drift_rmse_pN": float(
                    math.sqrt(np.mean(null_residual**2))
                ),
                "full_hydrodynamic_plus_drift_rmse_pN": float(
                    math.sqrt(np.mean(full_residual**2))
                ),
                "force_difference_robust_scale_pN": force_scale,
                "fit_rank": full_rank,
                "optimizer_success": full_success,
                "interpretation": "same_speed_history_control"
                if same_speed
                else "descriptive_shape_test_speed_confounded_with_map_history",
            }
        )
    return output


def build_branch_odd_even(maps: list[MapData]) -> tuple[list[dict], list[dict]]:
    rows: list[dict] = []
    summaries: list[dict] = []
    for item in maps:
        coordinates = sorted(item.points)
        viscosity = base.cheng_viscosity_mPa_s(
            item.concentration_wt_percent / 100.0,
            base.TEMPERATURE_C,
        )
        per_map: list[dict] = []
        for index, distance_nm in enumerate(BIN_CENTERS_NM):
            selected_coordinates: list[tuple[int, int]] = []
            even_values: list[float] = []
            odd_values: list[float] = []
            predicted_odd: list[float] = []
            approach_valid = 0
            for coordinate in coordinates:
                point = item.points[coordinate]
                app = point.force_linear_corrected_pN[index]
                ret = point.retract_free_force_corrected_pN[index]
                if np.isfinite(app):
                    approach_valid += 1
                if not np.isfinite(app) or not np.isfinite(ret):
                    continue
                app_speed = point.gap_closing_speed_um_per_s[index]
                ret_speed = point.retract_opening_speed_um_per_s[index]
                if not np.isfinite(app_speed) or not np.isfinite(ret_speed):
                    continue
                selected_coordinates.append(coordinate)
                even_values.append(0.5 * (app + ret))
                odd_values.append(0.5 * (app - ret))
                predicted_odd.append(
                    0.5
                    * float(
                        hydrodynamic_force_pN(distance_nm, app_speed, viscosity)
                        + hydrodynamic_force_pN(distance_nm, ret_speed, viscosity)
                    )
                )
            support = (
                len(selected_coordinates) / approach_valid
                if approach_valid > 0
                else 0.0
            )
            if not selected_coordinates:
                continue
            even_center, even_mad, _, blocks = spatial_block_statistics(
                selected_coordinates, np.asarray(even_values)
            )
            odd_center, odd_mad, _, _ = spatial_block_statistics(
                selected_coordinates, np.asarray(odd_values)
            )
            record = {
                "concentration_wt_percent": item.concentration_wt_percent,
                "primary_included": item.primary_included,
                "load_regime": item.load_regime,
                "source": item.source,
                "timestamp": item.timestamp,
                "speed_um_per_s": item.speed_um_per_s,
                "distance_nm": float(distance_nm),
                "approach_valid_pixel_count": approach_valid,
                "paired_free_branch_pixel_count": len(selected_coordinates),
                "paired_free_branch_support_fraction": support,
                "spatial_block_count": blocks,
                "even_conservative_candidate_pN": even_center,
                "even_spatial_mad_pN": even_mad,
                "odd_dissipative_candidate_pN": odd_center,
                "odd_spatial_mad_pN": odd_mad,
                "theoretical_no_slip_odd_pN": safe_median(predicted_odd),
                "qualified_50pct_support": support
                >= MIN_BRANCH_SUPPORT_FRACTION,
                "branch_processing": "per_branch_far_linear_corrected_free_retract_only",
            }
            rows.append(record)
            per_map.append(record)
        qualified = [
            row for row in per_map if bool(row["qualified_50pct_support"])
        ]
        joint_qualified = [
            row
            for row in qualified
            if JOINT_DMIN_NM <= float(row["distance_nm"]) <= JOINT_DMAX_NM
        ]
        summaries.append(
            {
                "concentration_wt_percent": item.concentration_wt_percent,
                "primary_included": item.primary_included,
                "load_regime": item.load_regime,
                "source": item.source,
                "timestamp": item.timestamp,
                "speed_um_per_s": item.speed_um_per_s,
                "minimum_distance_with_50pct_free_branch_support_nm": min(
                    (float(row["distance_nm"]) for row in qualified),
                    default=np.nan,
                ),
                "qualified_bins_25_250nm": len(joint_qualified),
                "odd_even_can_constrain_primary_fit_window": len(joint_qualified)
                >= 12,
                "interpretation": "free_retract_overlap_only_not_adhered_branch",
            }
        )
    return rows, summaries


def primary_joint_groups(maps: list[MapData]) -> list[tuple[int, str]]:
    groups = {
        (item.concentration_wt_percent, item.load_regime)
        for item in maps
        if item.primary_included
        and not (
            item.concentration_wt_percent == 30
            and item.load_regime == "low_load"
        )
    }
    return sorted(groups)


def fit_joint_model(
    map_curve_rows: list[dict],
    concentration: int,
    load_regime: str,
    velocity_model: str,
    dmin_nm: float,
    dmax_nm: float,
    scale_mode: str,
    variant: str,
    contact_shift_nm: float = 0.0,
    excluded_speed_label_um_per_s: float | None = None,
) -> tuple[dict, list[dict], list[dict]]:
    if velocity_model not in ("M0_no_hydrodynamic", "M1_no_slip_chi1", "M2_fitted_chi"):
        raise ValueError("unknown velocity model")
    if scale_mode not in ("contact_constrained", "fixed_unity"):
        raise ValueError("unknown scale mode")
    selected = [
        row
        for row in map_curve_rows
        if int(row["concentration_wt_percent"]) == concentration
        and str(row["load_regime"]) == load_regime
        and dmin_nm
        <= float(row["distance_nm"]) + float(contact_shift_nm)
        <= dmax_nm
        and (
            excluded_speed_label_um_per_s is None
            or abs(
                float(row["speed_label_um_per_s"])
                - float(excluded_speed_label_um_per_s)
            )
            >= 0.25
        )
        and int(row["valid_pixel_count"]) >= MIN_MAP_BIN_PIXELS
        and np.isfinite(float(row["force_raw_far_constant_referenced_pN"]))
        and np.isfinite(float(row["gap_closing_speed_um_per_s"]))
    ]
    sources = sorted({str(row["source"]) for row in selected})
    speed_labels = sorted(
        {float(row["speed_label_um_per_s"]) for row in selected}
    )
    if len(sources) < 2 or len(selected) < 24:
        return (
            {
                "concentration_wt_percent": concentration,
                "load_regime": load_regime,
                "variant": variant,
                "velocity_model": velocity_model,
                "scale_mode": scale_mode,
                "fit_valid": False,
                "invalid_reasons": "insufficient_map_data",
                "map_count": len(sources),
                "speed_group_count": len(speed_labels),
                "n_data": len(selected),
                "dmin_nm": dmin_nm,
                "dmax_nm": dmax_nm,
                "contact_shift_nm": contact_shift_nm,
                "excluded_speed_label_um_per_s": np.nan
                if excluded_speed_label_um_per_s is None
                else excluded_speed_label_um_per_s,
            },
            [],
            [],
        )

    source_index = {source: index for index, source in enumerate(sources)}
    observed_distance = np.asarray(
        [float(row["distance_nm"]) for row in selected]
    )
    distance = observed_distance + float(contact_shift_nm)
    force = np.asarray(
        [float(row["force_raw_far_constant_referenced_pN"]) for row in selected]
    )
    gap_speed = np.asarray(
        [float(row["gap_closing_speed_um_per_s"]) for row in selected]
    )
    map_index = np.asarray(
        [source_index[str(row["source"])] for row in selected], dtype=np.int64
    )
    sigma = np.asarray(
        [
            max(
                5.0,
                float(row["force_raw_spatial_block_se_pN"])
                if np.isfinite(float(row["force_raw_spatial_block_se_pN"]))
                else 5.0,
            )
            for row in selected
        ]
    )
    eps_r = base.EPSILON_R[concentration]
    viscosity = base.cheng_viscosity_mPa_s(
        concentration / 100.0, base.TEMPERATURE_C
    )
    hydro = hydrodynamic_force_pN(distance, gap_speed, viscosity)
    distance_center = 0.5 * (dmin_nm + dmax_nm)
    centered_distance = distance - distance_center

    response_centers = np.asarray(
        [
            safe_median(
                float(row["contact_response_center"])
                for row in selected
                if str(row["source"]) == source
            )
            for source in sources
        ]
    )
    response_mads = np.asarray(
        [
            safe_median(
                float(row["contact_response_mad"])
                for row in selected
                if str(row["source"]) == source
            )
            for source in sources
        ]
    )
    response_centers = np.where(
        np.isfinite(response_centers) & (response_centers > 0.0),
        response_centers,
        1.0,
    )
    response_sigma = np.maximum(
        np.where(np.isfinite(response_mads), response_mads, 0.0), 0.02
    )

    index = 0
    lambda_index = index
    index += 1
    zeta_index = index
    index += 1
    chi_index: int | None = None
    if velocity_model == "M2_fitted_chi":
        chi_index = index
        index += 1
    intercept_slice = slice(index, index + len(sources))
    index += len(sources)
    slope_slice = slice(index, index + len(sources))
    index += len(sources)
    scale_slice: slice | None = None
    if scale_mode == "contact_constrained":
        scale_slice = slice(index, index + len(sources))
        index += len(sources)
    parameter_count = index

    lower = np.full(parameter_count, -np.inf)
    upper = np.full(parameter_count, np.inf)
    lower[lambda_index] = math.log(1.0)
    upper[lambda_index] = math.log(1000.0)
    lower[zeta_index] = math.log(0.1)
    upper[zeta_index] = math.log(250.0)
    if chi_index is not None:
        lower[chi_index] = -20.0
        upper[chi_index] = 20.0
    lower[intercept_slice] = -5000.0
    upper[intercept_slice] = 5000.0
    lower[slope_slice] = -20.0
    upper[slope_slice] = 20.0
    if scale_slice is not None:
        lower[scale_slice] = math.log(0.75)
        upper[scale_slice] = math.log(1.25)

    def unpack(parameters: np.ndarray) -> tuple[float, float, float, np.ndarray, np.ndarray, np.ndarray]:
        lambda_nm = math.exp(float(parameters[lambda_index]))
        zeta_mV = math.exp(float(parameters[zeta_index]))
        if velocity_model == "M0_no_hydrodynamic":
            chi = 0.0
        elif velocity_model == "M1_no_slip_chi1":
            chi = 1.0
        else:
            assert chi_index is not None
            chi = float(parameters[chi_index])
        intercepts = np.asarray(parameters[intercept_slice], dtype=float)
        slopes = np.asarray(parameters[slope_slice], dtype=float)
        scales = (
            np.exp(np.asarray(parameters[scale_slice], dtype=float))
            if scale_slice is not None
            else np.ones(len(sources), dtype=float)
        )
        return lambda_nm, zeta_mV, chi, intercepts, slopes, scales

    def components(parameters: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        lambda_nm, zeta_mV, chi, intercepts, slopes, scales = unpack(parameters)
        equilibrium = base.total_equilibrium_force_pN(
            distance,
            lambda_nm,
            zeta_mV,
            0.0,
            eps_r,
            "nonlinear_pb_derjaguin",
        )
        hydrodynamic = chi * hydro
        physical = equilibrium + hydrodynamic
        nuisance = intercepts[map_index] + slopes[map_index] * centered_distance
        prediction = scales[map_index] * physical + nuisance
        return prediction, equilibrium, hydrodynamic, nuisance

    def residual(parameters: np.ndarray) -> np.ndarray:
        prediction, _, _, _ = components(parameters)
        result = [(prediction - force) / sigma]
        if scale_slice is not None:
            _, _, _, _, _, scales = unpack(parameters)
            result.append((scales - response_centers) / response_sigma)
        return np.concatenate(result)

    starts: list[tuple[float, float, float]] = [
        (20.0, 50.0, 1.0),
        (100.0, 100.0, 0.0),
    ]
    if velocity_model == "M2_fitted_chi":
        starts.append((200.0, 50.0, -1.0))
    best: tuple[float, scipy.optimize.OptimizeResult] | None = None
    for initial_lambda, initial_zeta, initial_chi in starts:
        parameters = np.zeros(parameter_count, dtype=float)
        parameters[lambda_index] = math.log(initial_lambda)
        parameters[zeta_index] = math.log(initial_zeta)
        if chi_index is not None:
            parameters[chi_index] = initial_chi
        if scale_slice is not None:
            parameters[scale_slice] = np.log(
                np.clip(response_centers, 0.76, 1.24)
            )
        # Initialize the map-specific line after removing the trial physical
        # force.  This materially improves conditioning without changing the
        # optimized model.
        trial_prediction, _, _, _ = components(parameters)
        for map_number in range(len(sources)):
            mask = map_index == map_number
            trial_without_line = trial_prediction[mask]
            line_target = force[mask] - trial_without_line
            try:
                slope, intercept, _, _ = base.robust_line(
                    centered_distance[mask], line_target
                )
            except (ValueError, np.linalg.LinAlgError):
                slope, intercept = 0.0, float(np.median(line_target))
            parameters[intercept_slice.start + map_number] = np.clip(
                intercept, -4900.0, 4900.0
            )
            parameters[slope_slice.start + map_number] = np.clip(
                slope, -19.0, 19.0
            )
        result = least_squares(
            residual,
            x0=parameters,
            bounds=(lower, upper),
            method="trf",
            loss="linear",
            xtol=1e-8,
            ftol=1e-8,
            gtol=1e-8,
            max_nfev=250,
        )
        objective = float(np.sum(residual(result.x) ** 2))
        if best is None or objective < best[0]:
            best = (objective, result)
    assert best is not None
    objective, result = best
    prediction, equilibrium, hydrodynamic_component, nuisance = components(result.x)
    lambda_nm, zeta_mV, chi, intercepts, slopes, scales = unpack(result.x)
    raw_residual = force - prediction
    normalized_data_residual = raw_residual / sigma
    weighted_rss = float(np.sum(normalized_data_residual**2))
    rmse = float(math.sqrt(np.mean(raw_residual**2)))
    ss_tot = float(np.sum((force - np.mean(force)) ** 2))
    r2 = float("nan") if ss_tot <= 0.0 else 1.0 - float(np.sum(raw_residual**2)) / ss_tot
    n_data = distance.size
    likelihood_observation_count = n_data + (
        len(sources) if scale_slice is not None else 0
    )
    joint_weighted_rss = objective
    # The residuals are already normalized by their separately estimated
    # uncertainties.  For a Gaussian model with those scales held fixed,
    # -2 log L differs from chi-square only by constants common to the three
    # primary nested models.  Using n*log(RSS/n) here would silently estimate
    # one extra common variance scale and would not match the stated weighted
    # likelihood.
    aic = joint_weighted_rss + 2.0 * parameter_count
    aicc = (
        aic
        + 2.0
        * parameter_count
        * (parameter_count + 1)
        / (likelihood_observation_count - parameter_count - 1)
        if likelihood_observation_count > parameter_count + 1
        else float("inf")
    )

    jacobian = np.asarray(result.jac, dtype=float)
    column_norm = np.linalg.norm(jacobian, axis=0)
    rank = 0
    condition = float("inf")
    se = np.full(parameter_count, np.nan)
    if np.all(np.isfinite(column_norm)) and np.all(column_norm > 0.0):
        normalized_jacobian = jacobian / column_norm
        singular = np.linalg.svd(normalized_jacobian, compute_uv=False)
        tolerance = max(normalized_jacobian.shape) * np.finfo(float).eps * singular[0]
        rank = int(np.count_nonzero(singular > tolerance))
        if singular[-1] > 0.0:
            condition = float(singular[0] / singular[-1])
        if rank == parameter_count and jacobian.shape[0] > parameter_count:
            try:
                covariance = np.linalg.inv(jacobian.T @ jacobian)
                covariance *= objective / (jacobian.shape[0] - parameter_count)
                se = np.sqrt(np.maximum(0.0, np.diag(covariance)))
            except np.linalg.LinAlgError:
                pass
    boundary_fraction = np.minimum(
        (result.x - lower) / (upper - lower),
        (upper - result.x) / (upper - lower),
    )
    finite_boundary = boundary_fraction[np.isfinite(boundary_fraction)]
    reasons: list[str] = []
    if not result.success:
        reasons.append("optimizer")
    if rank != parameter_count:
        reasons.append("jacobian_rank")
    if not np.isfinite(condition) or condition > 1e8:
        reasons.append("jacobian_condition")
    if finite_boundary.size and np.any(finite_boundary < 1e-4):
        reasons.append("parameter_boundary")
    if not np.isfinite(r2) or r2 < 0.50:
        reasons.append("r2")
    lambda_se_log = se[lambda_index]
    zeta_se_log = se[zeta_index]
    chi_se = se[chi_index] if chi_index is not None else np.nan
    lambda_ci = (
        lambda_nm * math.exp(-1.96 * lambda_se_log),
        lambda_nm * math.exp(1.96 * lambda_se_log),
    ) if np.isfinite(lambda_se_log) and lambda_se_log < 5.0 else (np.nan, np.nan)
    zeta_ci = (
        zeta_mV * math.exp(-1.96 * zeta_se_log),
        zeta_mV * math.exp(1.96 * zeta_se_log),
    ) if np.isfinite(zeta_se_log) and zeta_se_log < 5.0 else (np.nan, np.nan)
    chi_ci = (
        chi - 1.96 * chi_se,
        chi + 1.96 * chi_se,
    ) if np.isfinite(chi_se) else (np.nan, np.nan)

    fit_row = {
        "concentration_wt_percent": concentration,
        "load_regime": load_regime,
        "primary_included": concentration in PRIMARY_CONCENTRATIONS,
        "variant": variant,
        "velocity_model": velocity_model,
        "scale_mode": scale_mode,
        "fit_valid": not reasons,
        "invalid_reasons": ";".join(reasons),
        "experimental_velocity_identification": "confounded_with_sequential_map_history",
        "map_count": len(sources),
        "speed_group_count": len(speed_labels),
        "n_data": int(n_data),
        "contact_scale_observation_count": len(sources)
        if scale_slice is not None
        else 0,
        "likelihood_observation_count": int(likelihood_observation_count),
        "parameter_count": parameter_count,
        "dmin_nm": dmin_nm,
        "dmax_nm": dmax_nm,
        "contact_shift_nm": contact_shift_nm,
        "excluded_speed_label_um_per_s": np.nan
        if excluded_speed_label_um_per_s is None
        else excluded_speed_label_um_per_s,
        "lambda_D_nm": lambda_nm,
        "lambda_local_ci95_low_nm": lambda_ci[0],
        "lambda_local_ci95_high_nm": lambda_ci[1],
        "zeta_magnitude_mV": zeta_mV,
        "zeta_signed_silica_mV": -zeta_mV,
        "zeta_local_ci95_low_mV": zeta_ci[0],
        "zeta_local_ci95_high_mV": zeta_ci[1],
        "chi_observed_over_no_slip": chi,
        "chi_local_se": chi_se,
        "chi_local_ci95_low": chi_ci[0],
        "chi_local_ci95_high": chi_ci[1],
        "viscosity_mPa_s": viscosity,
        "force_weighted_rss": weighted_rss,
        "joint_force_plus_contact_scale_weighted_rss": joint_weighted_rss,
        "aic": aic,
        "aicc": aicc,
        "r2": r2,
        "rmse_pN": rmse,
        "jacobian_rank": rank,
        "jacobian_condition": condition,
        "optimizer_success": bool(result.success),
        "optimizer_status": int(result.status),
        "optimizer_nfev": int(result.nfev),
        "optimizer_message": str(result.message),
        "local_ci_interpretation": "conditional_bins_correlated_not_map_replication_ci",
    }
    prediction_rows: list[dict] = []
    for observation, row in enumerate(selected):
        prediction_rows.append(
            {
                "concentration_wt_percent": concentration,
                "load_regime": load_regime,
                "variant": variant,
                "velocity_model": velocity_model,
                "scale_mode": scale_mode,
                "source": row["source"],
                "timestamp": row["timestamp"],
                "speed_um_per_s": row["speed_um_per_s"],
                "observed_distance_nm": observed_distance[observation],
                "model_distance_nm": distance[observation],
                "gap_closing_speed_um_per_s": gap_speed[observation],
                "observed_force_pN": force[observation],
                "sigma_pN": sigma[observation],
                "predicted_force_pN": prediction[observation],
                "equilibrium_component_pN": equilibrium[observation],
                "hydrodynamic_component_pN": hydrodynamic_component[observation],
                "map_linear_nuisance_pN": nuisance[observation],
                "map_force_scale": scales[map_index[observation]],
                "residual_pN": raw_residual[observation],
            }
        )
    nuisance_rows = [
        {
            "concentration_wt_percent": concentration,
            "load_regime": load_regime,
            "variant": variant,
            "velocity_model": velocity_model,
            "scale_mode": scale_mode,
            "contact_shift_nm": contact_shift_nm,
            "excluded_speed_label_um_per_s": np.nan
            if excluded_speed_label_um_per_s is None
            else excluded_speed_label_um_per_s,
            "source": source,
            "map_intercept_pN_at_window_center": intercepts[number],
            "map_linear_slope_pN_per_nm": slopes[number],
            "fitted_map_force_scale": scales[number],
            "contact_response_constraint_center": response_centers[number],
            "contact_response_constraint_sigma": response_sigma[number],
        }
        for number, source in enumerate(sources)
    ]
    return fit_row, prediction_rows, nuisance_rows


def run_joint_models(
    maps: list[MapData],
    map_curve_rows: list[dict],
) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    fits: list[dict] = []
    predictions: list[dict] = []
    nuisance: list[dict] = []
    for concentration, regime in primary_joint_groups(maps):
        for model in (
            "M0_no_hydrodynamic",
            "M1_no_slip_chi1",
            "M2_fitted_chi",
        ):
            fit, prediction, current_nuisance = fit_joint_model(
                map_curve_rows,
                concentration,
                regime,
                model,
                JOINT_DMIN_NM,
                JOINT_DMAX_NM,
                "contact_constrained",
                "primary",
            )
            fits.append(fit)
            predictions.extend(prediction)
            nuisance.extend(current_nuisance)
        for dmin, dmax, scale_mode, variant in (
            (30.0, 250.0, "contact_constrained", "M2_dmin30"),
            (25.0, 150.0, "contact_constrained", "M2_dmax150"),
            (25.0, 350.0, "contact_constrained", "M2_dmax350"),
            (25.0, 250.0, "fixed_unity", "M2_fixed_unity_scale"),
        ):
            fit, prediction, current_nuisance = fit_joint_model(
                map_curve_rows,
                concentration,
                regime,
                "M2_fitted_chi",
                dmin,
                dmax,
                scale_mode,
                variant,
            )
            fits.append(fit)
            # Only primary predictions are saved at full resolution; variants
            # are summarized in the fit table to keep artifacts tractable.
            nuisance.extend(current_nuisance)
        for contact_shift_nm, variant in (
            (-2.5, "M2_contact_minus2p5nm"),
            (+2.5, "M2_contact_plus2p5nm"),
        ):
            fit, _, current_nuisance = fit_joint_model(
                map_curve_rows,
                concentration,
                regime,
                "M2_fitted_chi",
                JOINT_DMIN_NM,
                JOINT_DMAX_NM,
                "contact_constrained",
                variant,
                contact_shift_nm=contact_shift_nm,
            )
            fits.append(fit)
            nuisance.extend(current_nuisance)
        for excluded_speed in (1.0, 2.0, 4.0):
            fit, _, current_nuisance = fit_joint_model(
                map_curve_rows,
                concentration,
                regime,
                "M2_fitted_chi",
                JOINT_DMIN_NM,
                JOINT_DMAX_NM,
                "contact_constrained",
                f"M2_leave_out_{excluded_speed:g}um_s",
                excluded_speed_label_um_per_s=excluded_speed,
            )
            fits.append(fit)
            nuisance.extend(current_nuisance)

    for concentration, regime in primary_joint_groups(maps):
        primary = [
            row
            for row in fits
            if int(row["concentration_wt_percent"]) == concentration
            and str(row["load_regime"]) == regime
            and str(row["variant"]) == "primary"
            and np.isfinite(float(row.get("aicc", np.nan)))
        ]
        if primary:
            minimum = min(float(row["aicc"]) for row in primary)
            for row in primary:
                row["delta_aicc_within_primary_models"] = float(row["aicc"]) - minimum
                row["aicc_preferred_within_tested_models"] = (
                    abs(float(row["aicc"]) - minimum) < 1e-9
                )

    systematic_rows: list[dict] = []
    for concentration, regime in primary_joint_groups(maps):
        candidates = [
            row
            for row in fits
            if int(row["concentration_wt_percent"]) == concentration
            and str(row["load_regime"]) == regime
            and bool(row.get("fit_valid", False))
            and np.isfinite(float(row.get("lambda_D_nm", np.nan)))
            and np.isfinite(float(row.get("zeta_magnitude_mV", np.nan)))
        ]
        m2 = [
            row for row in candidates if row["velocity_model"] == "M2_fitted_chi"
        ]
        primary = [row for row in candidates if row["variant"] == "primary"]
        preferred = min(
            primary,
            key=lambda row: float(row.get("aicc", np.inf)),
            default=None,
        )
        systematic_rows.append(
            {
                "concentration_wt_percent": concentration,
                "load_regime": regime,
                "preferred_primary_model": ""
                if preferred is None
                else preferred["velocity_model"],
                "preferred_primary_delta_aicc": np.nan
                if preferred is None
                else preferred.get("delta_aicc_within_primary_models", np.nan),
                "lambda_all_tested_min_nm": min(
                    (float(row["lambda_D_nm"]) for row in candidates),
                    default=np.nan,
                ),
                "lambda_all_tested_max_nm": max(
                    (float(row["lambda_D_nm"]) for row in candidates),
                    default=np.nan,
                ),
                "zeta_all_tested_min_mV": min(
                    (float(row["zeta_magnitude_mV"]) for row in candidates),
                    default=np.nan,
                ),
                "zeta_all_tested_max_mV": max(
                    (float(row["zeta_magnitude_mV"]) for row in candidates),
                    default=np.nan,
                ),
                "M2_lambda_window_scale_min_nm": min(
                    (float(row["lambda_D_nm"]) for row in m2), default=np.nan
                ),
                "M2_lambda_window_scale_max_nm": max(
                    (float(row["lambda_D_nm"]) for row in m2), default=np.nan
                ),
                "M2_zeta_window_scale_min_mV": min(
                    (float(row["zeta_magnitude_mV"]) for row in m2),
                    default=np.nan,
                ),
                "M2_zeta_window_scale_max_mV": max(
                    (float(row["zeta_magnitude_mV"]) for row in m2),
                    default=np.nan,
                ),
                "M2_chi_min": min(
                    (float(row["chi_observed_over_no_slip"]) for row in m2),
                    default=np.nan,
                ),
                "M2_chi_max": max(
                    (float(row["chi_observed_over_no_slip"]) for row in m2),
                    default=np.nan,
                ),
                "claim_boundary": "exploratory_systematic_range_not_velocity_identified_parameter",
            }
        )
    return fits, predictions, nuisance, systematic_rows


def _symmetric_limit(values: Iterable[float], quantile: float = 0.98) -> float:
    array = np.abs(finite(values))
    if not array.size:
        return 1.0
    limit = float(np.quantile(array, quantile))
    return max(limit, np.finfo(float).eps)


def _pair_title(record: dict) -> str:
    concentration = int(record["concentration_wt_percent"])
    first = str(record["timestamp_first"])[11:19]
    second = str(record["timestamp_second"])[11:19]
    speed_first = float(record["speed_first_um_per_s"])
    speed_second = float(record["speed_second_um_per_s"])
    regime = str(record["load_regime"])
    return (
        f"{concentration}% {regime}\n"
        f"{first} {speed_first:.0f} -> {second} {speed_second:.0f} um/s"
    )


def plot_paired_metric_maps(
    rows: list[dict], figures_dir: Path
) -> list[str]:
    files: list[str] = []
    for metric in (
        "snap_distance_nm",
        "contact_invOLS_nm_per_V",
        "contact_height_nm",
        "terminal_load_nN",
        "far_slope_pN_per_100nm",
    ):
        selected = [
            row
            for row in rows
            if str(row["metric"]) == metric and bool(row["primary_included"])
        ]
        pair_ids = list(dict.fromkeys(str(row["pair_id"]) for row in selected))
        if not pair_ids:
            continue
        columns = 4
        panel_rows = math.ceil(len(pair_ids) / columns)
        figure, axes = plt.subplots(
            panel_rows,
            columns,
            figsize=(3.7 * columns, 3.35 * panel_rows),
            squeeze=False,
            constrained_layout=True,
        )
        limit = _symmetric_limit(
            float(row["difference_second_minus_first"]) for row in selected
        )
        image = None
        for axis, identifier in zip(axes.flat, pair_ids):
            records = [row for row in selected if str(row["pair_id"]) == identifier]
            matrix = np.full((16, 16), np.nan)
            for row in records:
                matrix[int(row["map_row"]), int(row["map_column"])] = float(
                    row["difference_second_minus_first"]
                )
            masked = np.ma.masked_invalid(matrix)
            image = axis.imshow(
                masked,
                origin="upper",
                cmap="coolwarm",
                vmin=-limit,
                vmax=limit,
                interpolation="nearest",
            )
            axis.set_title(_pair_title(records[0]), fontsize=8)
            axis.set_xticks((0, 5, 10, 15))
            axis.set_yticks((0, 5, 10, 15))
            axis.tick_params(labelsize=7)
        for axis in axes.flat[len(pair_ids) :]:
            axis.set_visible(False)
        if image is not None:
            unit = PAIR_METRICS[metric][1]
            figure.colorbar(
                image,
                ax=[axis for axis in axes.flat if axis.get_visible()],
                shrink=0.72,
                label=f"second - first ({unit})",
            )
        figure.suptitle(
            f"Same-pixel paired difference: {metric}", fontsize=13
        )
        filename = f"paired_{metric}_difference_maps.png"
        figure.savefig(figures_dir / filename, dpi=220)
        plt.close(figure)
        files.append(filename)
    return files


def plot_selected_force_difference_maps(
    maps: list[MapData], figures_dir: Path
) -> list[str]:
    files: list[str] = []
    primary_pairs = [
        pair for pair in map_pairs(maps) if pair[0].primary_included
    ]
    for target_distance in SELECTED_MAP_DISTANCES_NM:
        index = int(np.argmin(np.abs(BIN_CENTERS_NM - target_distance)))
        values: list[float] = []
        matrices: list[tuple[MapData, MapData, np.ndarray]] = []
        for first, second in primary_pairs:
            matrix = np.full((16, 16), np.nan)
            for coordinate in set(first.points) & set(second.points):
                a = first.points[coordinate].force_raw_reference_pN[index]
                b = second.points[coordinate].force_raw_reference_pN[index]
                if np.isfinite(a) and np.isfinite(b):
                    matrix[coordinate] = b - a
                    values.append(float(b - a))
            matrices.append((first, second, matrix))
        limit = _symmetric_limit(values)
        columns = 4
        panel_rows = math.ceil(len(matrices) / columns)
        figure, axes = plt.subplots(
            panel_rows,
            columns,
            figsize=(3.7 * columns, 3.35 * panel_rows),
            squeeze=False,
            constrained_layout=True,
        )
        image = None
        for axis, (first, second, matrix) in zip(axes.flat, matrices):
            image = axis.imshow(
                np.ma.masked_invalid(matrix),
                origin="upper",
                cmap="coolwarm",
                vmin=-limit,
                vmax=limit,
                interpolation="nearest",
            )
            axis.set_title(
                _pair_title(
                    {
                        "concentration_wt_percent": first.concentration_wt_percent,
                        "load_regime": first.load_regime,
                        "timestamp_first": first.timestamp,
                        "timestamp_second": second.timestamp,
                        "speed_first_um_per_s": first.speed_um_per_s,
                        "speed_second_um_per_s": second.speed_um_per_s,
                    }
                ),
                fontsize=8,
            )
            axis.set_xticks((0, 5, 10, 15))
            axis.set_yticks((0, 5, 10, 15))
            axis.tick_params(labelsize=7)
        for axis in axes.flat[len(matrices) :]:
            axis.set_visible(False)
        if image is not None:
            figure.colorbar(
                image,
                ax=[axis for axis in axes.flat if axis.get_visible()],
                shrink=0.72,
                label="raw force difference, second - first (pN)",
            )
        figure.suptitle(
            f"Same-pixel raw-force difference at D = {target_distance:.0f} nm",
            fontsize=13,
        )
        filename = f"paired_force_difference_maps_{target_distance:.0f}nm.png"
        figure.savefig(figures_dir / filename, dpi=220)
        plt.close(figure)
        files.append(filename)
    return files


def plot_force_differences(
    rows: list[dict], figures_dir: Path
) -> None:
    groups = ((0, "main"), (20, "main"), (30, "high_load"), (40, "main"))
    figure, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for axis, (concentration, regime) in zip(axes.flat, groups):
        group = [
            row
            for row in rows
            if int(row["concentration_wt_percent"]) == concentration
            and str(row["load_regime"]) == regime
            and PAIR_DMIN_NM <= float(row["distance_nm"]) <= PAIR_DMAX_NM
        ]
        identifiers = list(dict.fromkeys(str(row["pair_id"]) for row in group))
        for color_index, identifier in enumerate(identifiers):
            records = sorted(
                [row for row in group if str(row["pair_id"]) == identifier],
                key=lambda row: float(row["distance_nm"]),
            )
            if not records:
                continue
            distance = np.asarray([float(row["distance_nm"]) for row in records])
            force = np.asarray(
                [float(row["median_raw_force_difference_pN"]) for row in records]
            )
            hydro = np.asarray(
                [
                    float(row["theoretical_no_slip_hydrodynamic_difference_pN"])
                    for row in records
                ]
            )
            reference = records[0]
            label = (
                f"{float(reference['speed_first_um_per_s']):.0f}->"
                f"{float(reference['speed_second_um_per_s']):.0f} um/s"
            )
            color = plt.cm.tab10(color_index % 10)
            line_style = "--" if bool(reference["same_speed_control"]) else "-"
            axis.plot(distance, force, line_style, color=color, lw=1.6, label=label)
            if not bool(reference["same_speed_control"]):
                axis.plot(distance, hydro, ":", color=color, lw=1.0, alpha=0.8)
        axis.axhline(0.0, color="black", lw=0.7)
        axis.set_title(f"{concentration} wt%, {regime}")
        axis.set_xlabel("Contact-aligned separation D (nm)")
        axis.set_ylabel("second - first raw force (pN)")
        axis.grid(alpha=0.22)
        axis.legend(fontsize=7, ncol=2)
    figure.suptitle(
        "Same-pixel force differences; dotted curves are no-slip predictions",
        fontsize=13,
    )
    figure.savefig(figures_dir / "paired_force_differences_vs_distance.png", dpi=220)
    plt.close(figure)


def plot_hydrodynamic_coefficients(rows: list[dict], figures_dir: Path) -> None:
    selected = [
        row
        for row in rows
        if bool(row["primary_included"]) and not bool(row["same_speed_control"])
    ]
    if not selected:
        return
    figure, axis = plt.subplots(figsize=(13, 5.5), constrained_layout=True)
    x = np.arange(len(selected))
    chi = np.asarray([float(row["chi_observed_over_no_slip"]) for row in selected])
    low = np.asarray([float(row["chi_local_ci95_low"]) for row in selected])
    high = np.asarray([float(row["chi_local_ci95_high"]) for row in selected])
    error = np.vstack([chi - low, high - chi])
    valid_error = np.all(np.isfinite(error), axis=0) & np.all(error >= 0.0, axis=0)
    axis.scatter(x, chi, c=[int(row["concentration_wt_percent"]) for row in selected], cmap="viridis", s=45)
    if np.any(valid_error):
        axis.errorbar(
            x[valid_error],
            chi[valid_error],
            yerr=error[:, valid_error],
            fmt="none",
            color="0.35",
            capsize=2,
        )
    axis.axhline(0.0, color="black", lw=0.8)
    axis.axhline(1.0, color="#b22222", lw=1.0, ls="--", label="no-slip chi = 1")
    labels = [
        f"{int(row['concentration_wt_percent'])}%\n"
        f"{float(row['speed_first_um_per_s']):.0f}->{float(row['speed_second_um_per_s']):.0f}"
        for row in selected
    ]
    axis.set_xticks(x, labels, rotation=45, ha="right")
    axis.set_ylabel("Fitted chi = observed / no-slip hydrodynamic amplitude")
    axis.set_title("Pairwise 1/D shape test with a simultaneous linear-drift nuisance")
    axis.grid(axis="y", alpha=0.25)
    axis.legend()
    figure.savefig(figures_dir / "pairwise_hydrodynamic_chi.png", dpi=220)
    plt.close(figure)


def plot_branch_support(rows: list[dict], figures_dir: Path) -> None:
    primary = [row for row in rows if bool(row["primary_included"])]
    figure, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    sources = list(dict.fromkeys(str(row["source"]) for row in primary))
    for color_index, source in enumerate(sources):
        records = sorted(
            [row for row in primary if str(row["source"]) == source],
            key=lambda row: float(row["distance_nm"]),
        )
        if not records:
            continue
        distance = np.asarray([float(row["distance_nm"]) for row in records])
        support = np.asarray(
            [float(row["paired_free_branch_support_fraction"]) for row in records]
        )
        odd = np.asarray(
            [float(row["odd_dissipative_candidate_pN"]) for row in records]
        )
        predicted = np.asarray(
            [float(row["theoretical_no_slip_odd_pN"]) for row in records]
        )
        concentration = int(records[0]["concentration_wt_percent"])
        speed = float(records[0]["speed_um_per_s"])
        label = f"{concentration}% {speed:.0f} um/s"
        color = plt.cm.tab20(color_index % 20)
        axes[0].plot(distance, support, color=color, lw=1.2, label=label)
        qualified = support >= MIN_BRANCH_SUPPORT_FRACTION
        axes[1].plot(distance[qualified], odd[qualified], color=color, lw=1.2, label=label)
        axes[1].plot(distance[qualified], predicted[qualified], color=color, lw=0.8, ls=":")
    axes[0].axhline(0.5, color="black", lw=0.8, ls="--")
    axes[0].axvspan(JOINT_DMIN_NM, JOINT_DMAX_NM, color="0.8", alpha=0.25)
    axes[0].set_xlabel("D (nm)")
    axes[0].set_ylabel("free retract / approach paired support")
    axes[0].set_title("Verified free-retract overlap")
    axes[1].axhline(0.0, color="black", lw=0.8)
    axes[1].set_xlabel("D (nm)")
    axes[1].set_ylabel("odd component (pN)")
    axes[1].set_title("Odd component; dotted is no-slip theory")
    for axis in axes:
        axis.grid(alpha=0.22)
    axes[0].legend(fontsize=6, ncol=2)
    figure.savefig(figures_dir / "approach_retract_free_branch_diagnostics.png", dpi=220)
    plt.close(figure)


def plot_joint_fits(
    predictions: list[dict], figures_dir: Path
) -> None:
    selected = [
        row
        for row in predictions
        if str(row["variant"]) == "primary"
        and str(row["velocity_model"]) == "M2_fitted_chi"
        and str(row["scale_mode"]) == "contact_constrained"
    ]
    groups = ((0, "main"), (20, "main"), (30, "high_load"), (40, "main"))
    figure, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for axis, (concentration, regime) in zip(axes.flat, groups):
        group = [
            row
            for row in selected
            if int(row["concentration_wt_percent"]) == concentration
            and str(row["load_regime"]) == regime
        ]
        sources = list(dict.fromkeys(str(row["source"]) for row in group))
        for color_index, source in enumerate(sources):
            records = sorted(
                [row for row in group if str(row["source"]) == source],
                key=lambda row: float(row["model_distance_nm"]),
            )
            if not records:
                continue
            distance = np.asarray(
                [float(row["model_distance_nm"]) for row in records]
            )
            observed = np.asarray([float(row["observed_force_pN"]) for row in records])
            prediction = np.asarray([float(row["predicted_force_pN"]) for row in records])
            speed = float(records[0]["speed_um_per_s"])
            color = plt.cm.tab10(color_index % 10)
            axis.plot(distance, observed, "o", ms=2.2, color=color, alpha=0.55)
            axis.plot(distance, prediction, "-", color=color, lw=1.4, label=f"{speed:.2g} um/s")
        if group:
            reference = sorted(
                group, key=lambda row: float(row["model_distance_nm"])
            )
            distance = np.asarray(
                [float(row["model_distance_nm"]) for row in reference]
            )
            equilibrium = np.asarray(
                [float(row["equilibrium_component_pN"]) for row in reference]
            )
            unique_distance, indices = np.unique(distance, return_index=True)
            axis.plot(
                unique_distance,
                equilibrium[indices],
                "k--",
                lw=1.3,
                label="shared equilibrium component",
            )
        axis.axhline(0.0, color="black", lw=0.7)
        axis.set_title(f"{concentration} wt%, {regime}")
        axis.set_xlabel("D (nm)")
        axis.set_ylabel("raw constant-referenced force (pN)")
        axis.grid(alpha=0.22)
        axis.legend(fontsize=7)
    figure.suptitle("Joint PB + hydrodynamic + map-nuisance M2 fits", fontsize=13)
    figure.savefig(figures_dir / "joint_M2_force_fits.png", dpi=220)
    plt.close(figure)


def plot_joint_summary(fits: list[dict], figures_dir: Path) -> None:
    primary = [
        row
        for row in fits
        if str(row["variant"]) == "primary"
        and np.isfinite(float(row.get("lambda_D_nm", np.nan)))
    ]
    if not primary:
        return
    groups = list(
        dict.fromkeys(
            (int(row["concentration_wt_percent"]), str(row["load_regime"]))
            for row in primary
        )
    )
    models = (
        "M0_no_hydrodynamic",
        "M1_no_slip_chi1",
        "M2_fitted_chi",
    )
    colors = ("#4c78a8", "#f58518", "#54a24b")
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    x = np.arange(len(groups), dtype=float)
    for model_index, (model, color) in enumerate(zip(models, colors)):
        offset = (model_index - 1) * 0.22
        records = [
            next(
                (
                    row
                    for row in primary
                    if int(row["concentration_wt_percent"]) == concentration
                    and str(row["load_regime"]) == regime
                    and str(row["velocity_model"]) == model
                ),
                None,
            )
            for concentration, regime in groups
        ]
        lambdas = [float(row["lambda_D_nm"]) if row else np.nan for row in records]
        zetas = [float(row["zeta_magnitude_mV"]) if row else np.nan for row in records]
        delta_aicc = [
            float(row.get("delta_aicc_within_primary_models", np.nan))
            if row
            else np.nan
            for row in records
        ]
        axes[0].scatter(x + offset, lambdas, color=color, label=model, s=40)
        axes[1].scatter(x + offset, zetas, color=color, label=model, s=40)
        axes[2].scatter(x + offset, delta_aicc, color=color, label=model, s=40)
    labels = [f"{concentration}%\n{regime}" for concentration, regime in groups]
    for axis in axes:
        axis.set_xticks(x, labels)
        axis.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Debye length (nm)")
    axes[1].set_ylabel("Common |zeta| (mV)")
    axes[2].set_ylabel("Delta AICc")
    axes[2].axhline(2.0, color="0.4", ls="--", lw=0.8)
    axes[0].legend(fontsize=7)
    figure.suptitle("Primary joint-model comparison", fontsize=13)
    figure.savefig(figures_dir / "joint_model_parameter_comparison.png", dpi=220)
    plt.close(figure)


def self_test() -> dict[str, float | bool]:
    water_viscosity = base.cheng_viscosity_mPa_s(0.0, base.TEMPERATURE_C)
    force_100 = float(hydrodynamic_force_pN(100.0, 1.0, water_viscosity))
    force_50 = float(hydrodynamic_force_pN(50.0, 1.0, water_viscosity))
    force_double_speed = float(
        hydrodynamic_force_pN(100.0, 2.0, water_viscosity)
    )
    inverse_distance_error = abs(force_50 / force_100 - 2.0)
    linear_speed_error = abs(force_double_speed / force_100 - 2.0)
    if inverse_distance_error > 1e-14 or linear_speed_error > 1e-14:
        raise AssertionError("hydrodynamic limiting scalings failed")

    coordinate = np.linspace(1000.0, 0.0, 1001)
    measured_velocity = smoothed_velocity_um_per_s(
        coordinate, 0.5, closing=True
    )
    expected_velocity = 2.0 * coordinate.size / (coordinate.size - 1)
    velocity_error = float(
        np.max(np.abs(measured_velocity - expected_velocity))
    )
    if velocity_error > 1e-10:
        raise AssertionError("gap velocity recovery failed")

    rng = np.random.default_rng(12345)
    rows = np.repeat(np.arange(16), 16)
    columns = np.tile(np.arange(16), 16)
    constant = np.full(256, 7.25)
    ci_low, ci_high, blocks = block_bootstrap_median_ci(
        rows, columns, constant, rng, replicates=100
    )
    if ci_low != 7.25 or ci_high != 7.25 or blocks != 16:
        raise AssertionError("constant spatial-block bootstrap failed")

    synthetic_rows: list[dict] = []
    synthetic_distance = np.arange(25.0, 255.0, 5.0)
    synthetic_hydro = hydrodynamic_force_pN(
        synthetic_distance, np.full(synthetic_distance.shape, 2.0), water_viscosity
    ) - hydrodynamic_force_pN(
        synthetic_distance, np.full(synthetic_distance.shape, 1.0), water_viscosity
    )
    synthetic_force = (
        7.0 + 0.12 * (synthetic_distance - np.mean(synthetic_distance))
        + 2.5 * synthetic_hydro
    )
    for distance, force, hydro in zip(
        synthetic_distance, synthetic_force, synthetic_hydro
    ):
        synthetic_rows.append(
            {
                "pair_id": "synthetic",
                "concentration_wt_percent": 0,
                "primary_included": True,
                "load_regime": "main",
                "source_first": "a",
                "source_second": "b",
                "speed_first_um_per_s": 1.0,
                "speed_second_um_per_s": 2.0,
                "delta_speed_um_per_s": 1.0,
                "same_speed_control": False,
                "distance_nm": distance,
                "paired_pixel_count": 256,
                "median_raw_force_difference_pN": force,
                "raw_force_difference_block_se_pN": 3.0,
                "theoretical_no_slip_hydrodynamic_difference_pN": hydro,
            }
        )
    synthetic_fit = fit_pair_hydrodynamic_signature(synthetic_rows)[0]
    chi_error = abs(float(synthetic_fit["chi_observed_over_no_slip"]) - 2.5)
    if chi_error > 1e-7:
        raise AssertionError("synthetic hydrodynamic coefficient recovery failed")

    pb_checks = base.self_test()
    return {
        "hydrodynamic_force_water_1um_s_100nm_pN": force_100,
        "hydrodynamic_inverse_distance_relative_error": inverse_distance_error,
        "hydrodynamic_linear_speed_relative_error": linear_speed_error,
        "smoothed_gap_velocity_max_abs_error_um_per_s": velocity_error,
        "spatial_block_constant_ci_exact": True,
        "spatial_block_count_16x16_with_4x4_blocks": blocks,
        "synthetic_pair_chi_absolute_error": chi_error,
        "pb_linear_limit_max_relative_error": pb_checks[
            "pb_linear_limit_max_relative_error"
        ],
        "pb_far_asymptote_max_relative_error": pb_checks[
            "pb_far_asymptote_max_relative_error"
        ],
    }


def _format_number(value: object, digits: int = 3) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "n/a"
    return f"{number:.{digits}g}" if np.isfinite(number) else "n/a"


def measured_force_slice_median(
    map_force_rows: list[dict],
    concentration_wt_percent: int,
    speed_label_um_per_s: float,
    distance_nm: float,
    force_column: str,
) -> float:
    """Return the median eligible map-level force for a report slice.

    The 30 wt% comparison is restricted to the high-load stratum so that its
    1, 2 and 4 um/s entries share the same approximately 18.5 nN load regime.
    Each input row is already a valid-pixel median for one map and distance;
    this function takes the median across eligible same-speed maps when more
    than one exists.
    """

    values = []
    for row in map_force_rows:
        if not bool(row["primary_included"]):
            continue
        if int(row["concentration_wt_percent"]) != concentration_wt_percent:
            continue
        if concentration_wt_percent == 30 and str(row["load_regime"]) != "high_load":
            continue
        if not np.isclose(
            float(row["speed_label_um_per_s"]),
            float(speed_label_um_per_s),
            rtol=0.0,
            atol=1e-12,
        ):
            continue
        if not np.isclose(
            float(row["distance_nm"]),
            float(distance_nm),
            rtol=0.0,
            atol=1e-12,
        ):
            continue
        value = float(row[force_column])
        if np.isfinite(value):
            values.append(value)
    if not values:
        raise AssertionError(
            "missing measured-force report slice for "
            f"c={concentration_wt_percent}, U={speed_label_um_per_s}, "
            f"D={distance_nm}, column={force_column}"
        )
    return float(np.median(np.asarray(values, dtype=np.float64)))


def build_report(
    sources: list[base.SourceData],
    maps: list[MapData],
    map_force_rows: list[dict],
    paired_metric_summaries: list[dict],
    association_rows: list[dict],
    hydro_rows: list[dict],
    branch_summaries: list[dict],
    fits: list[dict],
    systematic_rows: list[dict],
    checks: dict,
) -> str:
    primary_maps = [item for item in maps if item.primary_included]
    same_speed_controls = [
        row
        for row in hydro_rows
        if bool(row["primary_included"]) and bool(row["same_speed_control"])
    ]
    speed_pairs = [
        row
        for row in hydro_rows
        if bool(row["primary_included"]) and not bool(row["same_speed_control"])
    ]
    history_floor = safe_median(
        float(row["null_linear_drift_rmse_pN"]) for row in same_speed_controls
    )
    positive_chi = sum(
        float(row["chi_observed_over_no_slip"]) > 0.0
        for row in speed_pairs
        if np.isfinite(float(row["chi_observed_over_no_slip"]))
    )
    usable_branch_maps = sum(
        bool(row["odd_even_can_constrain_primary_fit_window"])
        for row in branch_summaries
        if bool(row["primary_included"])
    )
    primary_fit_rows = [
        row
        for row in fits
        if str(row["variant"]) == "primary"
        and int(row["concentration_wt_percent"]) in PRIMARY_CONCENTRATIONS
    ]
    robustness_fit_rows = sorted(
        [
            row
            for row in fits
            if int(row["concentration_wt_percent"]) in PRIMARY_CONCENTRATIONS
            and (
                str(row["variant"]).startswith("M2_contact_")
                or str(row["variant"]).startswith("M2_leave_out_")
            )
        ],
        key=lambda row: (
            int(row["concentration_wt_percent"]),
            str(row["variant"]),
        ),
    )
    later_lower_speed = [
        row
        for row in speed_pairs
        if float(row["speed_second_um_per_s"])
        < float(row["speed_first_um_per_s"])
    ]
    later_four_speed = [
        row
        for row in speed_pairs
        if abs(float(row["speed_second_um_per_s"]) - 4.0) < 0.5
    ]
    positive_later_lower = sum(
        float(row["chi_observed_over_no_slip"]) > 0.0
        for row in later_lower_speed
    )
    negative_later_four = sum(
        float(row["chi_observed_over_no_slip"]) < 0.0
        for row in later_four_speed
    )
    lines = [
        "# Same-pixel velocity diagnostics and joint silica sphere-plane fit",
        "",
        "## Result status",
        "",
        f"The raw-only pipeline decoded {len(sources)} JPK sources and retained {len(maps)} maps; {len(primary_maps)} maps at 0, 20, 30 and 40 wt% form the primary analysis. Ten wt% is retained only as a measurement-regime QC control. Every force value uses a concentration-specific consensus calculated from all accepted raw hard-contact anchors and the calibrated cantilever-1 spring constant `k = {base.SPRING_CONSTANT_N_PER_M:.10g} N/m`; embedded JPK sensitivity and force are unused.",
        "",
        "The main conclusion is evidence-bounded: separation-resolved velocity differences are real, but the current sequential-map design does not identify a unique hydrodynamic coefficient independently of map history. Same-speed 30 wt% controls retain a median linear-nuisance residual scale of "
        + f"`{_format_number(history_floor)} pN`, while only `{positive_chi}/{len(speed_pairs)}` nonzero-speed pair fits have a positive fitted no-slip amplitude ratio. The joint PB parameters below are therefore model-conditioned exploratory results and systematic ranges, not a validated zero-velocity extrapolation.",
        "",
        "## Map inventory and load strata",
        "",
        "| wt% | time | nominal speed (um/s) | retained pixels | load regime | terminal load (nN) | contact response |",
        "|---:|:---|---:|---:|:---|---:|---:|",
    ]
    for item in maps:
        lines.append(
            f"| {item.concentration_wt_percent} | {item.timestamp[11:19]} | "
            f"{item.speed_label_um_per_s:g} | {len(item.points)} | {item.load_regime} | "
            f"{item.terminal_load_center_nN:.3g} | {item.contact_response_center:.4f} |"
        )
    lines += [
        "",
        "At 30 wt%, the first two 2 um/s maps are a roughly 6.7 nN low-load stratum. The later 2, 1, 1 and 4 um/s maps are a roughly 18.5 nN high-load stratum. Cross-stratum differences are not used as velocity contrasts.",
        "",
        "## Measured force-distance slices",
        "",
        "The following values are measured approach-force summaries, not PB/hydrodynamic fit components and not a zero-speed extrapolation. Each input value is the valid-pixel median of one map at the stated separation; where multiple eligible maps exist at the same concentration and speed, the table reports the median of those map medians. Ten wt% is excluded. The 30 wt% row uses only the approximately 18.5 nN high-load stratum; the previously quoted approximately 665 pN value at 2 um/s and 50 nm came from the incompatible low-load stratum, whereas the comparable high-load value is 579.4 pN.",
        "",
        "The primary table uses `force_linear_drift_corrected_pN`, obtained after subtracting a fitted far-field linear baseline:",
        "",
        "| D (nm) | wt% / load stratum | 1 um/s (pN) | 2 um/s (pN) | 4 um/s (pN) |",
        "|---:|:---|---:|---:|---:|",
    ]
    for distance_nm in (25.0, 50.0, 100.0, 150.0, 200.0):
        for concentration in PRIMARY_CONCENTRATIONS:
            load_label = (
                f"{concentration} high-load" if concentration == 30 else str(concentration)
            )
            values = [
                measured_force_slice_median(
                    map_force_rows,
                    concentration,
                    speed,
                    distance_nm,
                    "force_linear_drift_corrected_pN",
                )
                for speed in (1.0, 2.0, 4.0)
            ]
            lines.append(
                f"| {distance_nm:.0f} | {load_label} | "
                + " | ".join(f"{value:.1f}" for value in values)
                + " |"
            )
    lines += [
        "",
        "For zero-reference sensitivity, the same 50 nm slices after subtracting only a far-field constant, while retaining the far-field slope, are:",
        "",
        "| wt% / load stratum | 1 um/s (pN) | 2 um/s (pN) | 4 um/s (pN) |",
        "|:---|---:|---:|---:|",
    ]
    for concentration in PRIMARY_CONCENTRATIONS:
        load_label = (
            f"{concentration} high-load" if concentration == 30 else str(concentration)
        )
        values = [
            measured_force_slice_median(
                map_force_rows,
                concentration,
                speed,
                50.0,
                "force_raw_far_constant_referenced_pN",
            )
            for speed in (1.0, 2.0, 4.0)
        ]
        lines.append(
            f"| {load_label} | "
            + " | ".join(f"{value:.1f}" for value in values)
            + " |"
        )
    lines += [
        "",
        "Both tables use the reconstructed concentration-specific sensitivity and calibrated cantilever-1 spring constant. They still contain equilibrium interaction, residual hydrodynamic force, acquisition history and other systematics. The large difference between the line-corrected and constant-referenced branches shows that the far-field slope cannot be silently identified as purely instrumental drift. Because the primary table has had that slope removed, it is not an absolute force-distance relation.",
        "",
        "## Same-pixel event and contact comparisons",
        "",
        "All differences are `second map - first map` at the same physical raster coordinate after undoing serpentine row reversal. Spatial-block confidence intervals use 4x4 pixel blocks and quantify only within-area heterogeneity. Snap differences require a detected event in both maps, so their paired `n` and interval are detection-conditioned; missing/no-snap pixels are not silently assigned zero distance.",
        "",
        "| wt% | load | speed pair (um/s) | metric | paired n | median difference | block CI95 | acquisition-order r | plane R2 |",
        "|---:|:---|:---|:---|---:|---:|:---|---:|---:|",
    ]
    displayed_metrics = {
        "snap_distance_nm",
        "contact_invOLS_nm_per_V",
        "terminal_load_nN",
    }
    for row in paired_metric_summaries:
        if not bool(row["primary_included"]) or row["metric"] not in displayed_metrics:
            continue
        lines.append(
            f"| {row['concentration_wt_percent']} | {row['load_regime']} | "
            f"{float(row['speed_first_um_per_s']):.0f}->{float(row['speed_second_um_per_s']):.0f} | "
            f"{row['metric']} | {row['paired_pixel_count']} | "
            f"{_format_number(row['median_difference_second_minus_first'])} {row['unit']} | "
            f"[{_format_number(row['spatial_block_bootstrap_ci95_low'])}, {_format_number(row['spatial_block_bootstrap_ci95_high'])}] | "
            f"{_format_number(row['difference_correlation_with_acquisition_order'])} | "
            f"{_format_number(row['difference_plane_r2'])} |"
        )
    strongest_associations = sorted(
        [
            row
            for row in association_rows
            if bool(row["primary_included"])
            and np.isfinite(float(row["spearman_rho"]))
        ],
        key=lambda row: abs(float(row["spearman_rho"])),
        reverse=True,
    )[:8]
    lines += [
        "",
        "The strongest per-pair correlations between a snap-position change and a nuisance-metric change are diagnostic, not causal:",
        "",
        "| wt% | load | speed pair | predictor | n | Spearman rho |",
        "|---:|:---|:---|:---|---:|---:|",
    ]
    for row in strongest_associations:
        lines.append(
            f"| {row['concentration_wt_percent']} | {row['load_regime']} | "
            f"{float(row['speed_first_um_per_s']):.0f}->{float(row['speed_second_um_per_s']):.0f} | "
            f"{row['predictor']} | {row['paired_pixel_count']} | {_format_number(row['spearman_rho'])} |"
        )
    lines += [
        "",
        "## Separation-resolved hydrodynamic shape test",
        "",
        "For this test, each trace is referenced by only a far-field voltage constant; its far-field slope is retained. Each paired force difference is fitted as `Delta F = a + b(D-D0) + chi Delta[6 pi eta R^2 U/D]`. Thus the linear instrumental drift and the physical `1/D` shape compete in the same fit instead of the latter being removed beforehand.",
        "",
        "| wt% | load | speed pair | chi | local CI95 | weighted RSS improvement | detrended sign agreement | full RMSE (pN) | status |",
        "|---:|:---|:---|---:|:---|---:|---:|---:|:---|",
    ]
    for row in hydro_rows:
        if not bool(row["primary_included"]):
            continue
        lines.append(
            f"| {row['concentration_wt_percent']} | {row['load_regime']} | "
            f"{float(row['speed_first_um_per_s']):.0f}->{float(row['speed_second_um_per_s']):.0f} | "
            f"{_format_number(row['chi_observed_over_no_slip'])} | "
            f"[{_format_number(row['chi_local_ci95_low'])}, {_format_number(row['chi_local_ci95_high'])}] | "
            f"{_format_number(row['hydrodynamic_term_weighted_rss_improvement_fraction'])} | "
            f"{_format_number(row['detrended_shape_sign_agreement_with_hydrodynamic_difference'])} | "
            f"{_format_number(row['full_hydrodynamic_plus_drift_rmse_pN'])} | "
            f"{row['interpretation']} |"
        )
    lines += [
        "",
        f"The sign pattern is locked to acquisition order rather than a common velocity law: `{positive_later_lower}/{len(later_lower_speed)}` comparisons in which the later map is slower have positive `chi`, whereas `{negative_later_four}/{len(later_four_speed)}` comparisons ending at the chronologically later 4 um/s map have negative `chi`. A no-slip force cannot reverse its fitted sign merely because the speed order is changed; this is direct evidence for a superposed time/surface-history component.",
        "",
        "A physical no-slip result would require positive `chi` near one, a positive nested weighted-RSS improvement beyond the linear drift alone, viscosity/velocity scaling, and an approach/retract sign reversal. Sign agreement is evaluated only after projecting both the measured difference and hydrodynamic basis away from the same intercept/linear-distance nuisance. These pairwise local CIs treat separation bins as independent and are deliberately not promoted to experimental confidence intervals.",
        "",
        "## Approach/retract odd-even feasibility",
        "",
        f"Only `{usable_branch_maps}` primary maps retain at least 12 bins with 50% verified free-retract support inside 25-250 nm. The adhered retract branch is never substituted for a free non-contact branch; encoder-floor values and right-censored detachments remain excluded.",
        "",
        "| wt% | time | speed | minimum D with 50% free support (nm) | qualified bins 25-250 | usable for primary window |",
        "|---:|:---|---:|---:|---:|:---|",
    ]
    for row in branch_summaries:
        if not bool(row["primary_included"]):
            continue
        lines.append(
            f"| {row['concentration_wt_percent']} | {str(row['timestamp'])[11:19]} | "
            f"{float(row['speed_um_per_s']):.2g} | "
            f"{_format_number(row['minimum_distance_with_50pct_free_branch_support_nm'])} | "
            f"{row['qualified_bins_25_250nm']} | {row['odd_even_can_constrain_primary_fit_window']} |"
        )
    lines += [
        "",
        "## Joint same-surface PB fits",
        "",
        "The fitted observation model is `Fobs = gm [Feq(D; |zeta|, lambda_D) + chi Fhyd(D,U)] + am + bm(D-D0)`. `Feq` is nonlinear equal-constant-potential 1:1 PB plus sphere-plane van der Waals force, converted with Derjaguin. Each map's separate raw hard-contact response enters the joint likelihood as an observation of `gm` with a minimum 2% uncertainty; it is not a per-curve force recalibration. Zeta and Debye length are shared across speeds within a concentration/load stratum.",
        "",
        "| wt% | model | lambda_D (nm) | |zeta| (mV) | chi | R2 | RMSE (pN) | Delta AICc | numerical status |",
        "|---:|:---|---:|---:|---:|---:|---:|---:|:---|",
    ]
    for row in primary_fit_rows:
        lines.append(
            f"| {row['concentration_wt_percent']} | {row['velocity_model']} | "
            f"{_format_number(row.get('lambda_D_nm'))} | "
            f"{_format_number(row.get('zeta_magnitude_mV'))} | "
            f"{_format_number(row.get('chi_observed_over_no_slip'))} | "
            f"{_format_number(row.get('r2'), 4)} | {_format_number(row.get('rmse_pN'))} | "
            f"{_format_number(row.get('delta_aicc_within_primary_models'))} | "
            f"{'valid' if row.get('fit_valid') else 'invalid: ' + str(row.get('invalid_reasons',''))} |"
        )
    lines += [
        "",
        "The negative silica sign is assigned chemically; identical-surface normal-force data determine only `|zeta|`. These parameters are boundary-potential fit parameters. Equating them to electrokinetic zeta additionally assumes the boundary potential represents the slipping-plane potential.",
        "",
        "## Model and perturbation systematic ranges",
        "",
        "| wt% | preferred tested primary model | lambda all valid tests (nm) | |zeta| all valid tests (mV) | M2 chi range | claim boundary |",
        "|---:|:---|:---|:---|:---|:---|",
    ]
    for row in systematic_rows:
        lines.append(
            f"| {row['concentration_wt_percent']} | {row['preferred_primary_model']} | "
            f"{_format_number(row['lambda_all_tested_min_nm'])}-{_format_number(row['lambda_all_tested_max_nm'])} | "
            f"{_format_number(row['zeta_all_tested_min_mV'])}-{_format_number(row['zeta_all_tested_max_mV'])} | "
            f"{_format_number(row['M2_chi_min'])}-{_format_number(row['M2_chi_max'])} | "
            f"{row['claim_boundary']} |"
        )
    lines += [
        "",
        "The ranges above include M0/M1/M2, fit-window changes, fixed-unity versus contact-constrained force scale, a ±2.5 nm contact-zero shift, and leave-one-nominal-speed-out M2 fits, retaining only numerically valid results. The contact shift is half one 5 nm analysis bin and is a resolution-scale sensitivity test, not a calibrated contact-position confidence interval.",
        "",
        "### Contact-zero and leave-one-speed-out detail",
        "",
        "| wt% | perturbation | excluded speed (um/s) | contact shift (nm) | lambda_D (nm) | |zeta| (mV) | chi | status |",
        "|---:|:---|---:|---:|---:|---:|---:|:---|",
    ]
    for row in robustness_fit_rows:
        lines.append(
            f"| {row['concentration_wt_percent']} | {row['variant']} | "
            f"{_format_number(row.get('excluded_speed_label_um_per_s'))} | "
            f"{_format_number(row.get('contact_shift_nm'))} | "
            f"{_format_number(row.get('lambda_D_nm'))} | "
            f"{_format_number(row.get('zeta_magnitude_mV'))} | "
            f"{_format_number(row.get('chi_observed_over_no_slip'))} | "
            f"{'valid' if row.get('fit_valid') else 'invalid: ' + str(row.get('invalid_reasons',''))} |"
        )
    lines += [
        "",
        "AICc uses the Gaussian joint weighted likelihood of the force bins and separate map contact-scale observations: because all residuals are normalized by held-fixed uncertainty estimates, the model-dependent `-2 log L` term is the joint weighted chi-square. M2 is genuinely nested with M0 at `chi=0`. The comparison remains conditional on correlated map-median bins and does not remove the experimental speed/time confounding. A preferred model is therefore the best descriptor among M0/M1/M2, not proof that its velocity term is causal.",
        "",
        "## Numerical semantics and checks",
        "",
        f"- Sphere radius: `{base.PROBE_RADIUS_M * 1e6:.6g} um`; primary maximum `D/R = {JOINT_DMAX_NM * 1e-9 / base.PROBE_RADIUS_M:.4f}`, supporting the leading lubrication/Derjaguin small-gap approximation as a controlled diagnostic rather than an exact finite-separation law.",
        f"- Water no-slip force at 1 um/s and 100 nm: `{checks['hydrodynamic_force_water_1um_s_100nm_pN']:.6g} pN`; inverse-distance and linear-speed self-test errors are `{checks['hydrodynamic_inverse_distance_relative_error']:.3e}` and `{checks['hydrodynamic_linear_speed_relative_error']:.3e}`.",
        f"- Synthetic gap-speed maximum error: `{checks['smoothed_gap_velocity_max_abs_error_um_per_s']:.3e} um/s`; synthetic pairwise chi error: `{checks['synthetic_pair_chi_absolute_error']:.3e}`.",
        f"- Nonlinear PB linear-limit and far-asymptote maximum relative errors: `{checks['pb_linear_limit_max_relative_error']:.3e}` and `{checks['pb_far_asymptote_max_relative_error']:.3e}`.",
        "- Optimizer success is supplemented by Jacobian rank/condition, parameter-boundary and R2 checks. Local Jacobian intervals remain conditional because distance bins from one map are correlated.",
        "- All random resampling uses the recorded local seed; every raw ZIP passes CRC and is SHA-256 listed.",
        "",
        "## Files",
        "",
        "- `paired_metric_differences.csv`, `paired_metric_summary.csv`: same-pixel event/contact comparisons and block intervals.",
        "- `snap_nuisance_associations.csv`: paired snap changes versus contact, height, load and far-drift changes.",
        "- `map_force_curves.csv`, `paired_force_differences.csv`: map medians and same-pixel full-curve contrasts.",
        "- `pairwise_hydrodynamic_fits.csv`: `1/D` hydrodynamic shape test with simultaneous linear drift.",
        "- `branch_odd_even.csv`, `branch_odd_even_summary.csv`: verified free-retract support and odd/even candidates.",
        "- `joint_fit_results.csv`, `joint_fit_predictions.csv`, `joint_fit_map_nuisance.csv`: M0/M1/M2 results and components.",
        "- `parameter_systematic_summary.csv`, `map_inventory.csv`, `input_manifest.csv`, figures, provenance and SHA-256 manifest.",
    ]
    return "\n".join(lines) + "\n"


def run(results_dir: Path) -> None:
    print("[1/7] numerical self-tests", flush=True)
    checks = self_test()
    print("[2/7] raw JPK decode, hard-contact sensitivity, paired branches", flush=True)
    sources, maps, event_rows = load_raw_maps()

    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)

    input_rows = [
        {
            "source": str(source.path.relative_to(ROOT)),
            "concentration_wt_percent": source.concentration_wt_percent,
            "source_type": source.source_type,
            "timestamp": source.timestamp.isoformat(),
            "sha256": source.sha256,
            "size_bytes": source.path.stat().st_size,
            "zip_crc_pass": True,
        }
        for source in sources
    ]
    map_rows = [
        {
            "concentration_wt_percent": item.concentration_wt_percent,
            "primary_included": item.primary_included,
            "source": item.source,
            "timestamp": item.timestamp,
            "speed_um_per_s": item.speed_um_per_s,
            "speed_label_um_per_s": item.speed_label_um_per_s,
            "load_regime": item.load_regime,
            "retained_contact_valid_pixels": len(item.points),
            "sensitivity_used_nm_per_V": item.sensitivity_nm_per_V,
            "contact_response_center": item.contact_response_center,
            "contact_response_mad": item.contact_response_mad,
            "contact_invOLS_center_nm_per_V": item.contact_invOLS_center_nm_per_V,
            "terminal_load_center_nN": item.terminal_load_center_nN,
            "far_slope_center_pN_per_100nm": item.far_slope_center_pN_per_100nm,
        }
        for item in maps
    ]
    write_csv(results_dir / "input_manifest.csv", input_rows)
    write_csv(results_dir / "map_inventory.csv", map_rows)
    write_csv(results_dir / "raw_derived_event_metrics.csv", event_rows)

    print("[3/7] same-pixel event and full-force differences", flush=True)
    metric_differences, metric_summaries, association_rows = (
        build_paired_metric_differences(maps, rng)
    )
    map_force_rows = build_map_force_curves(maps)
    paired_force_rows = build_paired_force_differences(maps, rng)
    hydro_rows = fit_pair_hydrodynamic_signature(paired_force_rows)
    write_csv(results_dir / "paired_metric_differences.csv", metric_differences)
    write_csv(results_dir / "paired_metric_summary.csv", metric_summaries)
    write_csv(results_dir / "snap_nuisance_associations.csv", association_rows)
    write_csv(results_dir / "map_force_curves.csv", map_force_rows)
    write_csv(results_dir / "paired_force_differences.csv", paired_force_rows)
    write_csv(results_dir / "pairwise_hydrodynamic_fits.csv", hydro_rows)

    print("[4/7] verified-free-retract odd/even diagnostics", flush=True)
    branch_rows, branch_summaries = build_branch_odd_even(maps)
    write_csv(results_dir / "branch_odd_even.csv", branch_rows)
    write_csv(results_dir / "branch_odd_even_summary.csv", branch_summaries)

    print("[5/7] joint nonlinear-PB / hydrodynamic / nuisance fits", flush=True)
    fits, predictions, nuisance_rows, systematic_rows = run_joint_models(
        maps, map_force_rows
    )
    write_csv(results_dir / "joint_fit_results.csv", fits)
    write_csv(results_dir / "joint_fit_predictions.csv", predictions)
    write_csv(results_dir / "joint_fit_map_nuisance.csv", nuisance_rows)
    write_csv(results_dir / "parameter_systematic_summary.csv", systematic_rows)

    print("[6/7] figures", flush=True)
    plot_paired_metric_maps(metric_differences, figures_dir)
    plot_selected_force_difference_maps(maps, figures_dir)
    plot_force_differences(paired_force_rows, figures_dir)
    plot_hydrodynamic_coefficients(hydro_rows, figures_dir)
    plot_branch_support(branch_rows, figures_dir)
    plot_joint_fits(predictions, figures_dir)
    plot_joint_summary(fits, figures_dir)

    provenance = {
        "analysis_script": str(Path(__file__).relative_to(ROOT)),
        "analysis_script_sha256": sha256_file(Path(__file__)),
        "dependency_scripts": {
            "fit_glycerol_surface_forces.py": sha256_file(
                ROOT / "analysis" / "fit_glycerol_surface_forces.py"
            ),
            "analyze_velocity_systematics.py": sha256_file(
                ROOT / "analysis" / "analyze_velocity_systematics.py"
            ),
        },
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "random_seed": RANDOM_SEED,
        "bootstrap_replicates_event_metrics": BOOTSTRAP_REPLICATES,
        "bootstrap_replicates_force_bins": 500,
        "spatial_block_size_pixels": SPATIAL_BLOCK_SIZE,
        "temperature_C": base.TEMPERATURE_C,
        "temperature_K": base.TEMPERATURE_K,
        "spring_constant_N_per_m": base.SPRING_CONSTANT_N_PER_M,
        "probe_radius_m": base.PROBE_RADIUS_M,
        "hamaker_J": base.HAMAKER_J,
        "epsilon_r": base.EPSILON_R,
        "viscosity_model": "Cheng glycerol-water mass-fraction correlation",
        "jpk_embedded_sensitivity_used": False,
        "jpk_embedded_force_used": False,
        "sensitivity_method": "equal-source-weighted concentration consensus from all accepted raw hard contacts",
        "force_reconstruction": "F=k*S_concentration*(V-Vreference)",
        "distance_definition": "D=measuredHeight+delta-h_contact using the preliminary robust far-linear baseline for delta/contact alignment",
        "raw_force_pairing_baseline": "far-field median voltage constant only; slope retained and fitted jointly",
        "map_experimental_unit": True,
        "pixel_role": "same-location pairing and within-area spatial-block uncertainty only",
        "qc_excluded_concentration_wt_percent": QC_EXCLUDED_CONCENTRATION,
        "thirty_wt_load_split_nN": 12.0,
        "joint_fit_window_nm": [JOINT_DMIN_NM, JOINT_DMAX_NM],
        "contact_zero_sensitivity_shifts_nm": [-2.5, 2.5],
        "contact_zero_shift_convention": "D_model=D_observed+shift; fit-window selection is applied to D_model",
        "leave_one_speed_out_labels_um_per_s": [1.0, 2.0, 4.0],
        "hydrodynamic_model": "leading no-slip sphere-plane lubrication 6*pi*eta*R^2*U/D",
        "hydrodynamic_speed": "Savitzky-Golay derivative of contact-aligned D with implausible local derivatives replaced by robust measuredHeight speed",
        "joint_force_model": "g_map*(nonlinear equal-potential PB Derjaguin + van der Waals + chi*hydrodynamic)+a_map+b_map*(D-Dcenter)",
        "contact_scale_observation_sigma_floor_fraction": 0.02,
        "model_comparison_score": "Gaussian AICc with -2logL equal to joint weighted chi-square of force bins and separate map contact-scale observations; their uncertainty scales are held fixed",
        "retract_rule": "only verified unsaturated free side before detected snap-off; adhered and censored branches excluded",
        "local_interval_limit": "distance bins correlated and map-speed sequentially confounded; local Jacobian intervals are not experimental confidence intervals",
        "approximations": {
            "primary_max_D_over_R": JOINT_DMAX_NM * 1e-9 / base.PROBE_RADIUS_M,
            "max_reynolds_estimate_rho1000_U4um_s_R_eta_water": 1000.0
            * 4e-6
            * base.PROBE_RADIUS_M
            / (base.cheng_viscosity_mPa_s(0.0, base.TEMPERATURE_C) * 1e-3),
        },
        "self_checks": checks,
    }
    (results_dir / "provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    report = build_report(
        sources,
        maps,
        map_force_rows,
        metric_summaries,
        association_rows,
        hydro_rows,
        branch_summaries,
        fits,
        systematic_rows,
        checks,
    )
    (results_dir / "REPORT.md").write_text(report, encoding="utf-8")

    print("[7/7] provenance and artifact hashes", flush=True)
    artifacts = sorted(
        path
        for path in results_dir.rglob("*")
        if path.is_file() and path.name != "artifact_manifest.sha256"
    )
    manifest = "\n".join(
        f"{sha256_file(path)}  {path.relative_to(results_dir)}"
        for path in artifacts
    )
    (results_dir / "artifact_manifest.sha256").write_text(
        manifest + "\n", encoding="utf-8"
    )
    print(f"complete: {results_dir}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--self-test", action="store_true")
    arguments = parser.parse_args()
    if arguments.self_test:
        print(json.dumps(self_test(), indent=2))
        return
    run(arguments.results_dir.resolve())


if __name__ == "__main__":
    main()

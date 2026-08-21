#!/usr/bin/env python3
"""Diagnose velocity-correlated AFM systematics from paired raw branches.

This analysis deliberately starts again from raw JPK vDeflection and
measuredHeight channels.  It recomputes the concentration-specific InvOLS
consensus from raw approach contacts, applies the calibrated cantilever-1
spring constant, and pairs every approach trace with its retract trace.

The primary experimental unit for velocity comparisons is a complete map, not
an individual pixel.  The 10 wt% condition is retained as anomaly/QC evidence
but excluded from primary cross-velocity interpretation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.ndimage import median_filter

import fit_glycerol_surface_forces as base


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT / "20-08-26"
DEFAULT_RESULTS = ROOT / "analysis" / "velocity_systematics_results"
PRIMARY_CONCENTRATIONS = (0, 20, 30, 40)
QC_EXCLUDED_CONCENTRATION = 10

STEP_LAG = 2
SMOOTH_POINTS = 5
APPROACH_STEP_SIGMA = 8.0
RETRACT_STEP_SIGMA = 10.0
APPROACH_ABSOLUTE_STEP_PN = 50.0
RETRACT_ABSOLUTE_STEP_PN = 25.0
MIN_MAP_EVENT_SUPPORT_FRACTION = 0.50


@dataclass
class BranchCurve:
    point_index: int | None
    measured_height_m: np.ndarray
    deflection_V: np.ndarray
    duration_s: float
    raw_saturation_mask: np.ndarray
    raw_low_saturation_fraction: float
    raw_high_saturation_fraction: float
    far_raw_saturation_fraction: float
    first_sample_saturated: bool
    leading_unsaturated_points: int
    raw_min: float
    raw_max: float


@dataclass
class StepEvent:
    candidate_found: bool
    detected: bool
    left_index: int
    right_index: int
    step_pN: float
    threshold_pN: float
    noise_pN: float
    snr: float


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
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _branch_prefixes(
    names: Iterable[str],
    shared: Mapping[str, str],
    archive: zipfile.ZipFile,
    source_type: str,
    wanted_style: str,
) -> list[tuple[int | None, str]]:
    if source_type == "force":
        pattern = re.compile(r"^segments/(\d+)/segment-header\.properties$")
        candidates = [
            (None, int(match.group(1)), name.rsplit("/", 1)[0])
            for name in names
            if (match := pattern.match(name))
        ]
    else:
        pattern = re.compile(
            r"^index/(\d+)/segments/(\d+)/segment-header\.properties$"
        )
        candidates = [
            (
                int(match.group(1)),
                int(match.group(2)),
                name.rsplit("/", 1)[0],
            )
            for name in names
            if (match := pattern.match(name))
        ]
    output: list[tuple[int | None, str]] = []
    for point_index, segment_index, prefix in sorted(
        candidates,
        key=lambda item: (-1 if item[0] is None else item[0], item[1]),
    ):
        del segment_index
        segment = base.parse_properties(
            archive.read(f"{prefix}/segment-header.properties")
        )
        if base._segment_style(segment, shared) == wanted_style.lower():
            output.append((point_index, prefix))
    return output


def load_branch(path: Path, wanted_style: str) -> tuple[list[BranchCurve], int]:
    """Decode one branch and normalize it to far-to-contact orientation."""

    source_type = "map" if path.suffix == ".jpk-force-map" else "force"
    curves: list[BranchCurve] = []
    skipped = 0
    with zipfile.ZipFile(path, "r") as archive:
        corrupt = archive.testzip()
        if corrupt is not None:
            raise ValueError(f"CRC failure in {path}: {corrupt}")
        shared = base.parse_properties(
            archive.read("shared-data/header.properties")
        )
        prefixes = _branch_prefixes(
            archive.namelist(), shared, archive, source_type, wanted_style
        )
        for point_index, prefix in prefixes:
            try:
                segment = base.parse_properties(
                    archive.read(f"{prefix}/segment-header.properties")
                )
                deflection, _ = base.decode_base_channel(
                    archive, shared, segment, prefix, "vDeflection"
                )
                v_lcd = int(base._required(segment, "channel.vDeflection.lcd-info.*"))
                v_lcd_prefix = f"lcd-info.{v_lcd}"
                v_member = (
                    f"{prefix}/"
                    f"{base._required(segment, 'channel.vDeflection.data.file.name')}"
                )
                raw_dtype = base._raw_dtype(
                    base._required(shared, f"{v_lcd_prefix}.encoder.type")
                )
                raw_deflection = np.frombuffer(
                    archive.read(v_member), dtype=raw_dtype
                ).copy()
                measured_base, height_lcd = base.decode_base_channel(
                    archive, shared, segment, prefix, "measuredHeight"
                )
                measured_height = base.apply_linear_conversion(
                    measured_base, shared, height_lcd, "nominal"
                )
                duration = float(
                    base._required(segment, "force-segment-header.duration")
                )
                if (
                    deflection.shape != measured_height.shape
                    or deflection.size < 200
                    or not np.isfinite(duration)
                    or duration <= 0.0
                ):
                    raise ValueError("short or inconsistent branch")
                # The rest of the analysis uses one orientation: the first
                # sample is far from the surface and the last is in contact.
                if measured_height[0] < measured_height[-1]:
                    measured_height = measured_height[::-1].copy()
                    deflection = deflection[::-1].copy()
                    raw_deflection = raw_deflection[::-1].copy()
                if np.issubdtype(raw_dtype, np.integer):
                    limits = np.iinfo(raw_dtype)
                    # JPK's observed signed-int32 clamp is INT32_MIN + 48.
                    # A small integer margin catches that clamp without
                    # classifying merely large physical signals as saturated.
                    margin = max(
                        1024,
                        int((int(limits.max) - int(limits.min)) * 1e-6),
                    )
                    raw_i64 = raw_deflection.astype(np.int64)
                    low_saturation = raw_i64 <= int(limits.min) + margin
                    high_saturation = raw_i64 >= int(limits.max) - margin
                    saturation = low_saturation | high_saturation
                else:
                    low_saturation = np.zeros(raw_deflection.shape, dtype=bool)
                    high_saturation = np.zeros(raw_deflection.shape, dtype=bool)
                    saturation = np.zeros(raw_deflection.shape, dtype=bool)
                far_points = min(
                    saturation.size,
                    max(
                        base.FAR_FIELD_MIN_POINTS,
                        int(math.ceil(base.FAR_FIELD_FRACTION * saturation.size)),
                    ),
                )
                saturated_indices = np.flatnonzero(saturation)
                leading_unsaturated = (
                    int(saturated_indices[0])
                    if saturated_indices.size
                    else int(saturation.size)
                )
                curves.append(
                    BranchCurve(
                        point_index=point_index,
                        measured_height_m=measured_height,
                        deflection_V=deflection,
                        duration_s=duration,
                        raw_saturation_mask=saturation,
                        raw_low_saturation_fraction=float(np.mean(low_saturation)),
                        raw_high_saturation_fraction=float(np.mean(high_saturation)),
                        far_raw_saturation_fraction=float(
                            np.mean(saturation[:far_points])
                        ),
                        first_sample_saturated=bool(saturation[0]),
                        leading_unsaturated_points=leading_unsaturated,
                        raw_min=float(np.min(raw_deflection)),
                        raw_max=float(np.max(raw_deflection)),
                    )
                )
            except (KeyError, ValueError, zipfile.BadZipFile):
                skipped += 1
    return curves, skipped


def estimate_speed_um_per_s(height_m: np.ndarray, duration_s: float) -> float:
    count = height_m.size
    time_s = (np.arange(count, dtype=np.float64) + 0.5) * duration_s / count
    trim = max(5, int(math.ceil(0.10 * count)))
    slope, _, _, _ = base.robust_line(
        time_s[trim:-trim], height_m[trim:-trim]
    )
    return abs(float(slope)) * 1e6


def corrected_branch(
    curve: base.RawCurve | BranchCurve,
    sensitivity_m_per_V: float,
    far_fit: base.FarFieldDriftFit | None = None,
) -> tuple[base.FarFieldDriftFit, np.ndarray, np.ndarray]:
    selected_fit = (
        far_fit
        if far_fit is not None
        else base.fit_far_field_drift(curve.measured_height_m, curve.deflection_V)
    )
    corrected_V = curve.deflection_V - base.baseline_voltage(
        curve.measured_height_m, curve.deflection_V, selected_fit
    )
    delta_m = sensitivity_m_per_V * corrected_V
    force_pN = base.SPRING_CONSTANT_N_PER_M * delta_m * 1e12
    return selected_fit, delta_m, force_pN


def fit_far_field_free_prefix(
    curve: BranchCurve,
) -> base.FarFieldDriftFit:
    """Fit only the verified unsaturated free prefix of a retract trace."""

    default_points = min(
        curve.measured_height_m.size,
        max(
            base.FAR_FIELD_MIN_POINTS,
            int(math.ceil(base.FAR_FIELD_FRACTION * curve.measured_height_m.size)),
        ),
    )
    points = min(default_points, curve.leading_unsaturated_points)
    if points < base.FAR_FIELD_MIN_POINTS:
        raise ValueError("retract has no sufficiently long unsaturated free prefix")
    height = curve.measured_height_m[:points]
    voltage = curve.deflection_V[:points]
    slope, intercept, r2, n_inliers = base.robust_line(height, voltage)
    residual = voltage - (slope * height + intercept)
    return base.FarFieldDriftFit(
        slope_V_per_m=float(slope),
        intercept_V=float(intercept),
        r2=float(r2),
        n_points=int(points),
        n_inliers=int(n_inliers),
        height_min_m=float(np.min(height)),
        height_max_m=float(np.max(height)),
        residual_mad_V=base.robust_mad(residual),
    )


def contact_height_m(
    height_m: np.ndarray,
    corrected_V: np.ndarray,
    contact_fit: base.ContactFit,
) -> float:
    slope, intercept, _, _ = base.robust_line(
        height_m[contact_fit.start : contact_fit.stop],
        corrected_V[contact_fit.start : contact_fit.stop],
    )
    if not np.isfinite(slope) or slope == 0.0:
        raise ValueError("invalid contact slope")
    return float(-intercept / slope)


def detect_negative_step(
    force_pN: np.ndarray,
    search_mask: np.ndarray,
    far_points: int,
    residual_noise_pN: float,
    sigma_multiplier: float,
    absolute_floor_pN: float,
) -> StepEvent:
    values = np.asarray(force_pN, dtype=np.float64)
    mask = np.asarray(search_mask, dtype=bool)
    if values.ndim != 1 or mask.shape != (values.size - STEP_LAG,):
        raise ValueError("step-search dimensions are inconsistent")
    smooth = median_filter(values, size=SMOOTH_POINTS, mode="nearest")
    steps = smooth[STEP_LAG:] - smooth[:-STEP_LAG]
    far_stop = min(max(STEP_LAG + 10, int(far_points)), values.size)
    far_steps = (
        smooth[STEP_LAG:far_stop] - smooth[: far_stop - STEP_LAG]
    )
    noise = base.robust_mad(far_steps)
    if not np.isfinite(noise) or noise <= np.finfo(np.float64).tiny:
        noise = float(np.nanstd(far_steps))
    noise = max(float(noise), np.finfo(np.float64).tiny)
    residual_term = (
        5.0 * float(residual_noise_pN)
        if np.isfinite(residual_noise_pN)
        else 0.0
    )
    threshold = max(
        sigma_multiplier * noise, residual_term, float(absolute_floor_pN)
    )
    candidates = np.flatnonzero(mask & np.isfinite(steps))
    if not candidates.size:
        return StepEvent(False, False, -1, -1, np.nan, threshold, noise, np.nan)
    left = int(candidates[np.argmin(steps[candidates])])
    step = float(steps[left])
    right = left + STEP_LAG
    return StepEvent(
        candidate_found=True,
        detected=bool(step < -threshold),
        left_index=left,
        right_index=right,
        step_pN=step,
        threshold_pN=float(threshold),
        noise_pN=float(noise),
        snr=float(abs(step) / noise),
    )


def _nan_record() -> dict[str, float | bool]:
    fields = (
        "approach_contact_invOLS_nm_per_V",
        "approach_contact_r2",
        "approach_contact_response_ratio",
        "approach_contact_residual_pN",
        "approach_contact_height_um",
        "approach_terminal_load_nN",
        "approach_far_slope_pN_per_100nm",
        "approach_far_noise_pN",
        "approach_snap_step_pN",
        "approach_snap_threshold_pN",
        "approach_snap_snr",
        "approach_snap_distance_nm",
        "approach_snap_scanner_gap_nm",
        "approach_snap_force_before_nN",
        "approach_snap_force_after_nN",
        "approach_barrier_force_nN",
        "retract_contact_invOLS_nm_per_V",
        "retract_contact_r2",
        "retract_contact_response_ratio",
        "retract_contact_residual_pN",
        "retract_contact_height_um",
        "branch_contact_height_difference_nm",
        "retract_terminal_load_nN",
        "retract_far_slope_pN_per_100nm",
        "retract_far_noise_pN",
        "retract_snapoff_step_pN",
        "retract_snapoff_threshold_pN",
        "retract_snapoff_snr",
        "retract_detachment_piezo_travel_nm",
        "retract_detachment_free_side_distance_nm",
        "retract_detachment_adhesive_side_distance_nm",
        "retract_min_precontact_force_nN",
        "retract_pull_off_force_nN",
        "retract_pull_off_observed_floor_nN",
        "retract_retraction_range_nm",
    )
    record: dict[str, float | bool] = {field: np.nan for field in fields}
    record.update(
        {
            "approach_contact_valid": False,
            "approach_snap_candidate_found": False,
            "approach_snap_detected": False,
            "retract_branch_present": False,
            "retract_free_baseline_valid": False,
            "retract_contact_valid": False,
            "retract_snapoff_candidate_found": False,
            "retract_snapoff_detected": False,
            "retract_pull_off_censored": False,
            "retract_detachment_right_censored": False,
        }
    )
    return record


def analyze_pair(
    approach: base.RawCurve,
    retract: BranchCurve | None,
    sensitivity_m_per_V: float,
) -> dict:
    record = _nan_record()
    app_far, app_delta, app_force = corrected_branch(
        approach, sensitivity_m_per_V, approach.far_field_fit
    )
    app_corrected = approach.deflection_V - base.baseline_voltage(
        approach.measured_height_m, approach.deflection_V, app_far
    )
    record["approach_far_slope_pN_per_100nm"] = (
        base.SPRING_CONSTANT_N_PER_M
        * sensitivity_m_per_V
        * app_far.slope_V_per_m
        * 1e5
    )
    app_residual_pN = (
        base.SPRING_CONSTANT_N_PER_M
        * sensitivity_m_per_V
        * app_far.residual_mad_V
        * 1e12
    )
    record["approach_far_noise_pN"] = app_residual_pN
    record["approach_terminal_load_nN"] = float(
        np.median(app_force[-10:]) / 1e3
    )

    app_contact = approach.contact_fit
    app_contact_height = np.nan
    app_distance_nm = np.full(app_force.shape, np.nan)
    if app_contact is not None:
        app_contact_height = contact_height_m(
            approach.measured_height_m, app_corrected, app_contact
        )
        app_distance_nm = (
            approach.measured_height_m + app_delta - app_contact_height
        ) * 1e9
        contact_prediction = (
            app_contact.slope_V_per_m
            * approach.measured_height_m[app_contact.start : app_contact.stop]
            + app_contact.intercept_V
        )
        contact_residual = (
            app_corrected[app_contact.start : app_contact.stop]
            - contact_prediction
        )
        record.update(
            {
                "approach_contact_valid": True,
                "approach_contact_invOLS_nm_per_V": app_contact.sensitivity_m_per_V
                * 1e9,
                "approach_contact_r2": app_contact.r2,
                "approach_contact_response_ratio": sensitivity_m_per_V
                / app_contact.sensitivity_m_per_V,
                "approach_contact_residual_pN": base.SPRING_CONSTANT_N_PER_M
                * sensitivity_m_per_V
                * base.robust_mad(contact_residual)
                * 1e12,
                "approach_contact_height_um": app_contact_height * 1e6,
            }
        )
        left_distance = app_distance_nm[:-STEP_LAG]
        app_search = (
            np.arange(app_force.size - STEP_LAG) < app_contact.start - STEP_LAG
        ) & (left_distance >= 0.0) & (left_distance <= 350.0)
        app_event = detect_negative_step(
            app_force,
            app_search,
            app_far.n_points,
            app_residual_pN,
            APPROACH_STEP_SIGMA,
            APPROACH_ABSOLUTE_STEP_PN,
        )
        record.update(
            {
                "approach_snap_candidate_found": app_event.candidate_found,
                "approach_snap_detected": app_event.detected,
                "approach_snap_step_pN": app_event.step_pN,
                "approach_snap_threshold_pN": app_event.threshold_pN,
                "approach_snap_snr": app_event.snr,
            }
        )
        if app_event.candidate_found:
            left = app_event.left_index
            right = app_event.right_index
            record.update(
                {
                    "approach_snap_distance_nm": app_distance_nm[left],
                    "approach_snap_scanner_gap_nm": (
                        approach.measured_height_m[left] - app_contact_height
                    )
                    * 1e9,
                    "approach_snap_force_before_nN": app_force[left] / 1e3,
                    "approach_snap_force_after_nN": app_force[right] / 1e3,
                    "approach_barrier_force_nN": float(
                        np.max(app_force[: right + 1]) / 1e3
                    ),
                }
            )

    if retract is None:
        return record
    record["retract_branch_present"] = True
    record["retract_retraction_range_nm"] = float(
        np.ptp(retract.measured_height_m) * 1e9
    )
    if (
        retract.first_sample_saturated
        or retract.leading_unsaturated_points < base.FAR_FIELD_MIN_POINTS
    ):
        # There is no verified free-cantilever baseline at the end of the
        # physical retract.  Detachment therefore was not observed within the
        # available travel (or the signal was already outside encoder range).
        record["retract_detachment_right_censored"] = True
        record["retract_pull_off_censored"] = True
        return record
    try:
        ret_far_fit = fit_far_field_free_prefix(retract)
        ret_far, ret_delta, ret_force = corrected_branch(
            retract, sensitivity_m_per_V, ret_far_fit
        )
    except (ValueError, np.linalg.LinAlgError):
        return record
    record["retract_free_baseline_valid"] = True
    ret_corrected = retract.deflection_V - base.baseline_voltage(
        retract.measured_height_m, retract.deflection_V, ret_far
    )
    ret_contact = base.find_contact_fit(
        retract.measured_height_m, retract.deflection_V, ret_far
    )
    ret_residual_pN = (
        base.SPRING_CONSTANT_N_PER_M
        * sensitivity_m_per_V
        * ret_far.residual_mad_V
        * 1e12
    )
    record.update(
        {
            "retract_far_slope_pN_per_100nm": base.SPRING_CONSTANT_N_PER_M
            * sensitivity_m_per_V
            * ret_far.slope_V_per_m
            * 1e5,
            "retract_far_noise_pN": ret_residual_pN,
            "retract_terminal_load_nN": float(np.median(ret_force[-10:]) / 1e3),
        }
    )
    if ret_contact is None:
        return record
    ret_contact_height = contact_height_m(
        retract.measured_height_m, ret_corrected, ret_contact
    )
    contact_prediction = (
        ret_contact.slope_V_per_m
        * retract.measured_height_m[ret_contact.start : ret_contact.stop]
        + ret_contact.intercept_V
    )
    contact_residual = (
        ret_corrected[ret_contact.start : ret_contact.stop] - contact_prediction
    )
    record.update(
        {
            "retract_contact_valid": True,
            "retract_contact_invOLS_nm_per_V": ret_contact.sensitivity_m_per_V
            * 1e9,
            "retract_contact_r2": ret_contact.r2,
            "retract_contact_response_ratio": sensitivity_m_per_V
            / ret_contact.sensitivity_m_per_V,
            "retract_contact_residual_pN": base.SPRING_CONSTANT_N_PER_M
            * sensitivity_m_per_V
            * base.robust_mad(contact_residual)
            * 1e12,
            "retract_contact_height_um": ret_contact_height * 1e6,
        }
    )
    if np.isfinite(app_contact_height):
        record["branch_contact_height_difference_nm"] = (
            ret_contact_height - app_contact_height
        ) * 1e9
        ret_distance_nm = (
            retract.measured_height_m + ret_delta - app_contact_height
        ) * 1e9
        ret_piezo_travel_nm = (
            retract.measured_height_m - app_contact_height
        ) * 1e9
    else:
        ret_distance_nm = (
            retract.measured_height_m + ret_delta - ret_contact_height
        ) * 1e9
        ret_piezo_travel_nm = (
            retract.measured_height_m - ret_contact_height
        ) * 1e9
    left_travel = ret_piezo_travel_nm[:-STEP_LAG]
    ret_search = (
        np.arange(ret_force.size - STEP_LAG) < ret_contact.start - STEP_LAG
    ) & (left_travel >= -50.0) & (left_travel <= 1500.0)
    ret_event = detect_negative_step(
        ret_force,
        ret_search,
        ret_far.n_points,
        ret_residual_pN,
        RETRACT_STEP_SIGMA,
        RETRACT_ABSOLUTE_STEP_PN,
    )
    precontact = ret_force[: ret_contact.start]
    record["retract_min_precontact_force_nN"] = (
        float(np.min(precontact) / 1e3) if precontact.size else np.nan
    )
    record.update(
        {
            "retract_snapoff_candidate_found": ret_event.candidate_found,
            "retract_snapoff_detected": ret_event.detected,
            "retract_snapoff_step_pN": ret_event.step_pN,
            "retract_snapoff_threshold_pN": ret_event.threshold_pN,
            "retract_snapoff_snr": ret_event.snr,
        }
    )
    if ret_event.candidate_found:
        free = ret_event.left_index
        adhesive = ret_event.right_index
        record.update(
            {
                "retract_detachment_piezo_travel_nm": ret_piezo_travel_nm[
                    adhesive
                ],
                "retract_detachment_free_side_distance_nm": ret_distance_nm[free],
                "retract_detachment_adhesive_side_distance_nm": ret_distance_nm[
                    adhesive
                ],
            }
        )
        adhesive_branch = ret_force[adhesive : ret_contact.start]
        if adhesive_branch.size:
            observed_floor = float(np.min(adhesive_branch) / 1e3)
            record["retract_pull_off_observed_floor_nN"] = observed_floor
            adhesive_saturation = bool(
                np.any(retract.raw_saturation_mask[adhesive : ret_contact.start])
            )
            record["retract_pull_off_censored"] = adhesive_saturation
            if not adhesive_saturation:
                record["retract_pull_off_force_nN"] = observed_floor
    return record


def _curve_key(curve: base.RawCurve | BranchCurve, order: int) -> int:
    return int(curve.point_index) if curve.point_index is not None else int(order)


def pair_curves(
    approaches: list[base.RawCurve], retracts: list[BranchCurve]
) -> list[tuple[base.RawCurve, BranchCurve | None]]:
    retract_by_key = {
        _curve_key(curve, order): curve for order, curve in enumerate(retracts)
    }
    return [
        (curve, retract_by_key.get(_curve_key(curve, order)))
        for order, curve in enumerate(approaches)
    ]


SUMMARY_METRICS = (
    "approach_contact_invOLS_nm_per_V",
    "approach_contact_response_ratio",
    "approach_contact_residual_pN",
    "approach_contact_height_um",
    "approach_terminal_load_nN",
    "approach_far_slope_pN_per_100nm",
    "approach_far_noise_pN",
    "approach_snap_distance_nm",
    "approach_snap_step_pN",
    "approach_barrier_force_nN",
    "retract_contact_invOLS_nm_per_V",
    "retract_contact_response_ratio",
    "retract_terminal_load_nN",
    "retract_far_slope_pN_per_100nm",
    "retract_far_noise_pN",
    "branch_contact_height_difference_nm",
    "retract_detachment_piezo_travel_nm",
    "retract_detachment_free_side_distance_nm",
    "retract_pull_off_force_nN",
    "retract_pull_off_observed_floor_nN",
    "retract_min_precontact_force_nN",
    "retract_retraction_range_nm",
)


TREND_METRICS = (
    "approach_contact_invOLS_nm_per_V",
    "approach_contact_response_ratio",
    "approach_contact_height_um",
    "approach_terminal_load_nN",
    "approach_far_slope_pN_per_100nm",
    "approach_snap_distance_nm",
    "approach_barrier_force_nN",
    "retract_contact_invOLS_nm_per_V",
    "retract_detachment_piezo_travel_nm",
    "retract_pull_off_force_nN",
)


def finite_values(rows: list[dict], field: str) -> np.ndarray:
    return np.asarray(
        [float(row[field]) for row in rows if np.isfinite(float(row[field]))],
        dtype=np.float64,
    )


def source_summary(source: base.SourceData, rows: list[dict]) -> dict:
    speed_app = finite_values(rows, "approach_speed_um_per_s")
    speed_ret = finite_values(rows, "retract_speed_um_per_s")
    summary: dict = {
        "concentration_wt_percent": source.concentration_wt_percent,
        "primary_included": source.concentration_wt_percent
        in PRIMARY_CONCENTRATIONS,
        "source": str(source.path.relative_to(ROOT)),
        "source_type": source.source_type,
        "timestamp": source.timestamp.isoformat(),
        "curve_count": len(rows),
        "paired_retract_count": int(
            sum(bool(row["retract_branch_present"]) for row in rows)
        ),
        "approach_speed_um_per_s": float(np.median(speed_app))
        if speed_app.size
        else np.nan,
        "retract_speed_um_per_s": float(np.median(speed_ret))
        if speed_ret.size
        else np.nan,
        "speed_label_um_per_s": base._speed_label(float(np.median(speed_app)))
        if speed_app.size
        else np.nan,
        "sensitivity_used_nm_per_V": source.sensitivity_used_m_per_V * 1e9,
        "sensitivity_anchor_nm_per_V": source.sensitivity_anchor_m_per_V * 1e9,
        "sensitivity_anchor_mad_nm_per_V": source.sensitivity_anchor_mad_m_per_V
        * 1e9,
        "approach_contact_valid_fraction": float(
            np.mean([bool(row["approach_contact_valid"]) for row in rows])
        ),
        "approach_snap_detected_fraction": float(
            np.mean([bool(row["approach_snap_detected"]) for row in rows])
        ),
        "retract_contact_valid_fraction": float(
            np.mean([bool(row["retract_contact_valid"]) for row in rows])
        ),
        "approach_raw_saturation_fraction": float(
            np.mean([bool(row["approach_raw_saturated"]) for row in rows])
        ),
        "retract_raw_saturation_fraction": float(
            np.mean([bool(row["retract_raw_saturated"]) for row in rows])
        ),
        "retract_free_baseline_valid_fraction": float(
            np.mean([bool(row["retract_free_baseline_valid"]) for row in rows])
        ),
        "retract_snapoff_detected_fraction": float(
            np.mean([bool(row["retract_snapoff_detected"]) for row in rows])
        ),
        "retract_pull_off_censored_fraction": float(
            np.mean([bool(row["retract_pull_off_censored"]) for row in rows])
        ),
        "retract_detachment_right_censored_fraction": float(
            np.mean(
                [bool(row["retract_detachment_right_censored"]) for row in rows]
            )
        ),
    }
    for field in SUMMARY_METRICS:
        eligible_rows = rows
        if field.startswith("approach_"):
            eligible_rows = [
                row for row in rows if not bool(row["approach_raw_saturated"])
            ]
        values = finite_values(eligible_rows, field)
        # Event positions/amplitudes are summarized only when the event passed
        # the fixed detection threshold.
        if field.startswith("approach_snap") or field == "approach_barrier_force_nN":
            values = np.asarray(
                [
                    float(row[field])
                    for row in eligible_rows
                    if bool(row["approach_snap_detected"])
                    and np.isfinite(float(row[field]))
                ],
                dtype=np.float64,
            )
        if field.startswith("retract_detachment") or field in (
            "retract_pull_off_force_nN",
            "retract_pull_off_observed_floor_nN",
        ):
            values = np.asarray(
                [
                    float(row[field])
                    for row in eligible_rows
                    if bool(row["retract_snapoff_detected"])
                    and np.isfinite(float(row[field]))
                ],
                dtype=np.float64,
            )
        summary[f"median_{field}"] = (
            float(np.median(values)) if values.size else np.nan
        )
        summary[f"mad_{field}"] = base.robust_mad(values)
        summary[f"n_{field}"] = int(values.size)
    return summary


def within_map_trends(source: base.SourceData, rows: list[dict]) -> list[dict]:
    if source.source_type != "map":
        return []
    output: list[dict] = []
    order = np.asarray([float(row["acquisition_order"]) for row in rows])
    grid_row = np.asarray([float(row["map_row"]) for row in rows])
    grid_column = np.asarray([float(row["map_column"]) for row in rows])
    for field in TREND_METRICS:
        values = np.asarray([float(row[field]) for row in rows])
        if field in ("approach_snap_distance_nm", "approach_barrier_force_nN"):
            detected = np.asarray(
                [bool(row["approach_snap_detected"]) for row in rows]
            )
            values[~detected] = np.nan
        if field in (
            "retract_detachment_piezo_travel_nm",
            "retract_pull_off_force_nN",
        ):
            detected = np.asarray(
                [bool(row["retract_snapoff_detected"]) for row in rows]
            )
            values[~detected] = np.nan
        finite = np.isfinite(order) & np.isfinite(values)
        slope = intercept = correlation = total_change = np.nan
        if np.count_nonzero(finite) >= 8:
            slope, intercept, _, _ = base.robust_line(order[finite], values[finite])
            correlation = base._finite_correlation(order[finite], values[finite])
            total_change = slope * (
                float(np.max(order[finite])) - float(np.min(order[finite]))
            )
        plane = base._robust_spatial_plane(grid_column, grid_row, values)
        output.append(
            {
                "concentration_wt_percent": source.concentration_wt_percent,
                "primary_included": source.concentration_wt_percent
                in PRIMARY_CONCENTRATIONS,
                "source": str(source.path.relative_to(ROOT)),
                "timestamp": source.timestamp.isoformat(),
                "metric": field,
                "valid_pixels": int(np.count_nonzero(finite)),
                "median": float(np.nanmedian(values))
                if np.count_nonzero(finite)
                else np.nan,
                "mad": base.robust_mad(values),
                "robust_slope_per_acquisition": slope,
                "robust_change_across_map": total_change,
                "pearson_correlation_acquisition_order": correlation,
                "spatial_plane_column_slope_per_pixel": plane[1],
                "spatial_plane_row_slope_per_pixel": plane[2],
                "spatial_plane_r2": plane[3],
                "spatial_plane_inliers": plane[4],
            }
        )
    return output


def concentration_speed_summary(source_rows: list[dict]) -> list[dict]:
    output: list[dict] = []
    map_rows = [row for row in source_rows if row["source_type"] == "map"]
    for concentration in sorted(base.EPSILON_R):
        for speed in base.SPEED_LABELS_UM_PER_S:
            subset = [
                row
                for row in map_rows
                if int(row["concentration_wt_percent"]) == concentration
                and float(row["speed_label_um_per_s"]) == float(speed)
            ]
            if not subset:
                continue
            record: dict = {
                "concentration_wt_percent": concentration,
                "primary_included": concentration in PRIMARY_CONCENTRATIONS,
                "speed_label_um_per_s": float(speed),
                "map_count": len(subset),
            }
            for field in SUMMARY_METRICS:
                key = f"median_{field}"
                eligible = subset
                if field in (
                    "retract_pull_off_force_nN",
                    "retract_detachment_piezo_travel_nm",
                ):
                    eligible = [
                        row
                        for row in subset
                        if int(row[f"n_{field}"])
                        / max(1, int(row["curve_count"]))
                        >= MIN_MAP_EVENT_SUPPORT_FRACTION
                    ]
                values = finite_values(eligible, key)
                record[key] = float(np.median(values)) if values.size else np.nan
                record[f"map_mad_{field}"] = base.robust_mad(values)
                record[f"map_min_{field}"] = float(np.min(values)) if values.size else np.nan
                record[f"map_max_{field}"] = float(np.max(values)) if values.size else np.nan
            for fraction in (
                "approach_contact_valid_fraction",
                "approach_snap_detected_fraction",
                "retract_contact_valid_fraction",
                "approach_raw_saturation_fraction",
                "retract_raw_saturation_fraction",
                "retract_free_baseline_valid_fraction",
                "retract_snapoff_detected_fraction",
                "retract_pull_off_censored_fraction",
                "retract_detachment_right_censored_fraction",
            ):
                values = finite_values(subset, fraction)
                record[f"median_{fraction}"] = (
                    float(np.median(values)) if values.size else np.nan
                )
            output.append(record)
    return output


def metric_velocity_time_models(source_rows: list[dict]) -> list[dict]:
    """Descriptive source-level regressions; no pixel-level pseudoreplication."""

    output: list[dict] = []
    maps = [row for row in source_rows if row["source_type"] == "map"]
    for concentration in sorted(base.EPSILON_R):
        subset = [
            row
            for row in maps
            if int(row["concentration_wt_percent"]) == concentration
        ]
        subset.sort(key=lambda row: str(row["timestamp"]))
        for field in TREND_METRICS:
            key = f"median_{field}"
            values = np.asarray(
                [
                    float(row[key])
                    if field
                    not in (
                        "retract_pull_off_force_nN",
                        "retract_detachment_piezo_travel_nm",
                    )
                    or int(row[f"n_{field}"])
                    / max(1, int(row["curve_count"]))
                    >= MIN_MAP_EVENT_SUPPORT_FRACTION
                    else np.nan
                    for row in subset
                ]
            )
            speed = np.asarray([float(row["approach_speed_um_per_s"]) for row in subset])
            time_min = np.asarray(
                [
                    (
                        base.datetime.fromisoformat(str(row["timestamp"]))
                        - base.datetime.fromisoformat(str(subset[0]["timestamp"]))
                    ).total_seconds()
                    / 60.0
                    for row in subset
                ]
            )
            finite = np.isfinite(values) & np.isfinite(speed) & np.isfinite(time_min)
            x_speed = speed[finite]
            x_time = time_min[finite]
            y = values[finite]
            speed_corr = base._finite_correlation(x_speed, y)
            time_corr = base._finite_correlation(x_time, y)
            speed_coef = time_coef = model_r2 = np.nan
            if y.size >= 5 and np.std(x_speed) > 0.0 and np.std(x_time) > 0.0:
                speed_z = (x_speed - np.mean(x_speed)) / np.std(x_speed)
                time_z = (x_time - np.mean(x_time)) / np.std(x_time)
                design = np.column_stack([np.ones(y.size), speed_z, time_z])
                coefficients, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
                if rank == 3:
                    prediction = design @ coefficients
                    ss_res = float(np.sum((y - prediction) ** 2))
                    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
                    speed_coef = float(coefficients[1])
                    time_coef = float(coefficients[2])
                    model_r2 = np.nan if ss_tot <= 0.0 else 1.0 - ss_res / ss_tot
            output.append(
                {
                    "concentration_wt_percent": concentration,
                    "primary_included": concentration in PRIMARY_CONCENTRATIONS,
                    "metric": field,
                    "map_count": int(y.size),
                    "pearson_with_speed": speed_corr,
                    "pearson_with_time": time_corr,
                    "two_predictor_speed_coefficient_per_1sd_speed": speed_coef,
                    "two_predictor_time_coefficient_per_1sd_time": time_coef,
                    "two_predictor_r2": model_r2,
                    "interpretation": "descriptive_map_level_not_inferential",
                }
            )
    return output


METRIC_LABELS = {
    "median_approach_contact_invOLS_nm_per_V": "Approach contact InvOLS (nm/V)",
    "median_approach_contact_response_ratio": "Contact response ratio (global/local)",
    "median_approach_snap_distance_nm": "Snap-in distance (nm)",
    "median_approach_terminal_load_nN": "Approach terminal load (nN)",
    "median_retract_pull_off_force_nN": "Pull-off force (nN)",
    "median_retract_detachment_piezo_travel_nm": "Detachment piezo travel (nm)",
}


def supported_source_metric(row: dict, metric: str) -> float:
    value = float(row[metric])
    field = metric.removeprefix("median_")
    if field in (
        "retract_pull_off_force_nN",
        "retract_detachment_piezo_travel_nm",
    ):
        support = int(row[f"n_{field}"]) / max(1, int(row["curve_count"]))
        if row["source_type"] == "map" and support < MIN_MAP_EVENT_SUPPORT_FRACTION:
            return np.nan
    return value


def plot_source_time_series(source_rows: list[dict], path: Path) -> None:
    entries = list(source_rows)
    maps = [row for row in entries if row["source_type"] == "map"]
    start = min(base.datetime.fromisoformat(str(row["timestamp"])) for row in entries)
    metrics = list(METRIC_LABELS)
    colors = {0: "#1f77b4", 10: "#8c8c8c", 20: "#2ca02c", 30: "#d62728", 40: "#9467bd"}
    markers = {1.0: "o", 2.0: "s", 4.0: "^"}
    fig, axes = plt.subplots(3, 2, figsize=(13, 12), constrained_layout=True)
    for ax, metric in zip(axes.ravel(), metrics):
        for row in entries:
            concentration = int(row["concentration_wt_percent"])
            speed = float(row["speed_label_um_per_s"])
            elapsed = (
                base.datetime.fromisoformat(str(row["timestamp"])) - start
            ).total_seconds() / 60.0
            ax.scatter(
                elapsed,
                supported_source_metric(row, metric),
                facecolor=colors[concentration]
                if row["source_type"] == "map"
                else "none",
                edgecolor=colors[concentration]
                if concentration != 10
                else "black",
                marker=markers[speed],
                s=58 if row["source_type"] == "map" else 42,
                linewidth=1.0,
                zorder=3,
            )
        for concentration in sorted(base.EPSILON_R):
            subset = [row for row in maps if int(row["concentration_wt_percent"]) == concentration]
            subset.sort(key=lambda row: str(row["timestamp"]))
            x = [
                (base.datetime.fromisoformat(str(row["timestamp"])) - start).total_seconds() / 60.0
                for row in subset
            ]
            y = [supported_source_metric(row, metric) for row in subset]
            ax.plot(x, y, color=colors[concentration], alpha=0.45, linewidth=1.2)
        ax.set_xlabel("Elapsed experiment time (min)")
        ax.set_ylabel(METRIC_LABELS[metric])
        ax.grid(alpha=0.2)
    handles = [
        plt.Line2D([], [], color=colors[c], marker="o", linestyle="", label=f"{c} wt%")
        for c in sorted(base.EPSILON_R)
    ] + [
        plt.Line2D([], [], color="black", marker=markers[s], linestyle="", label=f"{s:g} um/s")
        for s in (1.0, 2.0, 4.0)
    ] + [
        plt.Line2D([], [], color="black", marker="o", markerfacecolor="black", linestyle="", label="map median"),
        plt.Line2D([], [], color="black", marker="o", markerfacecolor="none", linestyle="", label="independent force curve"),
    ]
    fig.legend(handles=handles, loc="outside upper center", ncol=6)
    fig.suptitle("Contact and adhesion metrics versus acquisition time")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_speed_summary(source_rows: list[dict], path: Path) -> None:
    maps = [row for row in source_rows if row["source_type"] == "map"]
    metrics = list(METRIC_LABELS)
    colors = {0: "#1f77b4", 10: "#8c8c8c", 20: "#2ca02c", 30: "#d62728", 40: "#9467bd"}
    fig, axes = plt.subplots(3, 2, figsize=(12, 12), constrained_layout=True)
    for ax, metric in zip(axes.ravel(), metrics):
        for concentration in sorted(base.EPSILON_R):
            subset = [row for row in maps if int(row["concentration_wt_percent"]) == concentration]
            by_speed: dict[float, list[float]] = {}
            for row in subset:
                value = supported_source_metric(row, metric)
                if np.isfinite(value):
                    by_speed.setdefault(
                        float(row["speed_label_um_per_s"]), []
                    ).append(value)
            speeds = sorted(by_speed)
            medians = [float(np.nanmedian(by_speed[speed])) for speed in speeds]
            ax.plot(
                speeds,
                medians,
                marker="o",
                color=colors[concentration],
                linestyle="--" if concentration == 10 else "-",
                linewidth=1.5,
                label=f"{concentration} wt%",
            )
            for speed in speeds:
                values = np.asarray(by_speed[speed], dtype=float)
                jitter = np.linspace(-0.05, 0.05, values.size) if values.size > 1 else np.zeros(1)
                ax.scatter(speed + jitter, values, color=colors[concentration], s=32, alpha=0.8)
        ax.set_xticks([1, 2, 4])
        ax.set_xlabel("Nominal approach speed (um/s)")
        ax.set_ylabel(METRIC_LABELS[metric])
        ax.grid(alpha=0.2)
    axes[0, 0].legend(ncol=2, fontsize=9)
    fig.suptitle("Velocity comparison uses complete maps as experimental units")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_heatmap_montage(
    curve_rows: list[dict],
    source_rows: list[dict],
    metric: str,
    detected_field: str | None,
    label: str,
    path: Path,
) -> None:
    maps = [row for row in source_rows if row["source_type"] == "map"]
    map_by_concentration = {
        concentration: [
            row for row in maps if int(row["concentration_wt_percent"]) == concentration
        ]
        for concentration in sorted(base.EPSILON_R)
    }
    values = np.asarray(
        [
            float(row[metric])
            for row in curve_rows
            if row["source_type"] == "map"
            and (detected_field is None or bool(row[detected_field]))
            and np.isfinite(float(row[metric]))
        ]
    )
    if not values.size:
        return
    low, high = np.nanpercentile(values, [2.0, 98.0])
    if low == high:
        low -= 1.0
        high += 1.0
    fig, axes = plt.subplots(5, 6, figsize=(16, 13), constrained_layout=True)
    image = None
    for row_index, concentration in enumerate(sorted(base.EPSILON_R)):
        entries = sorted(map_by_concentration[concentration], key=lambda row: str(row["timestamp"]))
        for column_index in range(6):
            ax = axes[row_index, column_index]
            if column_index >= len(entries):
                ax.axis("off")
                continue
            source = str(entries[column_index]["source"])
            subset = [row for row in curve_rows if row["source"] == source]
            grid = np.full((16, 16), np.nan)
            for item in subset:
                if detected_field is not None and not bool(item[detected_field]):
                    continue
                physical_row = item["map_row"]
                physical_column = item["map_column"]
                value = float(item[metric])
                if physical_row != "" and physical_column != "" and np.isfinite(value):
                    grid[int(physical_row), int(physical_column)] = value
            image = ax.imshow(grid, origin="upper", vmin=low, vmax=high, cmap="viridis", interpolation="nearest")
            ax.set_title(
                f"{str(entries[column_index]['timestamp'])[11:19]}\n{float(entries[column_index]['speed_label_um_per_s']):g} um/s",
                fontsize=8,
            )
            ax.set_xticks([])
            ax.set_yticks([])
            if column_index == 0:
                ax.set_ylabel(f"{concentration} wt%", fontsize=10)
    if image is not None:
        fig.colorbar(image, ax=axes, label=label, shrink=0.72)
    fig.suptitle(f"Per-pixel {label}; common 2--98 percentile color scale")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def self_test() -> dict:
    count = 1000
    force = np.zeros(count)
    force[700:] = -450.0
    force += 0.8 * np.sin(np.arange(count) * 0.17)
    mask = np.ones(count - STEP_LAG, dtype=bool)
    event = detect_negative_step(force, mask, 200, 1.0, 8.0, 50.0)
    if not event.detected or abs(event.left_index - 698) > 4:
        raise AssertionError("synthetic negative step was not recovered")
    smooth = 0.01 * np.arange(count) + 0.8 * np.sin(np.arange(count) * 0.17)
    no_event = detect_negative_step(smooth, mask, 200, 1.0, 8.0, 50.0)
    if no_event.detected:
        raise AssertionError("smooth synthetic trace produced a false event")
    if base.map_pixel_from_index(16, 16, True) != (1, 15):
        raise AssertionError("serpentine raster mapping failed")
    return {
        "synthetic_step_detected": event.detected,
        "synthetic_step_index_error_points": abs(event.left_index - 698),
        "smooth_trace_rejected": not no_event.detected,
        "serpentine_mapping_check": True,
    }


def build_report(
    sources: list[base.SourceData],
    source_rows: list[dict],
    speed_rows: list[dict],
    model_rows: list[dict],
    trend_rows: list[dict],
    self_checks: dict,
) -> str:
    del speed_rows, model_rows
    maps = [row for row in source_rows if row["source_type"] == "map"]
    maps.sort(key=lambda row: str(row["timestamp"]))
    force_sources = [
        row for row in source_rows if row["source_type"] == "force"
    ]

    def concentration_maps(concentration: int) -> list[dict]:
        return [
            row
            for row in maps
            if int(row["concentration_wt_percent"]) == concentration
        ]

    maps_30 = concentration_maps(30)
    same_speed_30 = [
        row for row in maps_30 if float(row["speed_label_um_per_s"]) == 2.0
    ]
    contact_trends = [
        row
        for row in trend_rows
        if row["metric"] == "approach_contact_invOLS_nm_per_V"
        and bool(row["primary_included"])
        and np.isfinite(float(row["robust_change_across_map"]))
    ]
    strongest_contact_trend = max(
        contact_trends,
        key=lambda row: abs(float(row["robust_change_across_map"])),
        default=None,
    )
    contact_height_trends = [
        row
        for row in trend_rows
        if row["metric"] == "approach_contact_height_um"
        and bool(row["primary_included"])
        and np.isfinite(float(row["robust_change_across_map"]))
    ]
    strongest_contact_height_trend = max(
        contact_height_trends,
        key=lambda row: abs(float(row["robust_change_across_map"])),
        default=None,
    )
    clean_zero_speed2 = sorted(
        [
            row
            for row in force_sources
            if int(row["concentration_wt_percent"]) == 0
            and float(row["speed_label_um_per_s"]) == 2.0
            and np.isfinite(float(row["median_retract_pull_off_force_nN"]))
        ],
        key=lambda row: str(row["timestamp"]),
    )
    lines = [
        "# Velocity-correlated AFM systematics",
        "",
        "## Scope and physical definitions",
        "",
        "This analysis re-decodes raw `vDeflection` and `measuredHeight` for paired approach/retract branches. Embedded JPK sensitivity and force are not used. One raw hard-contact sensitivity consensus is recomputed per glycerol concentration and force is rebuilt with the calibrated cantilever-1 spring constant.",
        "",
        "- **Approach snap-in distance**: contact-aligned `h + delta - h_contact` immediately before a statistically resolved downward force step. This is a tip-surface distance estimate.",
        "- **Retract detachment position**: scanner/piezo travel `h_detach - h_contact` immediately before snap-off. It is not called a true gap because the probe remains adhered before release.",
        "- **Pull-off force**: minimum reconstructed force on the adhered retract branch between snap-off and hard contact.",
        "- **Contact response ratio**: common concentration sensitivity divided by the local contact InvOLS. A rigid, unchanged optical/contact response is near one; systematic changes flag calibration/contact mechanics bias.",
        "- Velocity interpretation uses whole maps as experimental units. Pixel distributions characterize within-map heterogeneity only. The 10 wt% group is retained as anomaly QC and excluded from primary interpretation.",
        "",
        "## Main findings",
        "",
        "1. **Snap-in changes, but there is no concentration-independent monotonic velocity law.** The 4 um/s map has the smallest median detected snap-in distance in 0, 30 and 40 wt%, which is qualitatively compatible with hydrodynamic squeeze pressure delaying an instability. However, 20 wt% has nearly the same 2 and 4 um/s value and its 1 um/s value is largest; detection fractions also change strongly between maps. The result is therefore a possible dynamic contribution superposed on larger history/calibration effects, not an independently fitted velocity effect.",
        "",
    ]
    if len(same_speed_30) >= 3:
        first_30 = same_speed_30[0]
        last_30 = same_speed_30[-1]
        lines.append(
            "2. **A fixed-speed control directly shows time/history drift.** In the three consecutive 30 wt%, 2 um/s maps, median snap-in distance changes from "
            f"{float(first_30['median_approach_snap_distance_nm']):.3g} to "
            f"{float(last_30['median_approach_snap_distance_nm']):.3g} nm, while local contact InvOLS changes from "
            f"{float(first_30['median_approach_contact_invOLS_nm_per_V']):.3g} to "
            f"{float(last_30['median_approach_contact_invOLS_nm_per_V']):.3g} nm/V. The terminal load simultaneously changes from "
            f"{float(first_30['median_approach_terminal_load_nN']):.3g} to "
            f"{float(last_30['median_approach_terminal_load_nN']):.3g} nN. These changes cannot be caused by velocity because velocity is unchanged."
        )
        lines.append("")
    if strongest_contact_trend is not None:
        lines.append(
            "3. **The drift also occurs inside a single map.** The largest robust contact-InvOLS change across one primary map is "
            f"{float(strongest_contact_trend['robust_change_across_map']):.3g} nm/V "
            f"({strongest_contact_trend['concentration_wt_percent']} wt%, "
            f"{str(strongest_contact_trend['timestamp'])[11:19]}), with spatial/acquisition plane R2 "
            f"{float(strongest_contact_trend['spatial_plane_r2']):.3f}. A fixed-speed map therefore contains a systematic optical/contact-response gradient."
        )
        lines.append("")
    lines += [
        "4. **Retract adhesion is heavily censored by the raw encoder floor.** This is not a subtle statistical effect: many raw signed-int32 traces equal the lower encoding bound. Exact pull-off values are therefore available only for unclipped adhered branches. In particular, most 30 wt% traces do not contain a verified unsaturated free-cantilever prefix, so their detachment is not observed within the usable retract record; reporting near-zero pull-off for them would be wrong.",
        "",
    ]
    if len(clean_zero_speed2) >= 2:
        early = clean_zero_speed2[0]
        late = clean_zero_speed2[-1]
        lines += [
            "5. **Independent force curves also expose history dependence.** At 0 wt% and the same approximately 2 um/s speed, exact pull-off changes from "
            f"{float(early['median_retract_pull_off_force_nN']):.3g} nN at {str(early['timestamp'])[11:19]} to "
            f"{float(late['median_retract_pull_off_force_nN']):.3g} nN at {str(late['timestamp'])[11:19]}; detachment travel changes from "
            f"{float(early['median_retract_detachment_piezo_travel_nm']):.3g} to "
            f"{float(late['median_retract_detachment_piezo_travel_nm']):.3g} nm. This same-speed change is point heterogeneity, surface/probe conditioning, or time drift—not hydrodynamic velocity dependence.",
            "",
        ]
    lines += [
        "6. **The 10 wt% block is a different measurement regime.** Its approach hard-contact validity is low and terminal loads are about 1 nN rather than the roughly 7 or 18 nN regimes in the primary data. Retaining it as QC supports its exclusion from the primary surface-force fit.",
        "",
    ]
    if strongest_contact_height_trend is not None:
        lines += [
            "7. **Separation zero has a large map plane.** The largest robust hard-contact-height change across a primary map is "
            f"{float(strongest_contact_height_trend['robust_change_across_map']) * 1e3:.3g} nm "
            f"({strongest_contact_height_trend['concentration_wt_percent']} wt%, "
            f"{str(strongest_contact_height_trend['timestamp'])[11:19]}), with plane R2 "
            f"{float(strongest_contact_height_trend['spatial_plane_r2']):.3f}. The repeatable column component is compatible with sample tilt/topography, while changing row components can also contain piezo creep or time drift. Per-pixel contact alignment removes the first-order offset; a single common contact zero would not.",
            "",
        ]
    maps_20 = concentration_maps(20)
    if len(maps_20) >= 3:
        slopes_20 = [
            float(row["median_approach_far_slope_pN_per_100nm"])
            for row in maps_20
        ]
        lines += [
            "8. **Far-field photodiode drift is also map dependent and nonmonotonic.** At 20 wt% its chronological map medians are "
            f"{slopes_20[0]:.3g}, {slopes_20[1]:.3g} and {slopes_20[2]:.3g} pN per 100 nm for the 2, 1 and 4 um/s maps. The robust line is subtracted independently from every approach trace, so this is retained as instrumental QC rather than mistaken for a surface force.",
            "",
        ]
    lines += [
        "## Map-level results",
        "",
        "| wt% | time | speed | contact valid | local InvOLS (nm/V) | snap detected | snap D (nm) | terminal load (nN) | retract clipped | free baseline | snap-off detected | exact pull-off (nN) | exact n / total | detach piezo (nm) |",
        "|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in maps:
        def fmt(key: str, digits: int = 3) -> str:
            value = float(row[key])
            return "NA" if not np.isfinite(value) else f"{value:.{digits}g}"

        lines.append(
            f"| {row['concentration_wt_percent']} | {str(row['timestamp'])[11:19]} | "
            f"{float(row['approach_speed_um_per_s']):.3g} | "
            f"{float(row['approach_contact_valid_fraction']):.1%} | "
            f"{fmt('median_approach_contact_invOLS_nm_per_V')} | "
            f"{float(row['approach_snap_detected_fraction']):.1%} | "
            f"{fmt('median_approach_snap_distance_nm')} | "
            f"{fmt('median_approach_terminal_load_nN')} | "
            f"{float(row['retract_raw_saturation_fraction']):.1%} | "
            f"{float(row['retract_free_baseline_valid_fraction']):.1%} | "
            f"{float(row['retract_snapoff_detected_fraction']):.1%} | "
            f"{fmt('median_retract_pull_off_force_nN')} | "
            f"{int(row['n_retract_pull_off_force_nN'])}/{int(row['curve_count'])} | "
            f"{fmt('median_retract_detachment_piezo_travel_nm')} |"
        )
    lines += [
        "",
        "The table is descriptive. In conditions with one map at each speed, speed is confounded with acquisition time and surface history. The six 30 wt% maps are unbalanced (three at 2 um/s, two at 1 um/s and one at 4 um/s), so they provide fixed-speed drift controls but not balanced velocity replication.",
        "",
        "`Retract clipped` means that raw signed-int32 `vDeflection` reached the JPK encoder floor. Such pull-off amplitudes are censored and are not included as exact force values. `Free baseline` requires at least 80 leading unsaturated samples after normalizing retract to far-to-contact orientation. If the physical retract ends at the encoder floor, detachment was not verifiably observed within the available travel; its position and pull-off force are right/left censored rather than zero.",
        f"For velocity summaries/regressions, retract event metrics are used only when at least {MIN_MAP_EVENT_SUPPORT_FRACTION:.0%} of a map has an uncensored/detected value. The per-source table retains smaller subsets visibly, with their exact numerator and total, for QC only.",
        "",
        "## Physical interpretation",
        "",
        "Quasi-statically, snap-in is a cantilever instability when the attractive interaction-force gradient reaches the restoring stiffness. A finite approach rate can shift that event because hydrodynamic pressure and cantilever dynamics alter the instantaneous force balance; dynamic drive-velocity effects are treated explicitly by [Bowen and Cheneler, Langmuir 2012](https://doi.org/10.1021/la304009c), and hydrodynamic pressure delaying jump-to-contact in liquid is analysed by [Dai, Journal of Fluid Mechanics 2025](https://doi.org/10.1017/jfm.2025.10290). Thus a velocity-dependent snap position is physically possible, but only after contact-zero, optical response, load and history are controlled.",
        "",
        "Approach/retract comparison can isolate a hydrodynamic component when the non-hydrodynamic surface force is reproducible between branches and speeds; an example with silica sphere/plate is [McNamee et al., Colloids and Surfaces A 2007](https://doi.org/10.1016/j.colsurfa.2007.01.047). The present data violate that prerequisite through fixed-speed drift, load-regime changes and retract clipping. Large or gradual silica-silica adhesion in water has also been reported and attributed provisionally to bubbles/cavities by [Troncoso et al., Journal of Colloid and Interface Science 2014](https://doi.org/10.1016/j.jcis.2014.03.020); the current curves alone cannot identify that mechanism.",
        "",
        "## Event-detection and numerical semantics",
        "",
        f"Approach events require a two-sample median-filtered drop larger than max({APPROACH_STEP_SIGMA:g} x far-field step MAD, 5 x far-field force residual MAD, {APPROACH_ABSOLUTE_STEP_PN:g} pN). Retract snap-off uses max({RETRACT_STEP_SIGMA:g} x step MAD, 5 x residual MAD, {RETRACT_ABSOLUTE_STEP_PN:g} pN). Candidate values are retained per curve, but detected-event summaries include only threshold-passing events.",
        "",
        f"Decoded {len(sources)} raw sources and {sum(len(source.curves) for source in sources)} approach curves. Synthetic step recovery: `{self_checks['synthetic_step_detected']}`; smooth no-event rejection: `{self_checks['smooth_trace_rejected']}`; serpentine map-coordinate check: `{self_checks['serpentine_mapping_check']}`.",
        "",
        "## Files",
        "",
        "- `curve_event_metrics.csv`: every paired branch, event candidates/detections, contact response, loads, drift/noise and explicit coordinates.",
        "- `source_event_summary.csv`: robust per-source medians/MADs; map is the velocity-comparison unit.",
        "- `within_map_metric_trends.csv`: acquisition-order correlations and spatial-plane gradients.",
        "- `concentration_speed_summary.csv`: speed summaries of map medians, not pooled pixels.",
        "- `velocity_time_models.csv`: descriptive map-level speed/time correlations and, where possible, a two-predictor model.",
        "- `input_manifest.csv`, `provenance.json`, figures and SHA-256 artifact manifest.",
    ]
    return "\n".join(lines) + "\n"


def run(results_dir: Path) -> None:
    checks = self_test()
    sources: list[base.SourceData] = []
    for concentration in sorted(base.EPSILON_R):
        directory = DATA_ROOT / str(concentration)
        paths = sorted(directory.glob("*.jpk-force")) + sorted(
            directory.glob("*.jpk-force-map")
        )
        for path in paths:
            sources.append(base.load_source(path, concentration))
    sources.sort(key=lambda source: source.timestamp)
    base.calibrate_sensitivity(sources)

    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    curve_rows: list[dict] = []
    source_rows: list[dict] = []
    trend_rows: list[dict] = []
    input_rows: list[dict] = []
    for source in sources:
        raw_approaches, approach_raw_skipped = load_branch(source.path, "extend")
        retracts, retract_skipped = load_branch(source.path, "retract")
        pairs = pair_curves(source.curves, retracts)
        raw_approach_by_key = {
            _curve_key(curve, order): curve
            for order, curve in enumerate(raw_approaches)
        }
        current_rows: list[dict] = []
        for order, (approach, retract) in enumerate(pairs):
            raw_approach = raw_approach_by_key.get(_curve_key(approach, order))
            map_row: int | str = ""
            map_column: int | str = ""
            if (
                source.source_type == "map"
                and approach.point_index is not None
                and source.map_grid_i is not None
            ):
                map_row, map_column = base.map_pixel_from_index(
                    approach.point_index,
                    source.map_grid_i,
                    source.map_back_and_forth,
                )
            row: dict = {
                "concentration_wt_percent": source.concentration_wt_percent,
                "primary_included": source.concentration_wt_percent
                in PRIMARY_CONCENTRATIONS,
                "source": str(source.path.relative_to(ROOT)),
                "source_type": source.source_type,
                "timestamp": source.timestamp.isoformat(),
                "point_index": "" if approach.point_index is None else approach.point_index,
                "acquisition_order": order,
                "map_row": map_row,
                "map_column": map_column,
                "approach_points": approach.measured_height_m.size,
                "retract_points": 0 if retract is None else retract.measured_height_m.size,
                "approach_duration_s": approach.duration_s,
                "retract_duration_s": np.nan if retract is None else retract.duration_s,
                "approach_speed_um_per_s": estimate_speed_um_per_s(
                    approach.measured_height_m, approach.duration_s
                ),
                "retract_speed_um_per_s": np.nan
                if retract is None
                else estimate_speed_um_per_s(
                    retract.measured_height_m, retract.duration_s
                ),
                "sensitivity_used_nm_per_V": source.sensitivity_used_m_per_V * 1e9,
                "spring_constant_used_N_per_m": base.SPRING_CONSTANT_N_PER_M,
                "approach_raw_saturated": False
                if raw_approach is None
                else bool(np.any(raw_approach.raw_saturation_mask)),
                "approach_raw_low_saturation_fraction": np.nan
                if raw_approach is None
                else raw_approach.raw_low_saturation_fraction,
                "approach_raw_high_saturation_fraction": np.nan
                if raw_approach is None
                else raw_approach.raw_high_saturation_fraction,
                "approach_raw_min": np.nan
                if raw_approach is None
                else raw_approach.raw_min,
                "approach_raw_max": np.nan
                if raw_approach is None
                else raw_approach.raw_max,
                "retract_raw_saturated": False
                if retract is None
                else bool(np.any(retract.raw_saturation_mask)),
                "retract_raw_low_saturation_fraction": np.nan
                if retract is None
                else retract.raw_low_saturation_fraction,
                "retract_raw_high_saturation_fraction": np.nan
                if retract is None
                else retract.raw_high_saturation_fraction,
                "retract_far_raw_saturation_fraction": np.nan
                if retract is None
                else retract.far_raw_saturation_fraction,
                "retract_first_sample_saturated": False
                if retract is None
                else retract.first_sample_saturated,
                "retract_leading_unsaturated_points": 0
                if retract is None
                else retract.leading_unsaturated_points,
                "retract_raw_min": np.nan if retract is None else retract.raw_min,
                "retract_raw_max": np.nan if retract is None else retract.raw_max,
            }
            try:
                row.update(
                    analyze_pair(
                        approach,
                        retract,
                        source.sensitivity_used_m_per_V,
                    )
                )
            except (ValueError, np.linalg.LinAlgError, FloatingPointError):
                row.update(_nan_record())
            current_rows.append(row)
            curve_rows.append(row)
        contact_heights = finite_values(
            current_rows, "approach_contact_height_um"
        )
        contact_height_center_um = (
            float(np.median(contact_heights))
            if contact_heights.size
            else np.nan
        )
        for row in current_rows:
            height_um = float(row["approach_contact_height_um"])
            row["approach_contact_height_relative_source_median_nm"] = (
                (height_um - contact_height_center_um) * 1e3
                if np.isfinite(height_um) and np.isfinite(contact_height_center_um)
                else np.nan
            )
        source_rows.append(source_summary(source, current_rows))
        trend_rows.extend(within_map_trends(source, current_rows))
        input_rows.append(
            {
                "source": str(source.path.relative_to(ROOT)),
                "sha256": source.sha256,
                "size_bytes": source.path.stat().st_size,
                "zip_crc_pass": True,
                "approach_curves": len(source.curves),
                "raw_approach_curves": len(raw_approaches),
                "raw_approach_skipped": approach_raw_skipped,
                "retract_curves": len(retracts),
                "retract_skipped": retract_skipped,
            }
        )

    speed_rows = concentration_speed_summary(source_rows)
    model_rows = metric_velocity_time_models(source_rows)
    write_csv(results_dir / "curve_event_metrics.csv", curve_rows)
    write_csv(results_dir / "source_event_summary.csv", source_rows)
    write_csv(results_dir / "within_map_metric_trends.csv", trend_rows)
    write_csv(results_dir / "concentration_speed_summary.csv", speed_rows)
    write_csv(results_dir / "velocity_time_models.csv", model_rows)
    write_csv(results_dir / "input_manifest.csv", input_rows)

    plot_source_time_series(source_rows, figures_dir / "source_metrics_vs_time.png")
    plot_speed_summary(source_rows, figures_dir / "map_metrics_vs_speed.png")
    plot_heatmap_montage(
        curve_rows,
        source_rows,
        "approach_contact_invOLS_nm_per_V",
        None,
        "approach contact InvOLS (nm/V)",
        figures_dir / "contact_invOLS_maps.png",
    )
    plot_heatmap_montage(
        curve_rows,
        source_rows,
        "approach_contact_height_relative_source_median_nm",
        None,
        "contact height relative to source median (nm)",
        figures_dir / "relative_contact_height_maps.png",
    )
    plot_heatmap_montage(
        curve_rows,
        source_rows,
        "approach_snap_distance_nm",
        "approach_snap_detected",
        "approach snap-in distance (nm)",
        figures_dir / "snap_in_distance_maps.png",
    )
    plot_heatmap_montage(
        curve_rows,
        source_rows,
        "retract_detachment_piezo_travel_nm",
        "retract_snapoff_detected",
        "retract detachment piezo travel (nm)",
        figures_dir / "detachment_travel_maps.png",
    )
    plot_heatmap_montage(
        curve_rows,
        source_rows,
        "retract_raw_low_saturation_fraction",
        None,
        "retract encoder-floor sample fraction",
        figures_dir / "retract_encoder_floor_maps.png",
    )
    plot_heatmap_montage(
        curve_rows,
        source_rows,
        "retract_pull_off_force_nN",
        "retract_snapoff_detected",
        "retract pull-off force (nN)",
        figures_dir / "pull_off_force_maps.png",
    )

    provenance = {
        "analysis_script": str(Path(__file__).relative_to(ROOT)),
        "analysis_script_sha256": sha256_file(Path(__file__)),
        "raw_decoder_script": str(Path(base.__file__).relative_to(ROOT)),
        "raw_decoder_script_sha256": sha256_file(Path(base.__file__)),
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "temperature_C": base.TEMPERATURE_C,
        "spring_constant_N_per_m": base.SPRING_CONSTANT_N_PER_M,
        "jpk_embedded_sensitivity_used": False,
        "jpk_embedded_force_used": False,
        "sensitivity_method": "equal-source-weighted concentration consensus from raw approach hard-contact fits",
        "primary_concentrations_wt_percent": list(PRIMARY_CONCENTRATIONS),
        "qc_excluded_concentration_wt_percent": QC_EXCLUDED_CONCENTRATION,
        "experimental_unit_for_velocity": "one complete map",
        "pixel_role": "within-map heterogeneity and spatial/acquisition-order diagnostics only",
        "minimum_map_event_support_fraction_for_velocity_summary": MIN_MAP_EVENT_SUPPORT_FRACTION,
        "event_detector": {
            "orientation": "both branches normalized far-to-contact",
            "median_filter_points": SMOOTH_POINTS,
            "difference_lag_points": STEP_LAG,
            "approach_sigma_multiplier": APPROACH_STEP_SIGMA,
            "retract_sigma_multiplier": RETRACT_STEP_SIGMA,
            "approach_absolute_floor_pN": APPROACH_ABSOLUTE_STEP_PN,
            "retract_absolute_floor_pN": RETRACT_ABSOLUTE_STEP_PN,
            "far_residual_multiplier": 5.0,
        },
        "raw_encoder_floor_qc": {
            "integer_margin_counts": "max(1024, 1e-6 of integer range)",
            "observed_signed_int32_floor": -2147483600,
            "pull_off_policy": "exact only when adhered branch contains no floor-clipped samples",
            "free_retract_baseline_policy": "at least 80 leading unsaturated far-to-contact samples",
        },
        "randomness": "none",
        "self_checks": checks,
    }
    (results_dir / "provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    report = build_report(
        sources, source_rows, speed_rows, model_rows, trend_rows, checks
    )
    (results_dir / "REPORT.md").write_text(report, encoding="utf-8")
    artifacts = sorted(
        path
        for path in results_dir.rglob("*")
        if path.is_file() and path.name != "artifact_manifest.sha256"
    )
    manifest = "\n".join(
        f"{sha256_file(path)}  {path.relative_to(results_dir)}" for path in artifacts
    )
    (results_dir / "artifact_manifest.sha256").write_text(
        manifest + "\n", encoding="utf-8"
    )


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

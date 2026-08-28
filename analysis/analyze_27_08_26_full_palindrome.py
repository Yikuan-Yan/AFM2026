#!/usr/bin/env python3
"""Full raw analysis of the overlapping 27-08-26 palindrome experiments.

Experimental definitions supplied by the user:

* map1-map3: one complete three-block test, with a same-composition liquid
  refresh during the test.  The refresh boundary is represented between map2
  and map3.
* map3-map5: one complete no-refresh test.  Map3 is the shared bridge block.
* map5 is missing the final 0.1 um/s map.  Missingness is explicit; the primary
  velocity contrast uses the complete 0.05-to-0.2 um/s pair in every block.

Every force is rebuilt from raw vDeflection/measuredHeight, an all-map water
InvOLS consensus, and the independently air-calibrated D4 spring constant.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import re
import sys
import zipfile
from datetime import timedelta
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.stats


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import analyze_27_08_26_palindrome_pilot as pilot  # noqa: E402
import fit_glycerol_surface_forces as base  # noqa: E402


DATA_ROOT = ROOT / "27-08-26"
RESULTS = ROOT / "analysis" / "palindrome_27_08_26_full_results"
FIGURES = RESULTS / "figures"

TEMPERATURE_C = 25.6
SPEEDS = (0.05, 0.1, 0.2)
TARGET_DISTANCES_NM = (20.0, 50.0, 100.0, 200.0)
COLORS = {0.05: "#2a9d8f", 0.1: "#e9c46a", 0.2: "#e76f51"}
BLOCK_COLORS = {1: "#577590", 2: "#43aa8b", 3: "#f9c74f", 4: "#f9844a", 5: "#9b5de5"}
TESTS = {
    "map1_3_refresh_affected": {
        "blocks": (1, 2, 3),
        "primary": False,
        "description": "complete rotated palindrome; liquid refreshed between map2 and map3",
    },
    "map3_5_no_refresh": {
        "blocks": (3, 4, 5),
        "primary": True,
        "description": "primary complete rotated palindrome without liquid refresh; map5 late 0.1 map missing",
    },
}


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite_summary(values: np.ndarray) -> tuple[float, float, float, int]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not values.size:
        return float("nan"), float("nan"), float("nan"), 0
    return (
        float(np.median(values)),
        float(np.quantile(values, 0.25)),
        float(np.quantile(values, 0.75)),
        int(values.size),
    )


def finite_column_quantiles(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Column summaries without emitting warnings for unsupported distance bins."""

    matrix = np.asarray(matrix, dtype=float)
    median = np.full(matrix.shape[1], np.nan)
    q25 = np.full(matrix.shape[1], np.nan)
    q75 = np.full(matrix.shape[1], np.nan)
    supported = np.any(np.isfinite(matrix), axis=0)
    if np.any(supported):
        median[supported] = np.nanmedian(matrix[:, supported], axis=0)
        q25[supported] = np.nanquantile(matrix[:, supported], 0.25, axis=0)
        q75[supported] = np.nanquantile(matrix[:, supported], 0.75, axis=0)
    return median, q25, q75


def source_key(source: base.SourceData) -> str:
    return str(source.path.relative_to(ROOT))


def add_map_timing_metadata(
    sources: list[base.SourceData], map_rows: list[dict]
) -> None:
    """Add protocol-based map midpoint estimates for time-confounding models.

    JPK records a reliable map start time but repeats that value in the end-time
    field.  The midpoint below therefore uses the sum of the stored approach
    durations plus the retract duration declared in the map header.  XY motion
    and baseline-adjust overhead are not represented, so the field is explicitly
    labelled as an estimate rather than an acquisition timestamp.
    """

    rows = {row["source"]: row for row in map_rows}
    retract_duration_pattern = re.compile(
        r"force-scan-map\.settings\.force-settings\.segment\.(\d+)\.style=retract"
    )
    for source in sources:
        with zipfile.ZipFile(source.path) as archive:
            header = archive.read("header.properties").decode("iso-8859-1")
        match = retract_duration_pattern.search(header)
        if match is None:
            raise RuntimeError(f"no retract segment in {source.path}")
        segment = match.group(1)
        duration_match = re.search(
            rf"force-scan-map\.settings\.force-settings\.segment\.{segment}\.duration=([^\r\n]+)",
            header,
        )
        if duration_match is None:
            raise RuntimeError(f"no retract duration in {source.path}")
        retract_duration_s = float(duration_match.group(1))
        approach_duration_s = float(sum(curve.duration_s for curve in source.curves))
        protocol_duration_s = approach_duration_s + len(source.curves) * retract_duration_s
        midpoint = source.timestamp + timedelta(seconds=0.5 * protocol_duration_s)
        row = rows[source_key(source)]
        row["map_approach_duration_sum_s"] = approach_duration_s
        row["map_retract_duration_per_curve_s"] = retract_duration_s
        row["map_protocol_duration_estimate_s"] = protocol_duration_s
        row["map_protocol_midpoint_estimate"] = midpoint.isoformat()
        row["map_protocol_midpoint_epoch_s"] = midpoint.timestamp()


def build_pair_data(
    sources: list[base.SourceData],
    map_rows: list[dict],
    matrices: dict[str, np.ndarray],
) -> tuple[dict, list[dict], list[dict]]:
    """Build all available same-speed early/late pairs without imputation."""

    map_lookup = {row["source"]: row for row in map_rows}
    pair_data: dict[tuple[int, float, str], dict] = {}
    inventory: list[dict] = []
    curve_rows: list[dict] = []
    for block in range(1, 6):
        for speed in SPEEDS:
            members = sorted(
                [
                    source
                    for source in sources
                    if pilot.map_identity(source.path) == (block, speed)
                ],
                key=lambda source: source.timestamp,
            )
            status = "complete_pair" if len(members) == 2 else "missing_late_map"
            inventory.append(
                {
                    "block": block,
                    "nominal_speed_um_per_s": speed,
                    "map_count": len(members),
                    "pair_status": status,
                    "early_source": source_key(members[0]) if members else "",
                    "late_source": source_key(members[-1]) if len(members) == 2 else "",
                    "used_in_pair_analysis": len(members) == 2,
                }
            )
            if len(members) != 2:
                continue
            early, late = members
            early_key = source_key(early)
            late_key = source_key(late)
            actual_speed = float(
                np.median(
                    [
                        map_lookup[early_key]["gap_speed_20_200nm_median_um_per_s"],
                        map_lookup[late_key]["gap_speed_20_200nm_median_um_per_s"],
                    ]
                )
            )
            midpoint_order = 0.5 * (
                map_lookup[early_key]["acquisition_order"]
                + map_lookup[late_key]["acquisition_order"]
            )
            early_time_s = map_lookup[early_key]["map_protocol_midpoint_epoch_s"]
            late_time_s = map_lookup[late_key]["map_protocol_midpoint_epoch_s"]
            pair_center_time_s = 0.5 * (early_time_s + late_time_s)
            pair_half_span_min = 0.5 * (late_time_s - early_time_s) / 60.0
            for baseline, suffix in (
                ("linear_drift_corrected", "line"),
                ("far_constant_referenced", "constant"),
            ):
                early_matrix = matrices[f"{early_key}|{suffix}"]
                late_matrix = matrices[f"{late_key}|{suffix}"]
                sym = 0.5 * (early_matrix + late_matrix)
                history = 0.5 * (late_matrix - early_matrix)
                item = {
                    "block": block,
                    "speed": speed,
                    "actual_speed": actual_speed,
                    "midpoint_order": midpoint_order,
                    "pair_center_time_epoch_s": pair_center_time_s,
                    "pair_half_span_min": pair_half_span_min,
                    "early_source": early_key,
                    "late_source": late_key,
                    "sym": sym,
                    "history": history,
                }
                pair_data[(block, speed, baseline)] = item
                for index, distance in enumerate(pilot.BIN_CENTERS_NM):
                    for quantity, matrix in (("F_sym", sym), ("F_history", history)):
                        median, q25, q75, count = finite_summary(matrix[:, index])
                        if not count:
                            continue
                        curve_rows.append(
                            {
                                "block": block,
                                "nominal_speed_um_per_s": speed,
                                "actual_gap_speed_um_per_s": actual_speed,
                                "pair_midpoint_acquisition_order": midpoint_order,
                                "pair_center_time_epoch_s": pair_center_time_s,
                                "pair_half_span_min": pair_half_span_min,
                                "baseline_method": baseline,
                                "quantity": quantity,
                                "distance_nm": distance,
                                "median_pN": median,
                                "q25_pN": q25,
                                "q75_pN": q75,
                                "paired_pixel_count": count,
                                "early_source": early_key,
                                "late_source": late_key,
                            }
                        )
    return pair_data, inventory, curve_rows


def build_block_contrasts(pair_data: dict) -> list[dict]:
    """Within-block 0.05-to-0.2 slope, intercept, and midpoint residual."""

    rows: list[dict] = []
    for block in range(1, 6):
        for baseline in ("linear_drift_corrected", "far_constant_referenced"):
            low = pair_data[(block, 0.05, baseline)]
            high = pair_data[(block, 0.2, baseline)]
            delta_u = high["actual_speed"] - low["actual_speed"]
            if not np.isfinite(delta_u) or delta_u <= 0:
                raise RuntimeError(f"invalid speed contrast in block {block}")
            slope_matrix = (high["sym"] - low["sym"]) / delta_u
            intercept_matrix = low["sym"] - slope_matrix * low["actual_speed"]
            midpoint = pair_data.get((block, 0.1, baseline))
            if midpoint is not None:
                predicted_mid = low["sym"] + slope_matrix * (
                    midpoint["actual_speed"] - low["actual_speed"]
                )
                midpoint_residual = midpoint["sym"] - predicted_mid
            else:
                midpoint_residual = np.full_like(slope_matrix, np.nan)
            for index, distance in enumerate(pilot.BIN_CENTERS_NM):
                slope_median, slope_q25, slope_q75, count = finite_summary(
                    slope_matrix[:, index]
                )
                if not count:
                    continue
                f0_median, f0_q25, f0_q75, _ = finite_summary(
                    intercept_matrix[:, index]
                )
                mid_median, mid_q25, mid_q75, mid_count = finite_summary(
                    midpoint_residual[:, index]
                )
                theory = pilot.hydrodynamic_slope_pN_per_um_s(distance)
                rows.append(
                    {
                        "block": block,
                        "baseline_method": baseline,
                        "distance_nm": distance,
                        "low_actual_speed_um_per_s": low["actual_speed"],
                        "high_actual_speed_um_per_s": high["actual_speed"],
                        "speed_slope_median_pN_per_um_s": slope_median,
                        "speed_slope_q25_pN_per_um_s": slope_q25,
                        "speed_slope_q75_pN_per_um_s": slope_q75,
                        "zero_speed_intercept_median_pN": f0_median,
                        "zero_speed_intercept_q25_pN": f0_q25,
                        "zero_speed_intercept_q75_pN": f0_q75,
                        "mid_speed_residual_median_pN": mid_median,
                        "mid_speed_residual_q25_pN": mid_q25,
                        "mid_speed_residual_q75_pN": mid_q75,
                        "mid_speed_pair_present": midpoint is not None,
                        "mid_speed_pixel_count": mid_count,
                        "paired_pixel_count": count,
                        "theory_no_slip_hyd_slope_pN_per_um_s": theory,
                        "observed_to_theory_slope_ratio": slope_median / theory,
                    }
                )
    return rows


def summarize_tests(block_rows: list[dict]) -> list[dict]:
    """Treat blocks, not pixels, as the independent velocity units."""

    rows: list[dict] = []
    for test_name, definition in TESTS.items():
        for baseline in ("linear_drift_corrected", "far_constant_referenced"):
            for distance in pilot.BIN_CENTERS_NM:
                selected = [
                    row
                    for row in block_rows
                    if row["block"] in definition["blocks"]
                    and row["baseline_method"] == baseline
                    and row["distance_nm"] == distance
                ]
                if len(selected) != 3:
                    continue
                slopes = np.asarray(
                    [row["speed_slope_median_pN_per_um_s"] for row in selected]
                )
                intercepts = np.asarray(
                    [row["zero_speed_intercept_median_pN"] for row in selected]
                )
                if not np.all(np.isfinite(slopes)) or not np.all(np.isfinite(intercepts)):
                    continue
                slope_mean = float(np.mean(slopes))
                slope_sd = float(np.std(slopes, ddof=1))
                intercept_mean = float(np.mean(intercepts))
                intercept_sd = float(np.std(intercepts, ddof=1))
                tcrit = float(scipy.stats.t.ppf(0.975, df=2))
                slope_half = tcrit * slope_sd / math.sqrt(3)
                intercept_half = tcrit * intercept_sd / math.sqrt(3)
                theory = float(selected[0]["theory_no_slip_hyd_slope_pN_per_um_s"])
                mid_residuals = np.asarray(
                    [
                        row["mid_speed_residual_median_pN"]
                        for row in selected
                        if row["mid_speed_pair_present"]
                        and np.isfinite(row["mid_speed_residual_median_pN"])
                    ]
                )
                rows.append(
                    {
                        "test": test_name,
                        "primary_test": definition["primary"],
                        "blocks": ",".join(str(value) for value in definition["blocks"]),
                        "baseline_method": baseline,
                        "distance_nm": distance,
                        "block_count": 3,
                        "block_speed_slopes_pN_per_um_s": ";".join(
                            f"{value:.12g}" for value in slopes
                        ),
                        "speed_slope_mean_pN_per_um_s": slope_mean,
                        "speed_slope_sd_pN_per_um_s": slope_sd,
                        "speed_slope_median_pN_per_um_s": float(np.median(slopes)),
                        "speed_slope_95CI_low_pN_per_um_s": slope_mean - slope_half,
                        "speed_slope_95CI_high_pN_per_um_s": slope_mean + slope_half,
                        "positive_slope_block_fraction": float(np.mean(slopes > 0)),
                        "theory_no_slip_hyd_slope_pN_per_um_s": theory,
                        "mean_observed_to_theory_ratio": slope_mean / theory,
                        "zero_speed_force_mean_pN": intercept_mean,
                        "zero_speed_force_sd_pN": intercept_sd,
                        "zero_speed_force_median_pN": float(np.median(intercepts)),
                        "zero_speed_force_95CI_low_pN": intercept_mean - intercept_half,
                        "zero_speed_force_95CI_high_pN": intercept_mean + intercept_half,
                        "mid_speed_complete_blocks": int(mid_residuals.size),
                        "mid_speed_residual_median_across_blocks_pN": float(
                            np.median(mid_residuals)
                        )
                        if mid_residuals.size
                        else float("nan"),
                        "claim_status": "block_level_descriptive_with_n3",
                    }
                )
    return rows


def refresh_discontinuity(pair_data: dict) -> list[dict]:
    """Block3 minus block2 same-speed symmetric force across the refresh boundary."""

    rows: list[dict] = []
    for baseline in ("linear_drift_corrected", "far_constant_referenced"):
        for speed in SPEEDS:
            before = pair_data[(2, speed, baseline)]["sym"]
            after = pair_data[(3, speed, baseline)]["sym"]
            delta = after - before
            for index, distance in enumerate(pilot.BIN_CENTERS_NM):
                median, q25, q75, count = finite_summary(delta[:, index])
                if count:
                    rows.append(
                        {
                            "comparison": "block3_minus_block2_across_liquid_refresh",
                            "baseline_method": baseline,
                            "nominal_speed_um_per_s": speed,
                            "distance_nm": distance,
                            "force_jump_median_pN": median,
                            "force_jump_q25_pN": q25,
                            "force_jump_q75_pN": q75,
                            "paired_pixel_count": count,
                            "interpretation": "refresh_plus_elapsed_time_plus_block_order; not pure liquid causal effect",
                        }
                    )
    return rows


def force_lookup(map_force_rows: list[dict]) -> dict:
    return {
        (row["source"], row["baseline_method"], float(row["distance_nm"])): row[
            "force_median_pN"
        ]
        for row in map_force_rows
    }


def subset_trends(
    name: str,
    blocks: tuple[int, ...],
    map_rows: list[dict],
    map_force_rows: list[dict],
) -> list[dict]:
    selected_maps = [row for row in map_rows if row["block"] in blocks]
    selected_sources = {row["source"] for row in selected_maps}
    selected_force = [row for row in map_force_rows if row["source"] in selected_sources]
    _, trends = pilot.systematic_diagnostics(selected_maps, selected_force)
    for row in trends:
        row["subset"] = name
        row["blocks"] = ",".join(str(value) for value in blocks)
    return trends


PAIR_QC_FIELDS = (
    ("map_contact_invOLS_median_nm_per_V", "InvOLS", "nm/V"),
    ("terminal_load_median_nN", "terminal load", "nN"),
    ("far_slope_median_pN_per_100nm", "far slope", "pN/100 nm"),
    ("contact_height_median_um", "contact height", "um"),
    ("approach_snap_detected_fraction", "approach snap-detected fraction", "fraction"),
    ("approach_snap_distance_median_nm", "approach snap distance", "nm"),
    (
        "retract_detachment_piezo_travel_median_nm",
        "retract detachment travel",
        "nm",
    ),
    ("retract_pull_off_force_median_nN", "retract pull-off", "nN"),
    (
        "branch_contact_height_difference_median_nm",
        "approach-retract contact-height difference",
        "nm",
    ),
)


def build_pair_qc(map_rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """Palindrome-pair contact-state summaries and low-to-high contrasts."""

    pair_rows: list[dict] = []
    for block in range(1, 6):
        for speed in SPEEDS:
            members = sorted(
                [
                    row
                    for row in map_rows
                    if row["block"] == block
                    and row["nominal_speed_um_per_s"] == speed
                ],
                key=lambda row: row["acquisition_order"],
            )
            if len(members) != 2:
                continue
            early, late = members
            actual_speed = float(
                np.median(
                    [
                        early["gap_speed_20_200nm_median_um_per_s"],
                        late["gap_speed_20_200nm_median_um_per_s"],
                    ]
                )
            )
            for field, label, unit in PAIR_QC_FIELDS:
                early_value = float(early[field])
                late_value = float(late[field])
                if not np.isfinite(early_value) or not np.isfinite(late_value):
                    continue
                pair_rows.append(
                    {
                        "block": block,
                        "nominal_speed_um_per_s": speed,
                        "actual_gap_speed_um_per_s": actual_speed,
                        "metric": field,
                        "metric_label": label,
                        "unit": unit,
                        "early_value": early_value,
                        "late_value": late_value,
                        "pair_symmetric_value": 0.5 * (early_value + late_value),
                        "pair_history_half_difference": 0.5
                        * (late_value - early_value),
                        "early_source": early["source"],
                        "late_source": late["source"],
                    }
                )

    lookup = {
        (row["block"], row["nominal_speed_um_per_s"], row["metric"]): row
        for row in pair_rows
    }
    contrast_rows: list[dict] = []
    for block in range(1, 6):
        for field, label, unit in PAIR_QC_FIELDS:
            low = lookup[(block, 0.05, field)]
            high = lookup[(block, 0.2, field)]
            delta_u = (
                high["actual_gap_speed_um_per_s"]
                - low["actual_gap_speed_um_per_s"]
            )
            delta = high["pair_symmetric_value"] - low["pair_symmetric_value"]
            contrast_rows.append(
                {
                    "block": block,
                    "metric": field,
                    "metric_label": label,
                    "unit": unit,
                    "low_actual_speed_um_per_s": low[
                        "actual_gap_speed_um_per_s"
                    ],
                    "high_actual_speed_um_per_s": high[
                        "actual_gap_speed_um_per_s"
                    ],
                    "high_minus_low_pair_symmetric_value": delta,
                    "slope_per_um_per_s": delta / delta_u,
                    "claim_status": "paired_map_QC; diagnostic_not_force_correction",
                }
            )
    return pair_rows, contrast_rows


def _ols_fit(
    design: np.ndarray, response: np.ndarray, speed_column: int
) -> tuple[dict, np.ndarray, np.ndarray]:
    """OLS with classical and HC3 uncertainty plus leverage-aware LOOCV."""

    design = np.asarray(design, dtype=float)
    response = np.asarray(response, dtype=float)
    n, parameters = design.shape
    coefficients, _, rank, singular = np.linalg.lstsq(design, response, rcond=None)
    fitted = design @ coefficients
    residual = response - fitted
    rss = float(residual @ residual)
    residual_df = n - rank
    xtx_inverse = np.linalg.pinv(design.T @ design)
    sigma2 = rss / residual_df if residual_df > 0 else float("nan")
    covariance = sigma2 * xtx_inverse
    classical_se = float(math.sqrt(max(0.0, covariance[speed_column, speed_column])))
    leverage = np.einsum("ij,jk,ik->i", design, xtx_inverse, design)
    hc3_scale = residual**2 / np.maximum(1.0 - leverage, 1e-10) ** 2
    hc3_meat = design.T @ (hc3_scale[:, None] * design)
    hc3_covariance = xtx_inverse @ hc3_meat @ xtx_inverse
    hc3_se = float(math.sqrt(max(0.0, hc3_covariance[speed_column, speed_column])))
    tcrit = float(scipy.stats.t.ppf(0.975, residual_df))
    total = float(np.sum((response - np.mean(response)) ** 2))
    r_squared = float("nan") if total == 0 else 1.0 - rss / total
    aic = n * math.log(max(rss / n, np.finfo(float).tiny)) + 2 * parameters
    aicc = (
        aic + 2 * parameters * (parameters + 1) / (n - parameters - 1)
        if n > parameters + 1
        else float("nan")
    )
    loocv_residual = residual / np.maximum(1.0 - leverage, 1e-10)
    column_norm = np.linalg.norm(design, axis=0)
    scaled = design / np.where(column_norm > 0, column_norm, 1.0)
    result = {
        "map_count": n,
        "parameter_count": parameters,
        "design_rank": int(rank),
        "residual_df": int(residual_df),
        "unit_scaled_design_condition_number": float(np.linalg.cond(scaled)),
        "speed_slope_pN_per_um_s": float(coefficients[speed_column]),
        "speed_slope_classical_SE_pN_per_um_s": classical_se,
        "speed_slope_HC3_SE_pN_per_um_s": hc3_se,
        "speed_slope_HC3_95CI_low_pN_per_um_s": float(
            coefficients[speed_column] - tcrit * hc3_se
        ),
        "speed_slope_HC3_95CI_high_pN_per_um_s": float(
            coefficients[speed_column] + tcrit * hc3_se
        ),
        "residual_RMSE_pN": float(math.sqrt(rss / residual_df)),
        "LOOCV_RMSE_pN": float(math.sqrt(np.mean(loocv_residual**2))),
        "R_squared": r_squared,
        "AICc": aicc,
        "fit_status": "full_rank" if rank == parameters else "rank_deficient",
    }
    return result, fitted, residual


def primary_time_aware_models(
    map_rows: list[dict], map_force_rows: list[dict]
) -> tuple[list[dict], list[dict]]:
    """Separate a smooth sequential time term from the map-level speed term.

    These models are sensitivity analyses, not a randomized causal estimator.
    Their spread is deliberately retained because the strong force decay makes
    the inferred speed term depend on how time is represented.
    """

    maps = sorted(
        [row for row in map_rows if row["block"] in (3, 4, 5)],
        key=lambda row: row["map_protocol_midpoint_epoch_s"],
    )
    time_s = np.asarray([row["map_protocol_midpoint_epoch_s"] for row in maps])
    elapsed_min = (time_s - np.min(time_s)) / 60.0
    time_normalized = elapsed_min - np.mean(elapsed_min)
    time_normalized /= np.max(np.abs(time_normalized))
    speed = np.asarray(
        [row["gap_speed_20_200nm_median_um_per_s"] for row in maps]
    )
    speed_centered = speed - np.mean(speed)
    blocks = np.asarray([row["block"] for row in maps], dtype=int)
    block_matrix = np.column_stack([blocks == block for block in (3, 4, 5)]).astype(
        float
    )
    within_time = np.zeros_like(elapsed_min)
    for block in (3, 4, 5):
        selected = blocks == block
        within_time[selected] = elapsed_min[selected] - np.mean(elapsed_min[selected])
    within_time /= np.max(np.abs(within_time))
    within_time_squared = within_time**2
    for block in (3, 4, 5):
        selected = blocks == block
        within_time_squared[selected] -= np.mean(within_time_squared[selected])

    designs: list[tuple[str, np.ndarray, int]] = []
    for degree in range(1, 6):
        legendre = np.polynomial.legendre.legvander(time_normalized, degree)[:, 1:]
        design = np.column_stack([np.ones(len(maps)), speed_centered, legendre])
        designs.append((f"global_time_poly{degree}", design, 1))
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
                    [block_matrix, speed_centered, within_time, within_time_squared]
                ),
                3,
            ),
        ]
    )

    force = force_lookup(map_force_rows)
    result_rows: list[dict] = []
    prediction_rows: list[dict] = []
    distances = [
        float(value)
        for value in pilot.BIN_CENTERS_NM
        if 20.0 <= value <= 200.0
    ]
    for baseline in ("linear_drift_corrected", "far_constant_referenced"):
        for distance in distances:
            response = np.asarray(
                [force[(row["source"], baseline, distance)] for row in maps],
                dtype=float,
            )
            if not np.all(np.isfinite(response)):
                raise RuntimeError(
                    f"unsupported primary map-level force at {baseline}, {distance} nm"
                )
            theory = pilot.hydrodynamic_slope_pN_per_um_s(distance)
            for model, design, speed_column in designs:
                fit, fitted, residual = _ols_fit(design, response, speed_column)
                fit.update(
                    {
                        "test": "map3_5_no_refresh",
                        "baseline_method": baseline,
                        "distance_nm": distance,
                        "time_model": model,
                        "theory_no_slip_hyd_slope_pN_per_um_s": theory,
                        "observed_to_theory_slope_ratio": fit[
                            "speed_slope_pN_per_um_s"
                        ]
                        / theory,
                        "time_coordinate": "protocol_midpoint_estimate elapsed minutes",
                        "claim_status": "map_level_time_adjusted_sensitivity; sequential_not_randomized",
                    }
                )
                result_rows.append(fit)
                if distance in TARGET_DISTANCES_NM:
                    for row, time_value, speed_value, observed, prediction, error in zip(
                        maps,
                        elapsed_min,
                        speed,
                        response,
                        fitted,
                        residual,
                        strict=True,
                    ):
                        prediction_rows.append(
                            {
                                "source": row["source"],
                                "acquisition_order": row["acquisition_order"],
                                "block": row["block"],
                                "nominal_speed_um_per_s": row[
                                    "nominal_speed_um_per_s"
                                ],
                                "actual_gap_speed_um_per_s": speed_value,
                                "elapsed_protocol_midpoint_min": time_value,
                                "baseline_method": baseline,
                                "distance_nm": distance,
                                "time_model": model,
                                "observed_force_pN": observed,
                                "fitted_force_pN": prediction,
                                "residual_force_pN": error,
                            }
                        )
    return result_rows, prediction_rows


def make_chronology_figure(map_rows: list[dict], map_force_rows: list[dict]) -> None:
    rows = sorted(map_rows, key=lambda row: row["acquisition_order"])
    lookup = force_lookup(map_force_rows)
    x = np.asarray([row["acquisition_order"] for row in rows])
    fields = (
        ("map_contact_invOLS_median_nm_per_V", "Map InvOLS (nm/V)"),
        ("far_slope_median_pN_per_100nm", "Far slope (pN/100 nm)"),
        ("force_50", "F(50 nm) (pN)"),
        ("force_100", "F(100 nm) (pN)"),
        ("retract_pull_off_force_median_nN", "Pull-off (nN)"),
        ("approach_snap_detected_fraction", "Snap-detected fraction"),
    )
    fig, axes = plt.subplots(3, 2, figsize=(14.0, 11.0), sharex=True)
    for ax, (field, ylabel) in zip(axes.flat, fields, strict=True):
        if field == "force_50":
            y = np.asarray(
                [lookup[(row["source"], "linear_drift_corrected", 50.0)] for row in rows]
            )
        elif field == "force_100":
            y = np.asarray(
                [lookup[(row["source"], "linear_drift_corrected", 100.0)] for row in rows]
            )
        else:
            y = np.asarray([row[field] for row in rows])
        ax.plot(x, y, color="0.55", lw=1.0)
        ax.scatter(
            x,
            y,
            c=[COLORS[row["nominal_speed_um_per_s"]] for row in rows],
            s=36,
            zorder=3,
        )
        for boundary in (6.5, 12.5, 18.5, 24.5):
            ax.axvline(boundary, color="0.75", lw=0.7)
        ax.axvline(12.5, color="#d62828", lw=1.2, ls="--")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2)
    axes[0, 0].text(
        12.65,
        axes[0, 0].get_ylim()[1],
        "liquid refresh\nbefore map3",
        color="#d62828",
        va="top",
        fontsize=9,
    )
    labels = [f"M{row['block']} {row['nominal_speed_um_per_s']:g}" for row in rows]
    for ax in axes[-1]:
        ax.set_xticks(x, labels, rotation=70, ha="right", fontsize=7)
        ax.set_xlabel("Chronological map order")
    handles = [
        plt.Line2D([], [], marker="o", ls="", color=color, label=f"{speed:g} µm/s")
        for speed, color in COLORS.items()
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
        ncol=3,
        frameon=False,
    )
    fig.suptitle("All 29 maps: chronological response, force, and contact-state QC", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    fig.savefig(FIGURES / "all_maps_chronological_QC.png", dpi=220)
    plt.close(fig)


def make_pair_qc_figure(pair_qc_rows: list[dict]) -> None:
    selected_fields = (
        "map_contact_invOLS_median_nm_per_V",
        "terminal_load_median_nN",
        "far_slope_median_pN_per_100nm",
        "approach_snap_detected_fraction",
        "retract_pull_off_force_median_nN",
        "retract_detachment_piezo_travel_median_nm",
    )
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.0))
    for ax, field in zip(axes.flat, selected_fields, strict=True):
        label = next(row["metric_label"] for row in pair_qc_rows if row["metric"] == field)
        unit = next(row["unit"] for row in pair_qc_rows if row["metric"] == field)
        for block in (3, 4, 5):
            rows = sorted(
                [
                    row
                    for row in pair_qc_rows
                    if row["block"] == block and row["metric"] == field
                ],
                key=lambda row: row["actual_gap_speed_um_per_s"],
            )
            ax.plot(
                [row["actual_gap_speed_um_per_s"] for row in rows],
                [row["pair_symmetric_value"] for row in rows],
                marker="o",
                color=BLOCK_COLORS[block],
                lw=1.2,
                label=f"map{block}",
            )
        ax.set_title(label)
        ax.set_xlabel("Actual gap speed (µm/s)")
        ax.set_ylabel(unit)
        ax.grid(alpha=0.2)
    axes[0, 0].legend(ncol=3, fontsize=8)
    fig.suptitle(
        "Primary Test B: palindrome-paired contact-state diagnostics\n"
        "(retract protocol is common; map5 0.1 pair remains missing)"
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "primary_pair_contact_state_QC.png", dpi=220)
    plt.close(fig)


def make_time_aware_figures(
    time_rows: list[dict], prediction_rows: list[dict]
) -> None:
    model_order = (
        "global_time_poly1",
        "global_time_poly2",
        "global_time_poly3",
        "global_time_poly4",
        "global_time_poly5",
        "block_fixed_linear_within",
        "block_fixed_quadratic_within",
    )
    labels = ("time-1", "time-2", "time-3", "time-4", "time-5", "block+t", "block+t+t²")
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5))
    for ax, distance in zip(axes.flat, TARGET_DISTANCES_NM, strict=True):
        for baseline, color, marker, label in (
            ("linear_drift_corrected", "#264653", "o", "far-linear"),
            ("far_constant_referenced", "#8d99ae", "s", "far-constant"),
        ):
            lookup = {
                row["time_model"]: row
                for row in time_rows
                if row["baseline_method"] == baseline
                and row["distance_nm"] == distance
            }
            values = [lookup[model]["speed_slope_pN_per_um_s"] for model in model_order]
            ax.plot(
                np.arange(len(model_order)),
                values,
                marker=marker,
                color=color,
                lw=1.2,
                label=label,
            )
        theory = next(
            row["theory_no_slip_hyd_slope_pN_per_um_s"]
            for row in time_rows
            if row["distance_nm"] == distance
        )
        ax.axhline(theory, color="#d62828", ls="--", lw=1.3, label="no-slip hyd")
        ax.axhline(0.0, color="0.35", lw=0.6)
        ax.set_xticks(np.arange(len(model_order)), labels, rotation=35, ha="right")
        ax.set_title(f"D={distance:g} nm")
        ax.set_ylabel("dF/dU (pN per µm/s)")
        ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Primary Test B: inferred speed term versus time model and force baseline")
    fig.tight_layout()
    fig.savefig(FIGURES / "primary_time_model_speed_sensitivity.png", dpi=220)
    plt.close(fig)

    representative = "global_time_poly4"
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), sharex=True)
    for ax, distance in zip(axes.flat, TARGET_DISTANCES_NM, strict=True):
        rows = sorted(
            [
                row
                for row in prediction_rows
                if row["baseline_method"] == "linear_drift_corrected"
                and row["distance_nm"] == distance
                and row["time_model"] == representative
            ],
            key=lambda row: row["elapsed_protocol_midpoint_min"],
        )
        beta = next(
            row["speed_slope_pN_per_um_s"]
            for row in time_rows
            if row["baseline_method"] == "linear_drift_corrected"
            and row["distance_nm"] == distance
            and row["time_model"] == representative
        )
        elapsed = np.asarray([row["elapsed_protocol_midpoint_min"] for row in rows])
        observed = np.asarray([row["observed_force_pN"] for row in rows])
        speed = np.asarray([row["actual_gap_speed_um_per_s"] for row in rows])
        reference_fit = np.asarray([row["fitted_force_pN"] for row in rows]) + beta * (
            0.1 - speed
        )
        ax.plot(elapsed, reference_fit, color="black", lw=1.6, label="time trend at U=0.1")
        for nominal in SPEEDS:
            selected = np.asarray(
                [row["nominal_speed_um_per_s"] == nominal for row in rows]
            )
            ax.scatter(
                elapsed[selected],
                observed[selected],
                color=COLORS[nominal],
                s=38,
                label=f"{nominal:g} µm/s",
                zorder=3,
            )
        ax.set_title(f"D={distance:g} nm")
        ax.set_ylabel("Map-median force (pN)")
        ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=8, ncol=2)
    for ax in axes[-1]:
        ax.set_xlabel("Elapsed estimated map midpoint from map3 start (min)")
    fig.suptitle("Primary Test B: force relaxation dominates the chronological record")
    fig.tight_layout()
    fig.savefig(FIGURES / "primary_force_relaxation_time_model.png", dpi=220)
    plt.close(fig)


def make_pair_figures(pair_data: dict) -> None:
    for quantity, filename, ylabel in (
        ("sym", "symmetric_force_by_block.png", "F_sym (pN)"),
        ("history", "history_residual_by_block.png", "F_history (pN)"),
    ):
        fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.5))
        for block, ax in zip(range(1, 6), axes.flat, strict=False):
            for speed in SPEEDS:
                item = pair_data.get((block, speed, "linear_drift_corrected"))
                if item is None:
                    ax.text(
                        0.5,
                        0.88,
                        f"{speed:g} µm/s pair missing",
                        transform=ax.transAxes,
                        ha="center",
                        color=COLORS[speed],
                    )
                    continue
                medians, q25, q75 = finite_column_quantiles(item[quantity])
                mask = (
                    (pilot.BIN_CENTERS_NM >= 20)
                    & (pilot.BIN_CENTERS_NM <= 200)
                    & np.isfinite(medians)
                )
                ax.fill_between(
                    pilot.BIN_CENTERS_NM[mask],
                    q25[mask],
                    q75[mask],
                    color=COLORS[speed],
                    alpha=0.13,
                )
                ax.plot(
                    pilot.BIN_CENTERS_NM[mask],
                    medians[mask],
                    color=COLORS[speed],
                    lw=1.6,
                    label=f"{speed:g} µm/s",
                )
            ax.axhline(0, color="0.3", lw=0.6)
            ax.set_xlim(20, 200)
            ax.set_title(f"map{block} palindrome block")
            ax.set_xlabel("Separation D (nm)")
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.2)
            ax.legend(fontsize=8)
        axes.flat[-1].axis("off")
        fig.suptitle(
            "Same-pixel early/late palindrome decomposition; shaded band = pixel IQR"
        )
        fig.tight_layout()
        fig.savefig(FIGURES / filename, dpi=220)
        plt.close(fig)


def make_velocity_summary_figures(
    block_rows: list[dict], test_rows: list[dict]
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11.5, 9.0), sharex=True)
    for ax, (test_name, definition) in zip(axes, TESTS.items(), strict=True):
        for block in definition["blocks"]:
            selected = [
                row
                for row in block_rows
                if row["block"] == block
                and row["baseline_method"] == "linear_drift_corrected"
                and 20 <= row["distance_nm"] <= 200
            ]
            ax.plot(
                [row["distance_nm"] for row in selected],
                [row["speed_slope_median_pN_per_um_s"] for row in selected],
                color=BLOCK_COLORS[block],
                lw=1.0,
                alpha=0.75,
                label=f"map{block}",
            )
        summary = [
            row
            for row in test_rows
            if row["test"] == test_name
            and row["baseline_method"] == "linear_drift_corrected"
            and 20 <= row["distance_nm"] <= 200
        ]
        ax.plot(
            [row["distance_nm"] for row in summary],
            [row["speed_slope_mean_pN_per_um_s"] for row in summary],
            color="black",
            lw=2.0,
            label="3-block mean",
        )
        ax.plot(
            [row["distance_nm"] for row in summary],
            [row["theory_no_slip_hyd_slope_pN_per_um_s"] for row in summary],
            color="#d62828",
            lw=1.5,
            ls="--",
            label="no-slip hyd theory",
        )
        ax.axhline(0, color="0.3", lw=0.6)
        ax.set_ylabel("dF/dU (pN per µm/s)")
        ax.set_title(f"{test_name}: {definition['description']}")
        ax.grid(alpha=0.2)
        ax.legend(ncol=3, fontsize=8)
    axes[-1].set_xlabel("Separation D (nm)")
    fig.suptitle("Within-block 0.05→0.2 µm/s velocity contrast")
    fig.tight_layout()
    fig.savefig(FIGURES / "velocity_slope_by_test.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(11.5, 9.0), sharex=True)
    for ax, (test_name, definition) in zip(axes, TESTS.items(), strict=True):
        for block in definition["blocks"]:
            selected = [
                row
                for row in block_rows
                if row["block"] == block
                and row["baseline_method"] == "linear_drift_corrected"
                and 20 <= row["distance_nm"] <= 200
            ]
            ax.plot(
                [row["distance_nm"] for row in selected],
                [row["zero_speed_intercept_median_pN"] for row in selected],
                color=BLOCK_COLORS[block],
                lw=1.1,
                label=f"map{block}",
            )
        summary = [
            row
            for row in test_rows
            if row["test"] == test_name
            and row["baseline_method"] == "linear_drift_corrected"
            and 20 <= row["distance_nm"] <= 200
        ]
        ax.plot(
            [row["distance_nm"] for row in summary],
            [row["zero_speed_force_mean_pN"] for row in summary],
            color="black",
            lw=2.0,
            label="3-block mean",
        )
        ax.axhline(0, color="0.3", lw=0.6)
        ax.set_ylabel("Endpoint-linear U→0 intercept (pN)")
        ax.set_title(test_name)
        ax.grid(alpha=0.2)
        ax.legend()
    axes[-1].set_xlabel("Separation D (nm)")
    fig.suptitle("Provisional blockwise zero-speed intercepts")
    fig.tight_layout()
    fig.savefig(FIGURES / "zero_speed_intercepts_by_test.png", dpi=220)
    plt.close(fig)


def make_midpoint_residual_figure(block_rows: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    for block in range(1, 5):
        rows = [
            row
            for row in block_rows
            if row["block"] == block
            and row["baseline_method"] == "linear_drift_corrected"
            and 20 <= row["distance_nm"] <= 200
        ]
        ax.plot(
            [row["distance_nm"] for row in rows],
            [row["mid_speed_residual_median_pN"] for row in rows],
            color=BLOCK_COLORS[block],
            lw=1.4,
            label=f"map{block}",
        )
    ax.axhline(0, color="0.25", lw=0.7)
    ax.set_xlabel("Separation D (nm)")
    ax.set_ylabel("F(0.1) − endpoint-linear prediction (pN)")
    ax.set_title("Mid-speed linearity diagnostic; map5 missing 0.1 pair")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "mid_speed_linearity_residual.png", dpi=220)
    plt.close(fig)


def make_refresh_figure(refresh_rows: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    for speed in SPEEDS:
        rows = [
            row
            for row in refresh_rows
            if row["baseline_method"] == "linear_drift_corrected"
            and row["nominal_speed_um_per_s"] == speed
            and 20 <= row["distance_nm"] <= 200
        ]
        ax.plot(
            [row["distance_nm"] for row in rows],
            [row["force_jump_median_pN"] for row in rows],
            color=COLORS[speed],
            lw=1.6,
            label=f"{speed:g} µm/s",
        )
    ax.axhline(0, color="0.25", lw=0.7)
    ax.set_xlabel("Separation D (nm)")
    ax.set_ylabel("F_sym(map3) − F_sym(map2) (pN)")
    ax.set_title("Observed discontinuity across same-composition liquid refresh")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "liquid_refresh_discontinuity.png", dpi=220)
    plt.close(fig)


def make_spatial_slope_figures(curve_rows: list[dict], map_rows: list[dict]) -> None:
    all_values = np.asarray([row["far_slope_pN_per_100nm"] for row in curve_rows])
    limit = float(np.nanquantile(np.abs(all_values), 0.98))
    for block in range(1, 6):
        maps = sorted(
            [row for row in map_rows if row["block"] == block],
            key=lambda row: row["timestamp"],
        )
        fig, axes = plt.subplots(2, 3, figsize=(11.3, 7.5), layout="constrained")
        image = None
        for ax in axes.flat:
            ax.axis("off")
        for ax, map_row in zip(axes.flat, maps, strict=False):
            grid = np.full((8, 8), np.nan)
            selected = [row for row in curve_rows if row["source"] == map_row["source"]]
            for row in selected:
                grid[int(row["row"]), int(row["column"])] = row[
                    "far_slope_pN_per_100nm"
                ]
            image = ax.imshow(
                grid, cmap="coolwarm", vmin=-limit, vmax=limit, origin="lower"
            )
            ax.axis("on")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(
                f"#{map_row['acquisition_order']}  {map_row['nominal_speed_um_per_s']:g} µm/s"
            )
        if image is not None:
            fig.colorbar(
                image,
                ax=axes.ravel().tolist(),
                shrink=0.78,
                pad=0.02,
                label="pN per 100 nm",
            )
        fig.suptitle(f"map{block}: 8×8 far-field linear slope maps")
        fig.savefig(FIGURES / f"far_field_slope_maps_map{block}.png", dpi=220)
        plt.close(fig)


def fixed_pixel_outputs(
    sources: list[base.SourceData], matrices: dict[str, np.ndarray]
) -> list[dict]:
    rows: list[dict] = []
    fig, axes = plt.subplots(5, 6, figsize=(17.0, 13.2), sharex=True, sharey=True)
    for ax in axes.flat:
        ax.axis("off")
    visible: list[float] = []
    records = []
    for order, source in enumerate(sources, start=1):
        block, speed = pilot.map_identity(source.path)
        point = next(
            index
            for index in range(64)
            if base.map_pixel_from_index(index, 8, source.map_back_and_forth) == (3, 3)
        )
        force = matrices[f"{source_key(source)}|line"][point]
        records.append((order, block, speed, source, point, force))
        mask = (pilot.BIN_CENTERS_NM >= 20) & (pilot.BIN_CENTERS_NM <= 200)
        visible.extend(force[mask][np.isfinite(force[mask])].tolist())
        for distance in TARGET_DISTANCES_NM:
            index = int(np.flatnonzero(pilot.BIN_CENTERS_NM == distance)[0])
            rows.append(
                {
                    "source": source_key(source),
                    "timestamp": source.timestamp.isoformat(),
                    "acquisition_order": order,
                    "block": block,
                    "nominal_speed_um_per_s": speed,
                    "physical_row_zero_based": 3,
                    "physical_column_zero_based": 3,
                    "point_index_in_archive": point,
                    "distance_nm": distance,
                    "force_linear_drift_corrected_pN": force[index],
                }
            )
    y_min = min(-30.0, float(np.quantile(visible, 0.005)))
    y_max = float(np.quantile(visible, 0.995) * 1.08)
    for ax, (order, block, speed, source, point, force) in zip(
        axes.flat, records, strict=False
    ):
        mask = (
            (pilot.BIN_CENTERS_NM >= 20)
            & (pilot.BIN_CENTERS_NM <= 200)
            & np.isfinite(force)
        )
        ax.axis("on")
        ax.plot(
            pilot.BIN_CENTERS_NM[mask], force[mask], color=COLORS[speed], lw=1.2
        )
        ax.axhline(0, color="0.4", lw=0.5)
        ax.set_xlim(20, 200)
        ax.set_ylim(y_min, y_max)
        ax.set_title(f"#{order} M{block} {speed:g}", fontsize=8)
        ax.grid(alpha=0.15)
        ax.tick_params(labelsize=6)
    fig.suptitle("Fixed pixel row3/column3: all maps in chronological order")
    fig.supxlabel("Separation D (nm)")
    fig.supylabel("Force (pN)")
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.97))
    fig.savefig(FIGURES / "fixed_pixel_all_29_maps_chronological_FD.png", dpi=220)
    plt.close(fig)
    return rows


def report_text(
    calibration: list[dict],
    water_sensitivity: float,
    map_rows: list[dict],
    pair_inventory: list[dict],
    pair_curve_rows: list[dict],
    block_rows: list[dict],
    test_rows: list[dict],
    pair_qc_contrast_rows: list[dict],
    time_rows: list[dict],
    refresh_rows: list[dict],
    trend_rows: list[dict],
    baseline_rows: list[dict],
) -> str:
    d4 = next(row for row in calibration if row["cantilever"] == "D4")
    primary = {
        row["distance_nm"]: row
        for row in test_rows
        if row["test"] == "map3_5_no_refresh"
        and row["baseline_method"] == "linear_drift_corrected"
        and row["distance_nm"] in TARGET_DISTANCES_NM
    }
    primary_constant = {
        row["distance_nm"]: row
        for row in test_rows
        if row["test"] == "map3_5_no_refresh"
        and row["baseline_method"] == "far_constant_referenced"
        and row["distance_nm"] in TARGET_DISTANCES_NM
    }
    block_lookup = {
        (row["block"], row["distance_nm"]): row
        for row in block_rows
        if row["baseline_method"] == "linear_drift_corrected"
        and row["distance_nm"] in TARGET_DISTANCES_NM
    }
    trend_lookup = {
        (row["subset"], row["metric"]): row for row in trend_rows
    }
    time_sensitivity_models = {
        "global_time_poly2",
        "global_time_poly3",
        "global_time_poly4",
        "global_time_poly5",
        "block_fixed_linear_within",
        "block_fixed_quadratic_within",
    }

    def time_slope_range(distance: float, baseline: str) -> tuple[float, float]:
        values = np.asarray(
            [
                row["speed_slope_pN_per_um_s"]
                for row in time_rows
                if row["distance_nm"] == distance
                and row["baseline_method"] == baseline
                and row["time_model"] in time_sensitivity_models
            ]
        )
        return float(np.min(values)), float(np.max(values))

    primary_maps = sorted(
        [row for row in map_rows if row["block"] in (3, 4, 5)],
        key=lambda row: row["acquisition_order"],
    )
    baseline_lookup = {
        (row["source"], row["distance_nm"]): row for row in baseline_rows
    }
    history_50 = np.asarray(
        [
            row["median_pN"]
            for row in pair_curve_rows
            if row["block"] in (3, 4, 5)
            and row["baseline_method"] == "linear_drift_corrected"
            and row["quantity"] == "F_history"
            and row["distance_nm"] == 50.0
        ]
    )
    history_negative = int(np.count_nonzero(history_50 < 0))
    history_sign_p = float(
        scipy.stats.binomtest(history_negative, history_50.size, 0.5).pvalue
    )
    qc_contrast_lookup = {
        (row["block"], row["metric"]): row
        for row in pair_qc_contrast_rows
        if row["block"] in (3, 4, 5)
    }

    def qc_deltas(metric: str, scale: float = 1.0) -> list[float]:
        return [
            qc_contrast_lookup[(block, metric)][
                "high_minus_low_pair_symmetric_value"
            ]
            * scale
            for block in (3, 4, 5)
        ]
    speed_timing = {}
    for speed in SPEEDS:
        members = [row for row in map_rows if row["nominal_speed_um_per_s"] == speed]
        gap_speeds = np.asarray(
            [row["gap_speed_20_200nm_median_um_per_s"] for row in members]
        )
        durations = np.asarray(
            [row["map_protocol_duration_estimate_s"] for row in members]
        )
        speed_timing[speed] = {
            "median": float(np.median(gap_speeds)),
            "minimum": float(np.min(gap_speeds)),
            "maximum": float(np.max(gap_speeds)),
            "duration_min": float(np.median(durations) / 60.0),
        }
    lines = [
        "# 27-08-26 D4 纯水回文测试：完整分析",
        "",
        "## 实验分组与缺失数据",
        "",
        "- **Test A / map1–3**：三个轮换回文block构成完整测试；同浓度液体在map2与map3之间刷新。因此它用于量化换液/时间不连续性，不作为无干预速度因果的primary test。",
        "- **Test B / map3–5**：map3、map4、map5构成无中途换液的primary test；map3同时是两套测试的共享bridge block。",
        "- map5少保存最后一张0.1 µm/s map。分析没有补值：map5的0.1早/晚对标为missing；所有block都完整具有0.05与0.2早/晚对，因此primary速度估计使用block内0.05→0.2差分。map3/map4的0.1对用于检验线性。",
        "",
        "回文顺序：",
        "",
    ]
    for block in range(1, 6):
        members = sorted(
            [row for row in map_rows if row["block"] == block],
            key=lambda row: row["timestamp"],
        )
        order = " → ".join(f"{row['nominal_speed_um_per_s']:g}" for row in members)
        lines.append(f"- map{block}: `{order} µm/s`。")
    lines.extend(["", "实际20–200 nm gap speed和单张map protocol时长为：", ""])
    for speed in SPEEDS:
        item = speed_timing[speed]
        lines.append(
            f"- nominal {speed:g} µm/s：actual median {item['median']:.4f} µm/s，"
            f"range {item['minimum']:.4f}–{item['maximum']:.4f} µm/s；protocol约 {item['duration_min']:.2f} min/map。"
        )
    lines.extend(
        [
            "",
            "## 力标定",
            "",
            f"温度采用用户给定 **{TEMPERATURE_C:.1f} °C**。D4空气热标定给出 `k={d4['spring_constant_N_per_m']:.6f} N/m`；全部29张水中map的有效硬接触共同给出 global `InvOLS={water_sensitivity*1e9:.3f} nm/V`。文件内写入的sensitivity和force conversion未用于最终力。",
            "",
            "| Cantilever | Air InvOLS (nm/V) | k (N/m) | f0 (kHz) | Q | unique TND / files |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in calibration:
        lines.append(
            f"| {row['cantilever']} | {row['sensitivity_nm_per_V']:.2f} ± {row['sensitivity_repeatability_sd_nm_per_V']:.2f} | "
            f"{row['spring_constant_N_per_m']:.4f} ± {row['spring_constant_repeatability_sd_N_per_m']:.4f} | "
            f"{row['resonance_frequency_kHz']:.4f} | {row['quality_factor']:.2f} | "
            f"{row['thermal_spectra']} / {row['thermal_files']} |"
        )
    map_inv_values = np.asarray(
        [row["map_contact_invOLS_median_nm_per_V"] for row in map_rows]
    )
    lines.extend(
        [
            "",
            "D4有两个byte-identical TND文件；inventory保留两者，但SHA-256 dedup后只用4个unique spectra计算thermal mean和repeatability。该修正使D4 k相对重复计权结果改变约0.12%。",
            f"表中k的repeatability SD对D4为 {100*d4['spring_constant_repeatability_sd_N_per_m']/d4['spring_constant_N_per_m']:.1f}%；它是所有force共享的乘法scale uncertainty，未并入每个map/pixel的IQR。水中local map InvOLS范围为 {np.min(map_inv_values):.3f}–{np.max(map_inv_values):.3f} nm/V；按既定规则仍对全部map使用同一个global值，local值只用于QC。",
        ]
    )
    lines.extend(
        [
            "",
            "## Primary Test B：无换液block内速度结果",
            "",
            "每个map block先对相同速度的早/晚同pixel曲线计算 `F_sym=(F_early+F_late)/2` 和 `F_history=(F_late−F_early)/2`。随后在每个block内用0.05与0.2 µm/s的实际gap speed计算slope和endpoint-linear U→0 intercept。下表的±SD来自map3/map4/map5三个block，不来自64个pixels。",
            "",
            "| D (nm) | slopes map3 / map4 / map5 (pN/(µm/s)) | 3-block mean ± SD | hyd theory | ratio | U→0 force mean ± SD (pN) | mid-speed residual* (pN) |",
            "|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for distance in TARGET_DISTANCES_NM:
        row = primary[distance]
        slopes = [
            block_lookup[(block, distance)]["speed_slope_median_pN_per_um_s"]
            for block in (3, 4, 5)
        ]
        lines.append(
            f"| {distance:.0f} | {slopes[0]:.1f} / {slopes[1]:.1f} / {slopes[2]:.1f} | "
            f"{row['speed_slope_mean_pN_per_um_s']:.1f} ± {row['speed_slope_sd_pN_per_um_s']:.1f} | "
            f"{row['theory_no_slip_hyd_slope_pN_per_um_s']:.2f} | {row['mean_observed_to_theory_ratio']:.1f} | "
            f"{row['zero_speed_force_mean_pN']:.1f} ± {row['zero_speed_force_sd_pN']:.1f} | "
            f"{row['mid_speed_residual_median_across_blocks_pN']:.1f} |"
        )
    lines.extend(
        [
            "",
            "注：mid-speed residual仅来自具有完整0.1 pair的map3/map4，定义为实际F_sym(0.1)减去0.05–0.2 endpoint-linear prediction。",
            "",
            "如果是真正no-slip drainage，dF/dU应为正、跨block接近一致，并按1/D衰减。三个block的斜率在20 nm出现两种符号，在50–100 nm也有map5负斜率；因此上面的U→0 intercept是provisional calculation，不是validated equilibrium force。",
            "",
            "### 时间模型和baseline敏感性",
            "",
            "回文平均只严格抵消关于pair中心近似线性的时间项。本数据的force relaxation明显弯曲，因此又在17张原始map层级拟合 `F=平滑时间项+βU`；下表给出二至五次时间多项式以及block-fixed线性/二次时间模型得到的β范围。范围不是置信区间，而是model sensitivity。",
            "",
            "| D (nm) | block-pair mean β: far-linear | block-pair mean β: far-constant | time-aware β range: far-linear | time-aware β range: far-constant | no-slip hyd |",
            "|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for distance in TARGET_DISTANCES_NM:
        linear_range = time_slope_range(distance, "linear_drift_corrected")
        constant_range = time_slope_range(distance, "far_constant_referenced")
        lines.append(
            f"| {distance:.0f} | {primary[distance]['speed_slope_mean_pN_per_um_s']:.1f} | "
            f"{primary_constant[distance]['speed_slope_mean_pN_per_um_s']:.1f} | "
            f"{linear_range[0]:.1f}…{linear_range[1]:.1f} | "
            f"{constant_range[0]:.1f}…{constant_range[1]:.1f} | "
            f"{primary[distance]['theory_no_slip_hyd_slope_pN_per_um_s']:.2f} |"
        )
    lines.extend(
        [
            "",
            "20 nm的time-aware速度项在两种baseline和全部时间模型下均为负，与drainage方向相反；50 nm随时间模型从略负到远高于hyd theory，跨过零。100–200 nm又对far-field零点定义发生量级甚至符号变化。所有target distance、两种baseline和全部time model的HC3 95% CI都包含零。因此当前数据没有任何距离区间满足可靠hyd subtraction和U→0 extrapolation所需的统计、model和baseline稳定性。",
            "",
            "## 换液边界",
            "",
            "以下为相同速度的F_sym(map3)−F_sym(map2)，包含换液、elapsed time和不同回文block order，不能称为纯换液因果效应：",
            "",
            "| D (nm) | 0.05 µm/s | 0.1 µm/s | 0.2 µm/s |",
            "|---:|---:|---:|---:|",
        ]
    )
    for distance in TARGET_DISTANCES_NM:
        values = []
        for speed in SPEEDS:
            row = next(
                item
                for item in refresh_rows
                if item["baseline_method"] == "linear_drift_corrected"
                and item["nominal_speed_um_per_s"] == speed
                and item["distance_nm"] == distance
            )
            values.append(row["force_jump_median_pN"])
        lines.append(
            f"| {distance:.0f} | {values[0]:+.1f} pN | {values[1]:+.1f} pN | {values[2]:+.1f} pN |"
        )
    no_refresh_pull = trend_lookup[("map3_5_no_refresh", "retract_pull_off_force_nN")]
    no_refresh_inv = trend_lookup[("map3_5_no_refresh", "water_InvOLS_nm_per_V")]
    no_refresh_far = trend_lookup[("map3_5_no_refresh", "far_slope_pN_per_100nm")]
    lines.extend(
        [
            "",
            "## 无换液阶段仍存在的时间漂移",
            "",
            f"在map3–5的17张连续map中，map-median InvOLS与order的Spearman `ρ={no_refresh_inv['spearman_rho_vs_order']:.3f}`；far slope `ρ={no_refresh_far['spearman_rho_vs_order']:.3f}`；pull-off `ρ={no_refresh_pull['spearman_rho_vs_order']:.3f}`。这些是顺序采集的描述性趋势，p值不等于随机化速度因果检验。",
            "",
            "从Test B第一张map3到最后一张map5，far-linear map median变化为：",
            "",
            "| D (nm) | first (pN) | last (pN) | relative change |",
            "|---:|---:|---:|---:|",
        ]
    )
    for distance in TARGET_DISTANCES_NM:
        first = baseline_lookup[(primary_maps[0]["source"], distance)][
            "linear_drift_corrected_force_pN"
        ]
        last = baseline_lookup[(primary_maps[-1]["source"], distance)][
            "linear_drift_corrected_force_pN"
        ]
        relative = 100.0 * (last - first) / abs(first)
        lines.append(
            f"| {distance:.0f} | {first:.1f} | {last:.1f} | {relative:+.1f}% |"
        )
    inv_first = primary_maps[0]["map_contact_invOLS_median_nm_per_V"]
    inv_last = primary_maps[-1]["map_contact_invOLS_median_nm_per_V"]
    lines.extend(
        [
            "",
            f"同期local contact InvOLS只从 {inv_first:.3f} 变到 {inv_last:.3f} nm/V（{100*(inv_last-inv_first)/inv_first:+.2f}%），terminal load也维持约21.8 nN；因此20–100 nm的巨大force relaxation不能由sensitivity gain或加载力变化解释。50 nm处完整的8个early/late pair全部为 `F_history<0`，即later map的force更小；在假定pair-sign独立时exact two-sided sign test `p={history_sign_p:.4f}`。这些pair共享同一条顺序relaxation，因此p值只作方向性描述；8/8同号本身已经说明无换液并不等于stationary。",
            "",
            "## Approach速度还改变了contact/retract状态",
            "",
        ]
    )
    inv_deltas = qc_deltas("map_contact_invOLS_median_nm_per_V")
    load_deltas = qc_deltas("terminal_load_median_nN")
    contact_deltas_nm = qc_deltas("contact_height_median_um", 1000.0)
    pull_deltas = qc_deltas("retract_pull_off_force_median_nN")
    snap_fraction_deltas = qc_deltas("approach_snap_detected_fraction")
    lines.extend(
        [
            "这里比较的是每个block内0.2 pair-symmetric值减去0.05 pair-symmetric值；map3/map4/map5分别为：",
            "",
            f"- InvOLS：{inv_deltas[0]:+.3f} / {inv_deltas[1]:+.3f} / {inv_deltas[2]:+.3f} nm/V；不足以解释force差。",
            f"- terminal load：{load_deltas[0]:+.3f} / {load_deltas[1]:+.3f} / {load_deltas[2]:+.3f} nN；相对约21.8 nN很小。",
            f"- apparent contact height：{contact_deltas_nm[0]:+.1f} / {contact_deltas_nm[1]:+.1f} / {contact_deltas_nm[2]:+.1f} nm；符号不一致，不支持一个固定的speed-dependent contact-zero shift。",
            f"- approach snap-detected fraction：{snap_fraction_deltas[0]:+.3f} / {snap_fraction_deltas[1]:+.3f} / {snap_fraction_deltas[2]:+.3f}。高速度在map3/map4显著降低检测率，而且被检测事件的apparent distance可到约90–110 nm；这更像branch-shape/threshold QC异常，不能直接当真实snap-in位置。",
            f"- retract pull-off：{pull_deltas[0]:+.3f} / {pull_deltas[1]:+.3f} / {pull_deltas[2]:+.3f} nN，三个block都变得更不负。retract protocol在所有map中相同，所以approach condition与后续接触/脱离状态明确相关；因速度顺序不是随机化的，尚不能区分真实approach-history因果与residual time confounding，但两者都不是可直接相减的approach hydrodynamic force。",
            "",
            "## Far-field零点敏感性",
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
            f"- D={distance:.0f} nm：linear-corrected − constant-referenced的29-map median为 {np.median(differences):+.2f} pN，范围 {np.min(differences):+.2f}到{np.max(differences):+.2f} pN。"
        )
    lines.extend(
        [
            "",
            "## 结论边界",
            "",
            "- map是速度/时间实验单位；8×8 pixels是同位置paired observations，未被当作独立速度重复。",
            "- map5缺失0.1 late map不会影响三block的0.05→0.2 endpoint contrast；但它破坏了完整三速度回文的linearity/curvature诊断，且没有补值。lone early 0.1 map只进入chronology/time-aware model。",
            "- Test A在同浓度refresh边界观察到force/contact-state不连续，但它与elapsed time和block order共变，不能量化纯refresh因果；map1–3因此不用于无干预速度因果的primary结论。",
            "- Test B虽是目前最强的速度证据，但force随时间强烈、非线性衰减；block斜率不一致，100–200 nm又对baseline发生符号变化。当前结论是 **hyd slope不可唯一识别，不能据此做hyd subtraction或U→0 equilibrium recovery**。",
            "- 20 nm的time-aware斜率方向与hyd相反，50 nm则强烈依赖时间模型；同时fixed-retract pull-off随approach speed改变。这支持额外的approach-history/contact-state systematic，而不是把速度差解释成单一drainage项。",
            "- 当前输出不做PB/zeta/Debye拟合；先判定哪些距离的U→0 force可辨识，之后才应进入表面电势拟合。",
            "",
            "## 输出",
            "",
            "- `map_inventory_QC.csv`, `pixel_QC.csv`, `map_force_curves.csv`: 29-map raw reconstruction and QC.",
            "- `pair_inventory.csv`, `palindrome_pair_curves.csv`: explicit missingness and early/late decomposition.",
            "- `block_velocity_contrasts.csv`, `test_velocity_summary.csv`: block-level slopes and provisional U→0 intercepts.",
            "- `primary_time_aware_speed_models.csv`, `primary_time_aware_model_predictions.csv`: map-level time/speed sensitivity models and fitted chronology.",
            "- `pair_contact_state_QC.csv`, `pair_contact_state_velocity_contrasts.csv`: paired InvOLS/load/contact/snap/retract diagnostics.",
            "- `liquid_refresh_discontinuity.csv`, `chronological_trends.csv`, `baseline_sensitivity_slices.csv`: intervention/time/baseline systematics.",
            "- `fixed_pixel_row3_col3_slices.csv`, `pixel_force_slices_20_50_100_200nm.csv`: same-position examples and spatial target slices.",
            "- `figures/`: chronology, block pair curves, velocity slopes, zero-speed intercepts, refresh jump, midpoint residual, fixed pixel, and 8×8 far-slope maps.",
            "- `provenance.json`, `artifact_manifest.sha256`: complete raw hashes, parameters, software, and artifact identities.",
            "",
        ]
    )
    return "\n".join(lines)


def create_manifest(paths: list[Path], destination: Path) -> None:
    lines = [
        f"{sha256_file(path)}  {path.relative_to(ROOT)}"
        for path in sorted(paths, key=lambda path: str(path.relative_to(ROOT)))
    ]
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    # Re-route reusable plot functions to the full-analysis bundle.
    pilot.RESULTS = RESULTS
    pilot.FIGURES = FIGURES

    calibration, thermal_rows, force_cal_rows, calibration_plots = (
        pilot.calibration_analysis()
    )
    d4 = next(row for row in calibration if row["cantilever"] == "D4")
    paths = sorted((DATA_ROOT / "0").glob("map[1-5]-*.jpk-force-map"))
    sources = [base.load_source(path.resolve(), 0) for path in paths]
    sources.sort(key=lambda source: source.timestamp)
    if len(sources) != 29:
        raise RuntimeError(f"expected 29 maps including incomplete map5, got {len(sources)}")
    sensitivity_rows, water_sensitivity = pilot.set_map_contact_fits(sources)
    curve_rows, map_rows, map_force_rows, matrices = pilot.map_analysis(
        d4["spring_constant_N_per_m"], sources
    )
    add_map_timing_metadata(sources, map_rows)

    pair_data, pair_inventory, pair_curve_rows = build_pair_data(
        sources, map_rows, matrices
    )
    block_rows = build_block_contrasts(pair_data)
    test_rows = summarize_tests(block_rows)
    pair_qc_rows, pair_qc_contrast_rows = build_pair_qc(map_rows)
    time_rows, time_prediction_rows = primary_time_aware_models(
        map_rows, map_force_rows
    )
    refresh_rows = refresh_discontinuity(pair_data)
    baseline_rows, full_trends = pilot.systematic_diagnostics(map_rows, map_force_rows)
    for row in full_trends:
        row["subset"] = "all_maps"
        row["blocks"] = "1,2,3,4,5"
    trend_rows = full_trends
    trend_rows += subset_trends(
        "map1_3_refresh_affected", (1, 2, 3), map_rows, map_force_rows
    )
    trend_rows += subset_trends(
        "map3_5_no_refresh", (3, 4, 5), map_rows, map_force_rows
    )
    fixed_pixel_rows = fixed_pixel_outputs(sources, matrices)

    pixel_slice_rows: list[dict] = []
    for source in sources:
        block, speed = pilot.map_identity(source.path)
        matrix = matrices[f"{source_key(source)}|line"]
        for distance in TARGET_DISTANCES_NM:
            index = int(np.flatnonzero(pilot.BIN_CENTERS_NM == distance)[0])
            for point in range(64):
                row, column = base.map_pixel_from_index(
                    point, 8, source.map_back_and_forth
                )
                pixel_slice_rows.append(
                    {
                        "source": source_key(source),
                        "timestamp": source.timestamp.isoformat(),
                        "block": block,
                        "nominal_speed_um_per_s": speed,
                        "point_index": point,
                        "row": row,
                        "column": column,
                        "distance_nm": distance,
                        "force_linear_drift_corrected_pN": matrix[point, index],
                    }
                )

    csv_outputs = {
        "calibration_summary.csv": calibration,
        "thermal_refits.csv": thermal_rows,
        "calibration_force_contact_fits.csv": force_cal_rows,
        "water_contact_sensitivity_curves.csv": sensitivity_rows,
        "map_inventory_QC.csv": map_rows,
        "pixel_QC.csv": curve_rows,
        "map_force_curves.csv": map_force_rows,
        "pair_inventory.csv": pair_inventory,
        "palindrome_pair_curves.csv": pair_curve_rows,
        "block_velocity_contrasts.csv": block_rows,
        "test_velocity_summary.csv": test_rows,
        "pair_contact_state_QC.csv": pair_qc_rows,
        "pair_contact_state_velocity_contrasts.csv": pair_qc_contrast_rows,
        "primary_time_aware_speed_models.csv": time_rows,
        "primary_time_aware_model_predictions.csv": time_prediction_rows,
        "liquid_refresh_discontinuity.csv": refresh_rows,
        "chronological_trends.csv": trend_rows,
        "baseline_sensitivity_slices.csv": baseline_rows,
        "fixed_pixel_row3_col3_slices.csv": fixed_pixel_rows,
        "pixel_force_slices_20_50_100_200nm.csv": pixel_slice_rows,
    }
    for filename, rows in csv_outputs.items():
        write_csv(RESULTS / filename, rows)

    np.savez_compressed(
        RESULTS / "pixel_force_curves.npz",
        distance_nm=pilot.BIN_CENTERS_NM,
        source=np.asarray([source_key(source) for source in sources]),
        force_linear_drift_corrected_pN=np.asarray(
            [matrices[f"{source_key(source)}|line"] for source in sources]
        ),
        force_far_constant_referenced_pN=np.asarray(
            [matrices[f"{source_key(source)}|constant"] for source in sources]
        ),
    )

    pilot.make_calibration_figure(calibration, calibration_plots)
    make_chronology_figure(map_rows, map_force_rows)
    make_pair_qc_figure(pair_qc_rows)
    make_time_aware_figures(time_rows, time_prediction_rows)
    make_pair_figures(pair_data)
    make_velocity_summary_figures(block_rows, test_rows)
    make_midpoint_residual_figure(block_rows)
    make_refresh_figure(refresh_rows)
    make_spatial_slope_figures(curve_rows, map_rows)

    report = report_text(
        calibration,
        water_sensitivity,
        map_rows,
        pair_inventory,
        pair_curve_rows,
        block_rows,
        test_rows,
        pair_qc_contrast_rows,
        time_rows,
        refresh_rows,
        trend_rows,
        baseline_rows,
    )
    (RESULTS / "REPORT.md").write_text(report, encoding="utf-8")

    raw_paths = sorted(DATA_ROOT.rglob("*.jpk-force"))
    raw_paths += sorted(DATA_ROOT.rglob("*.jpk-force-map"))
    raw_paths += sorted(DATA_ROOT.rglob("*.tnd"))
    provenance = {
        "analysis": "full overlapping map1-3 refresh-affected and map3-5 no-refresh palindrome analysis",
        "temperature_C": TEMPERATURE_C,
        "calibration_environment": "air",
        "calibration_exact_duplicate_handling": "SHA-256-identical TND files retained in inventory but counted once in thermal calibration",
        "map_environment": "nominally pure water, no added salt",
        "geometry": "silica sphere against silica plane",
        "same_region": True,
        "test_definitions": TESTS,
        "liquid_refresh_boundary": "between map2 and map3 per experiment interpretation",
        "missing_data": "map5 final 0.1 um/s map absent; no imputation",
        "primary_velocity_contrast": "within-block paired Fsym 0.05-to-0.2 um/s",
        "time_aware_sensitivity_models": [
            "global Legendre-polynomial time degree 1 through 5 plus speed",
            "block-fixed within-block linear time plus speed",
            "block-fixed within-block quadratic time plus speed",
        ],
        "time_coordinate": "JPK map start plus half stored approach+retract protocol duration; excludes XY/baseline overhead",
        "statistical_unit": "map/block; 8x8 pixels retained as paired spatial observations, not independent speed replicates",
        "probe_radius_m": pilot.PROBE_RADIUS_M,
        "water_viscosity_mPa_s": pilot.WATER_VISCOSITY_MPA_S,
        "water_InvOLS_nm_per_V": water_sensitivity * 1e9,
        "water_InvOLS_method": "global median of retained hard-contact fits from all 29 maps",
        "D4_spring_constant_N_per_m": d4["spring_constant_N_per_m"],
        "force_baselines_retained": [
            "linear_drift_corrected",
            "far_constant_referenced",
        ],
        "bin_centers_nm": pilot.BIN_CENTERS_NM.tolist(),
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256_file(path) for path in raw_paths
        },
    }
    code_paths = sorted(
        {
            Path(__file__).resolve(),
            Path(pilot.__file__).resolve(),
            Path(pilot.cal.__file__).resolve(),
            Path(pilot.events.__file__).resolve(),
            Path(base.__file__).resolve(),
        }
    )
    provenance["code_hashes"] = {
        str(path.relative_to(ROOT)): sha256_file(path) for path in code_paths
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
    artifacts.extend(code_paths)
    create_manifest(artifacts, RESULTS / "artifact_manifest.sha256")
    print(f"Wrote {RESULTS}")
    print(f"maps={len(sources)}, water InvOLS={water_sensitivity*1e9:.6f} nm/V")


if __name__ == "__main__":
    main()

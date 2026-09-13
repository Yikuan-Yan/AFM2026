#!/usr/bin/env python3
"""Speed-conditioned dynamic-baseline analysis for 99.7 wt% glycerol/D3.

The approach motion has a speed-dependent start-up transient.  This script
therefore detects a settled, constant-velocity baseline window in physical
scanner travel separately for each speed.  Candidate windows are selected on
map blocks 1 and 3 and evaluated on held-out blocks 2 and 4.  The baseline of
each individual curve is the median voltage in the selected window; no slope,
contact, hydrodynamic, PB, or surface-force model is fitted.

The result is an incremental signal relative to the in-motion plateau, not an
absolute zero-force signal.  Retract curves are intentionally not corrected.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import platform
import sys
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import analyze_12_09_26_997glycerol_observational as observational  # noqa: E402
import fit_glycerol_surface_forces as base  # noqa: E402


RAW_ROOT = observational.RAW_ROOT
DATA = observational.DATA
DOWNLOAD_MANIFEST = observational.DOWNLOAD_MANIFEST
OUT = ROOT / "analysis" / "glycerol_99p7_D3_adaptive_baseline_results"
FIG = OUT / "figures"

SPEEDS = observational.SPEEDS_UM_PER_S
COLORS = observational.COLORS
GRID_SIZE = observational.GRID_SIZE
MAP_COUNT = observational.MAP_COUNT
MAPS_PER_BLOCK = observational.MAPS_PER_BLOCK
D3_SCALE_NN_PER_V = observational.D3_PRIOR_SCALE_NN_PER_V

BIN_WIDTH_NM = 5.0
BIN_EDGES_NM = np.arange(0.0, 495.0, BIN_WIDTH_NM)
BIN_CENTERS_NM = BIN_EDGES_NM[:-1] + BIN_WIDTH_NM / 2.0

# Candidate geometry is fixed before seeing the held-out blocks.  A common
# 120 nm width prevents the optimizer from winning by choosing a short window.
CANDIDATE_STARTS_NM = np.arange(20.0, 180.0 + 0.1, 5.0)
WINDOW_WIDTH_NM = 120.0
GUARD_WIDTH_NM = 20.0
TRAIN_BLOCKS = (1, 3)
VALIDATION_BLOCKS = (2, 4)
FIXED_REFERENCE_START_NM = 100.0

# Score terms are all in raw volts.  The weights express the operational goal:
# prioritize a flat window, then continuity into the adjacent guards, with a
# small penalty for bin-to-bin roughness.  This is a detector objective, not a
# statistical uncertainty or a physical model likelihood.
DRIFT_WEIGHT = 1.0
GUARD_WEIGHT = 0.5
ROUGHNESS_WEIGHT = 0.25
NEAR_OPTIMAL_RELATIVE_TOLERANCE = 0.01


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
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def finite(values: Iterable[float]) -> np.ndarray:
    array = np.asarray(list(values), dtype=np.float64)
    return array[np.isfinite(array)]


def median(values: Iterable[float]) -> float:
    array = finite(values)
    return float(np.median(array)) if array.size else float("nan")


def quantiles(values: Iterable[float]) -> tuple[float, float, float]:
    array = finite(values)
    if not array.size:
        return float("nan"), float("nan"), float("nan")
    return tuple(float(value) for value in np.percentile(array, [25, 50, 75]))


def serpentine_coordinates(point_index: int, back_and_forth: bool) -> tuple[int, int]:
    return base.map_pixel_from_index(point_index, GRID_SIZE, back_and_forth)


def binned_raw(travel_nm: np.ndarray, deflection_v: np.ndarray) -> np.ndarray:
    output = np.full(BIN_CENTERS_NM.shape, np.nan, dtype=np.float64)
    indices = np.floor(travel_nm / BIN_WIDTH_NM).astype(int)
    for index in np.unique(indices):
        if 0 <= index < output.size:
            output[index] = float(np.median(deflection_v[indices == index]))
    return output


def window_values(x: np.ndarray, y: np.ndarray, start_nm: float, stop_nm: float) -> np.ndarray:
    values = y[(x >= start_nm) & (x < stop_nm)]
    values = values[np.isfinite(values)]
    if values.size < 3:
        raise RuntimeError(f"Window {start_nm:g}-{stop_nm:g} nm has too few samples")
    return values


def candidate_components(binned_v: np.ndarray, start_nm: float) -> dict[str, float]:
    """Return physically interpretable stability terms for one candidate."""

    stop_nm = start_nm + WINDOW_WIDTH_NM
    window = (BIN_CENTERS_NM >= start_nm) & (BIN_CENTERS_NM < stop_nm)
    pre = (BIN_CENTERS_NM >= start_nm - GUARD_WIDTH_NM) & (BIN_CENTERS_NM < start_nm)
    post = (BIN_CENTERS_NM >= stop_nm) & (BIN_CENTERS_NM < stop_nm + GUARD_WIDTH_NM)
    y = binned_v[window]
    pre_y = binned_v[pre]
    post_y = binned_v[post]
    if np.count_nonzero(np.isfinite(y)) < 20 or min(
        np.count_nonzero(np.isfinite(pre_y)), np.count_nonzero(np.isfinite(post_y))
    ) < 3:
        raise RuntimeError("Candidate or guard window has insufficient finite bins")
    y = y[np.isfinite(y)]
    pre_y = pre_y[np.isfinite(pre_y)]
    post_y = post_y[np.isfinite(post_y)]
    centre = float(np.median(y))
    quarter = max(3, y.size // 4)
    drift = abs(float(np.median(y[-quarter:]) - np.median(y[:quarter])))
    pre_mismatch = abs(float(np.median(pre_y) - centre))
    post_mismatch = abs(float(np.median(post_y) - centre))
    roughness = float(np.median(np.abs(np.diff(y))))
    score = (
        DRIFT_WEIGHT * drift
        + GUARD_WEIGHT * (pre_mismatch + post_mismatch)
        + ROUGHNESS_WEIGHT * roughness
    )
    return {
        "drift_V": drift,
        "pre_guard_mismatch_V": pre_mismatch,
        "post_guard_mismatch_V": post_mismatch,
        "roughness_V": roughness,
        "score_V": score,
    }


def aggregate_map_medians(rows: list[dict], field: str, blocks: tuple[int, ...]) -> float:
    map_values: list[float] = []
    for order in sorted({int(row["acquisition_order"]) for row in rows if int(row["block"]) in blocks}):
        selected = [float(row[field]) for row in rows if int(row["acquisition_order"]) == order]
        map_values.append(median(selected))
    return median(map_values)


def select_stable_candidate(candidates: list[dict], score_field: str) -> tuple[dict, set[float]]:
    """Choose the centre of the contiguous flat minimum, not a noisy single bin."""

    ordered = sorted(candidates, key=lambda row: float(row["candidate_start_nm"]))
    scores = np.asarray([float(row[score_field]) for row in ordered])
    minimum_index = int(np.argmin(scores))
    threshold = float(scores[minimum_index] * (1.0 + NEAR_OPTIMAL_RELATIVE_TOLERANCE))
    left = minimum_index
    right = minimum_index
    while left > 0 and scores[left - 1] <= threshold:
        left -= 1
    while right + 1 < len(ordered) and scores[right + 1] <= threshold:
        right += 1
    flat_bottom = ordered[left : right + 1]
    target = 0.5 * (
        float(flat_bottom[0]["candidate_start_nm"])
        + float(flat_bottom[-1]["candidate_start_nm"])
    )
    winner = min(
        flat_bottom,
        key=lambda row: (
            abs(float(row["candidate_start_nm"]) - target),
            float(row["candidate_start_nm"]),
        ),
    )
    return winner, {float(row["candidate_start_nm"]) for row in flat_bottom}


def load_curves() -> tuple[list[dict], list[dict]]:
    paths = sorted(DATA.glob("*.jpk-force-map"))
    if len(paths) != MAP_COUNT:
        raise RuntimeError(f"Expected {MAP_COUNT} maps, found {len(paths)}")
    download = json.loads(DOWNLOAD_MANIFEST.read_text(encoding="utf-8"))
    expected_hashes = {Path(row["local_path"]).name: row["sha256"] for row in download["files"]}
    records: list[dict] = []
    inventory: list[dict] = []
    for acquisition_order, path in enumerate(paths, start=1):
        digest = sha256_file(path)
        if digest != expected_hashes.get(path.name):
            raise RuntimeError(f"Download hash mismatch for {path.name}")
        metadata = observational.map_metadata(path)
        source = base.load_source(path, 100)
        if source.skipped_curves or len(source.curves) != GRID_SIZE**2:
            raise RuntimeError(f"Incomplete approach inventory in {path.name}")
        speed = observational.speed_from_source(source)
        if not np.isclose(speed, metadata["approach_nominal_speed_um_per_s"], atol=0.002):
            raise RuntimeError(f"Speed mismatch in {path.name}")
        by_index = {int(curve.point_index): curve for curve in source.curves}
        if set(by_index) != set(range(GRID_SIZE**2)):
            raise RuntimeError(f"Unexpected point indices in {path.name}")
        block = (acquisition_order - 1) // MAPS_PER_BLOCK + 1
        position_in_block = (acquisition_order - 1) % MAPS_PER_BLOCK + 1
        for point_index in range(GRID_SIZE**2):
            curve = by_index[point_index]
            travel_nm = (curve.measured_height_m[0] - curve.measured_height_m) * 1e9
            if travel_nm[-1] <= travel_nm[0] or float(np.min(np.diff(travel_nm))) < -2.0:
                raise RuntimeError("Unexpected scanner-travel direction")
            row, column = serpentine_coordinates(point_index, source.map_back_and_forth)
            records.append(
                {
                    "curve_row_index": len(records),
                    "acquisition_order": acquisition_order,
                    "instrument_scan_number": metadata["instrument_scan_number"],
                    "block": block,
                    "position_in_block": position_in_block,
                    "source": path.relative_to(ROOT).as_posix(),
                    "map_start_time": metadata["map_start_time"],
                    "speed_um_per_s": speed,
                    "point_index": point_index,
                    "line_index": point_index // GRID_SIZE + 1,
                    "position_within_line": point_index % GRID_SIZE + 1,
                    "physical_row": row,
                    "physical_column": column,
                    "duration_s": float(curve.duration_s),
                    "sample_count": int(curve.deflection_V.size),
                    "travel_nm": travel_nm,
                    "raw_V": np.asarray(curve.deflection_V, dtype=np.float64),
                    "binned_raw_V": binned_raw(travel_nm, curve.deflection_V),
                }
            )
        inventory.append(
            {
                "acquisition_order": acquisition_order,
                "instrument_scan_number": metadata["instrument_scan_number"],
                "block": block,
                "position_in_block": position_in_block,
                "source": path.relative_to(ROOT).as_posix(),
                "sha256": digest,
                "map_start_time": metadata["map_start_time"],
                "speed_um_per_s": speed,
                "approach_points": int(np.median([curve.deflection_V.size for curve in source.curves])),
                "approach_duration_s": float(np.median([curve.duration_s for curve in source.curves])),
                "baseline_adjust_enabled": metadata["baseline_adjust_enabled"],
                "baseline_adjust_begin_of_line": metadata["baseline_adjust_begin_of_line"],
                "baseline_adjust_deadtime_samples": metadata["baseline_adjust_deadtime_samples"],
                "baseline_adjust_average_samples": metadata["baseline_adjust_average_samples"],
            }
        )
        print(f"[{acquisition_order:02d}/{MAP_COUNT}] {speed:g} µm/s · 64 approaches", flush=True)
    counts = {speed: sum(record["speed_um_per_s"] == speed for record in records) for speed in SPEEDS}
    if counts != {speed: 512 for speed in SPEEDS}:
        raise RuntimeError(f"Unexpected per-speed curve counts: {counts}")
    return records, inventory


def build_candidate_table(records: list[dict]) -> tuple[list[dict], dict[float, float]]:
    curve_rows: list[dict] = []
    for record in records:
        for start_nm in CANDIDATE_STARTS_NM:
            curve_rows.append(
                {
                    "acquisition_order": record["acquisition_order"],
                    "block": record["block"],
                    "speed_um_per_s": record["speed_um_per_s"],
                    "candidate_start_nm": float(start_nm),
                    **candidate_components(record["binned_raw_V"], float(start_nm)),
                }
            )
    output: list[dict] = []
    selected: dict[float, float] = {}
    fields = (
        "score_V",
        "drift_V",
        "pre_guard_mismatch_V",
        "post_guard_mismatch_V",
        "roughness_V",
    )
    for speed in SPEEDS:
        candidates: list[dict] = []
        for start_nm in CANDIDATE_STARTS_NM:
            subset = [
                row
                for row in curve_rows
                if row["speed_um_per_s"] == speed and row["candidate_start_nm"] == start_nm
            ]
            result: dict = {
                "speed_um_per_s": speed,
                "candidate_start_nm": float(start_nm),
                "candidate_stop_nm": float(start_nm + WINDOW_WIDTH_NM),
                "candidate_start_time_s": float(start_nm / (speed * 1000.0)),
                "candidate_stop_time_s": float((start_nm + WINDOW_WIDTH_NM) / (speed * 1000.0)),
                "training_blocks": "+".join(map(str, TRAIN_BLOCKS)),
                "validation_blocks": "+".join(map(str, VALIDATION_BLOCKS)),
            }
            for field in fields:
                result[f"training_map_median_{field}"] = aggregate_map_medians(subset, field, TRAIN_BLOCKS)
                result[f"validation_map_median_{field}"] = aggregate_map_medians(
                    subset, field, VALIDATION_BLOCKS
                )
            candidates.append(result)
        winner, flat_bottom = select_stable_candidate(
            candidates, "training_map_median_score_V"
        )
        selected[speed] = float(winner["candidate_start_nm"])
        for row in candidates:
            row["selected_from_training"] = row is winner
            row["within_1pct_contiguous_training_minimum"] = (
                float(row["candidate_start_nm"]) in flat_bottom
            )
            row["fixed_100_220nm_reference"] = bool(
                np.isclose(row["candidate_start_nm"], FIXED_REFERENCE_START_NM)
            )
            output.append(row)
    return output, selected


def build_window_stability(records: list[dict]) -> list[dict]:
    """Leave one acquisition block out and reselect the window."""

    rows: list[dict] = []
    for speed in SPEEDS:
        speed_records = [record for record in records if record["speed_um_per_s"] == speed]
        for held_out in range(1, 5):
            training_blocks = tuple(block for block in range(1, 5) if block != held_out)
            candidates: list[dict] = []
            for start_nm in CANDIDATE_STARTS_NM:
                per_curve = []
                for record in speed_records:
                    components = candidate_components(record["binned_raw_V"], float(start_nm))
                    per_curve.append(
                        {
                            "acquisition_order": record["acquisition_order"],
                            "block": record["block"],
                            "score_V": components["score_V"],
                        }
                    )
                candidates.append(
                    {
                        "candidate_start_nm": float(start_nm),
                        "training_map_median_score_V": aggregate_map_medians(
                            per_curve, "score_V", training_blocks
                        ),
                    }
                )
            winner, flat_bottom = select_stable_candidate(
                candidates, "training_map_median_score_V"
            )
            score = float(winner["training_map_median_score_V"])
            start_nm = float(winner["candidate_start_nm"])
            rows.append(
                {
                    "speed_um_per_s": speed,
                    "held_out_block": held_out,
                    "training_blocks": "+".join(map(str, training_blocks)),
                    "selected_start_nm": start_nm,
                    "selected_stop_nm": start_nm + WINDOW_WIDTH_NM,
                    "training_map_median_score_V": score,
                    "contiguous_near_optimal_start_min_nm": min(flat_bottom),
                    "contiguous_near_optimal_start_max_nm": max(flat_bottom),
                }
            )
    return rows


def apply_baseline(
    records: list[dict], selected: dict[float, float]
) -> tuple[list[dict], np.ndarray, np.ndarray]:
    rows: list[dict] = []
    corrected_bins: list[np.ndarray] = []
    endpoint_bins: list[np.ndarray] = []
    for record in records:
        speed = float(record["speed_um_per_s"])
        start_nm = selected[speed]
        stop_nm = start_nm + WINDOW_WIDTH_NM
        x = record["travel_nm"]
        raw = record["raw_V"]
        baseline_values = window_values(x, raw, start_nm, stop_nm)
        baseline_v = float(np.median(baseline_values))
        corrected = raw - baseline_v
        components = candidate_components(record["binned_raw_V"], start_nm)
        initial_v = float(np.median(window_values(x, raw, 0.0, 10.0)))
        corrected_binned = record["binned_raw_V"] - baseline_v
        corrected_bins.append(corrected_binned)
        endpoint_bins.append(record["binned_raw_V"] - initial_v)
        pre = float(np.median(window_values(x, corrected, start_nm - GUARD_WIDTH_NM, start_nm)))
        post = float(np.median(window_values(x, corrected, stop_nm, stop_nm + GUARD_WIDTH_NM)))
        near_320_360 = float(np.median(window_values(x, corrected, 320.0, 360.0)))
        near_380_420 = float(np.median(window_values(x, corrected, 380.0, 420.0)))
        terminal = float(np.median(corrected[x >= float(np.max(x) - 25.0)]))
        serializable = {
            key: value
            for key, value in record.items()
            if key not in {"travel_nm", "raw_V", "binned_raw_V"}
        }
        rows.append(
            {
                **serializable,
                "baseline_mode": "speed_conditioned_settled_window_median",
                "baseline_start_nm": start_nm,
                "baseline_stop_nm": stop_nm,
                "baseline_start_time_s": start_nm / (speed * 1000.0),
                "baseline_stop_time_s": stop_nm / (speed * 1000.0),
                "baseline_sample_count": int(baseline_values.size),
                "dynamic_baseline_raw_V": baseline_v,
                "initial_0_10nm_raw_V": initial_v,
                "baseline_minus_initial_V": baseline_v - initial_v,
                "pre_guard_corrected_V": pre,
                "post_guard_corrected_V": post,
                "near_320_360nm_corrected_V": near_320_360,
                "near_380_420nm_corrected_V": near_380_420,
                "terminal_last_25nm_corrected_V": terminal,
                **{f"detector_{key}": value for key, value in components.items()},
            }
        )
    return rows, np.vstack(corrected_bins), np.vstack(endpoint_bins)


def build_map_summary(curve_rows: list[dict]) -> list[dict]:
    output: list[dict] = []
    fields = (
        "dynamic_baseline_raw_V",
        "baseline_minus_initial_V",
        "pre_guard_corrected_V",
        "post_guard_corrected_V",
        "near_320_360nm_corrected_V",
        "near_380_420nm_corrected_V",
        "terminal_last_25nm_corrected_V",
        "detector_drift_V",
        "detector_roughness_V",
        "detector_score_V",
    )
    for order in range(1, MAP_COUNT + 1):
        subset = [row for row in curve_rows if int(row["acquisition_order"]) == order]
        first = subset[0]
        result: dict = {
            "acquisition_order": order,
            "instrument_scan_number": first["instrument_scan_number"],
            "block": first["block"],
            "position_in_block": first["position_in_block"],
            "source": first["source"],
            "map_start_time": first["map_start_time"],
            "speed_um_per_s": first["speed_um_per_s"],
            "baseline_start_nm": first["baseline_start_nm"],
            "baseline_stop_nm": first["baseline_stop_nm"],
        }
        for field in fields:
            q25, q50, q75 = quantiles(float(row[field]) for row in subset)
            result[f"{field}_q25"] = q25
            result[f"{field}_median"] = q50
            result[f"{field}_q75"] = q75
            result[f"{field}_spatial_IQR"] = q75 - q25
        output.append(result)
    return output


def build_speed_summary(
    map_rows: list[dict], candidate_rows: list[dict], stability_rows: list[dict]
) -> list[dict]:
    output: list[dict] = []
    for speed in SPEEDS:
        maps = [row for row in map_rows if row["speed_um_per_s"] == speed]
        selected = next(
            row
            for row in candidate_rows
            if row["speed_um_per_s"] == speed and row["selected_from_training"]
        )
        fixed = next(
            row
            for row in candidate_rows
            if row["speed_um_per_s"] == speed and row["fixed_100_220nm_reference"]
        )
        selected_validation = float(selected["validation_map_median_score_V"])
        fixed_validation = float(fixed["validation_map_median_score_V"])
        lobo = [float(row["selected_start_nm"]) for row in stability_rows if row["speed_um_per_s"] == speed]
        result: dict = {
            "speed_um_per_s": speed,
            "map_count": len(maps),
            "selected_start_nm": selected["candidate_start_nm"],
            "selected_stop_nm": selected["candidate_stop_nm"],
            "selected_start_time_s": selected["candidate_start_time_s"],
            "selected_stop_time_s": selected["candidate_stop_time_s"],
            "training_map_median_score_V": selected["training_map_median_score_V"],
            "validation_map_median_score_V": selected_validation,
            "fixed_100_220nm_validation_map_median_score_V": fixed_validation,
            "validation_score_improvement_vs_fixed_fraction": (
                (fixed_validation - selected_validation) / fixed_validation
                if fixed_validation > 0
                else float("nan")
            ),
            "leave_one_block_out_start_nm_min": min(lobo),
            "leave_one_block_out_start_nm_median": median(lobo),
            "leave_one_block_out_start_nm_max": max(lobo),
        }
        for field in (
            "dynamic_baseline_raw_V_median",
            "baseline_minus_initial_V_median",
            "pre_guard_corrected_V_median",
            "post_guard_corrected_V_median",
            "near_320_360nm_corrected_V_median",
            "near_380_420nm_corrected_V_median",
            "terminal_last_25nm_corrected_V_median",
            "detector_drift_V_median",
            "detector_roughness_V_median",
            "detector_score_V_median",
        ):
            values = [float(row[field]) for row in maps]
            q25, q50, q75 = quantiles(values)
            result[f"map_{field}_q25"] = q25
            result[f"map_{field}_median"] = q50
            result[f"map_{field}_q75"] = q75
        output.append(result)
    return output


def build_binned_summary(records: list[dict], corrected: np.ndarray, endpoint: np.ndarray) -> list[dict]:
    rows: list[dict] = []
    speeds = np.asarray([record["speed_um_per_s"] for record in records], dtype=float)
    for speed in SPEEDS:
        selected = speeds == speed
        for index, centre in enumerate(BIN_CENTERS_NM):
            corrected_q = quantiles(corrected[selected, index])
            endpoint_q = quantiles(endpoint[selected, index])
            rows.append(
                {
                    "speed_um_per_s": speed,
                    "scanner_travel_bin_center_nm": centre,
                    "curve_count": int(np.count_nonzero(np.isfinite(corrected[selected, index]))),
                    "endpoint_delta_V_q25": endpoint_q[0],
                    "endpoint_delta_V_median": endpoint_q[1],
                    "endpoint_delta_V_q75": endpoint_q[2],
                    "dynamic_baseline_delta_V_q25": corrected_q[0],
                    "dynamic_baseline_delta_V_median": corrected_q[1],
                    "dynamic_baseline_delta_V_q75": corrected_q[2],
                    "dynamic_baseline_delta_nN_prior_D3_scale_median": corrected_q[1]
                    * D3_SCALE_NN_PER_V,
                }
            )
    return rows


def build_line_position_summary(curve_rows: list[dict]) -> list[dict]:
    output: list[dict] = []
    for speed in SPEEDS:
        for position in range(1, GRID_SIZE + 1):
            selected = [
                row
                for row in curve_rows
                if row["speed_um_per_s"] == speed and row["position_within_line"] == position
            ]
            result: dict = {
                "speed_um_per_s": speed,
                "position_within_line": position,
                "curve_count": len(selected),
            }
            for field in (
                "pre_guard_corrected_V",
                "post_guard_corrected_V",
                "baseline_minus_initial_V",
            ):
                q25, q50, q75 = quantiles(float(row[field]) for row in selected)
                result[f"{field}_q25"] = q25
                result[f"{field}_median"] = q50
                result[f"{field}_q75"] = q75
            output.append(result)
    return output


def plot_candidate_optimization(candidate_rows: list[dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0), sharex=True)
    for ax, speed in zip(axes.flat, SPEEDS, strict=True):
        rows = [row for row in candidate_rows if row["speed_um_per_s"] == speed]
        x = np.asarray([row["candidate_start_nm"] for row in rows])
        train = np.asarray([row["training_map_median_score_V"] for row in rows]) * 1e3
        validation = np.asarray([row["validation_map_median_score_V"] for row in rows]) * 1e3
        winner = next(row for row in rows if row["selected_from_training"])
        ax.plot(x, train, color=COLORS[speed], linewidth=2.0, label="training blocks 1+3")
        ax.plot(x, validation, color="#444444", linewidth=1.5, linestyle="--", label="held-out blocks 2+4")
        ax.axvline(
            float(winner["candidate_start_nm"]),
            color=COLORS[speed],
            linewidth=1.2,
            alpha=0.8,
            label="selected start" if speed == SPEEDS[0] else None,
        )
        ax.axvline(
            FIXED_REFERENCE_START_NM,
            color="#999999",
            linewidth=1.0,
            linestyle=":",
            label="fixed start 100 nm" if speed == SPEEDS[0] else None,
        )
        ax.set_title(f"{speed:g} µm/s · selected {winner['candidate_start_nm']:g}–{winner['candidate_stop_nm']:g} nm")
        ax.set_ylabel("detector score (mV)")
        ax.grid(alpha=0.22)
    for ax in axes[-1, :]:
        ax.set_xlabel("candidate window start (nm scanner travel)")
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle("Speed-conditioned baseline-window selection\n120 nm windows; selection uses training maps only", fontsize=14)
    fig.tight_layout()
    fig.savefig(FIG / "baseline_candidate_optimization.png", dpi=190)
    plt.close(fig)


def plot_detection_effect(binned_rows: list[dict], speed_rows: list[dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0), sharex=True)
    for ax, speed in zip(axes.flat, SPEEDS, strict=True):
        rows = [row for row in binned_rows if row["speed_um_per_s"] == speed]
        x = np.asarray([row["scanner_travel_bin_center_nm"] for row in rows])
        endpoint = np.asarray([row["endpoint_delta_V_median"] for row in rows]) * D3_SCALE_NN_PER_V
        corrected = np.asarray([row["dynamic_baseline_delta_V_median"] for row in rows]) * D3_SCALE_NN_PER_V
        q25 = np.asarray([row["dynamic_baseline_delta_V_q25"] for row in rows]) * D3_SCALE_NN_PER_V
        q75 = np.asarray([row["dynamic_baseline_delta_V_q75"] for row in rows]) * D3_SCALE_NN_PER_V
        summary = next(row for row in speed_rows if row["speed_um_per_s"] == speed)
        ax.axvspan(summary["selected_start_nm"], summary["selected_stop_nm"], color=COLORS[speed], alpha=0.10)
        ax.plot(x, endpoint, color="#777777", linewidth=1.4, label="0–10 nm endpoint reference")
        ax.fill_between(x, q25, q75, color=COLORS[speed], alpha=0.18, linewidth=0)
        ax.plot(x, corrected, color=COLORS[speed], linewidth=2.0, label="dynamic-baseline corrected")
        ax.axhline(0.0, color="black", linewidth=0.7, alpha=0.5)
        ax.set_title(f"{speed:g} µm/s")
        ax.set_ylabel("Δ signal (nN, prior D3 scale)")
        ax.grid(alpha=0.20)
    for ax in axes[-1, :]:
        ax.set_xlabel("scanner travel from approach start (nm)")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("Raw start-up reference versus in-motion dynamic baseline\nshading marks the selected baseline window", fontsize=14)
    fig.tight_layout()
    fig.savefig(FIG / "adaptive_baseline_detection.png", dpi=190)
    plt.close(fig)


def plot_validation(
    binned_rows: list[dict], speed_rows: list[dict], map_rows: list[dict], line_rows: list[dict]
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2))
    ax = axes[0, 0]
    for speed in SPEEDS:
        rows = [row for row in binned_rows if row["speed_um_per_s"] == speed]
        x = np.asarray([row["scanner_travel_bin_center_nm"] for row in rows])
        y = np.asarray([row["dynamic_baseline_delta_V_median"] for row in rows]) * D3_SCALE_NN_PER_V
        ax.plot(x, y, color=COLORS[speed], linewidth=1.8, label=f"{speed:g} µm/s")
    ax.axhline(0.0, color="black", linewidth=0.7, alpha=0.5)
    ax.set_xlim(0, 430)
    ax.set_xlabel("scanner travel (nm)")
    ax.set_ylabel("corrected median (nN, prior D3 scale)")
    ax.set_title("Increment above the dynamic plateau")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.20)

    ax = axes[0, 1]
    x = np.arange(len(SPEEDS))
    selected = np.asarray([next(row for row in speed_rows if row["speed_um_per_s"] == speed)["validation_map_median_score_V"] for speed in SPEEDS]) * 1e3
    fixed = np.asarray([next(row for row in speed_rows if row["speed_um_per_s"] == speed)["fixed_100_220nm_validation_map_median_score_V"] for speed in SPEEDS]) * 1e3
    ax.bar(x - 0.18, fixed, width=0.36, color="#bdbdbd", label="fixed 100–220 nm")
    ax.bar(x + 0.18, selected, width=0.36, color=[COLORS[speed] for speed in SPEEDS], label="optimized")
    ax.set_xticks(x, [f"{speed:g}" for speed in SPEEDS])
    ax.set_xlabel("approach speed (µm/s)")
    ax.set_ylabel("held-out detector score (mV)")
    ax.set_title("Held-out map-block comparison")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="y", alpha=0.20)

    ax = axes[1, 0]
    for speed in SPEEDS:
        rows = [row for row in map_rows if row["speed_um_per_s"] == speed]
        ax.plot(
            [row["acquisition_order"] for row in rows],
            [row["dynamic_baseline_raw_V_median"] for row in rows],
            marker="o",
            color=COLORS[speed],
            linewidth=1.3,
            label=f"{speed:g} µm/s",
        )
    ax.set_xlabel("acquisition order")
    ax.set_ylabel("map-median dynamic baseline (raw V)")
    ax.set_title("Absolute detector level remains history dependent")
    ax.grid(alpha=0.20)

    ax = axes[1, 1]
    for speed in SPEEDS:
        rows = [row for row in line_rows if row["speed_um_per_s"] == speed]
        ax.plot(
            [row["position_within_line"] for row in rows],
            [row["post_guard_corrected_V_median"] * D3_SCALE_NN_PER_V for row in rows],
            marker="o",
            color=COLORS[speed],
            linewidth=1.3,
            label=f"{speed:g} µm/s",
        )
    ax.axhline(0.0, color="black", linewidth=0.7, alpha=0.5)
    ax.set_xlabel("position within scan line")
    ax.set_ylabel("post-window guard residual (nN, prior D3 scale)")
    ax.set_title("Adjacent-window residual by line position")
    ax.grid(alpha=0.20)

    fig.suptitle("Adaptive baseline validation and retained limitations", fontsize=14)
    fig.tight_layout()
    fig.savefig(FIG / "adaptive_baseline_validation.png", dpi=190)
    plt.close(fig)


def build_report(speed_rows: list[dict], stability_rows: list[dict]) -> str:
    lines = [
        "# 12-09-26 · 99.7 wt% glycerol · D3：速度条件化动态 baseline",
        "",
        "## 物理图像",
        "",
        "approach 开始后，cantilever/液体/控制回路先经历速度相关瞬态，随后进入一段近似稳定的恒速运动平台，最后才出现靠近表面的快速变化。因此这里的 baseline 不是开头的静止端点，也不是绝对零力；它被定义为**启动瞬态结束后、近表面上升前的 in-motion 动态平台**。减去它以后得到的是相对该运动平台的增量信号。",
        "",
        "## 选中的 baseline 模式",
        "",
        "- 所有窗口都用 scanner travel 的 nm 定义，绝不使用固定点数；这避免 0.1–2.7 µm/s 不同采样长度造成不同物理范围。",
        "- 每个候选窗口宽 120 nm，并带前后各 20 nm guard。检测分数由窗口首末四分位差、两侧 guard 与窗口中位数的不连续量、5 nm bin 相邻粗糙度组成。",
        "- blocks 1+3 只负责选择，blocks 2+4 完全留出用于评价；分数先在每张 8×8 map 内取中位数，再跨 map 取中位数，避免把 64 pixels 当成独立重复。为避免追逐 5 nm 网格上的单点噪声，在包含全局最小值、且不高于最小值 1% 的连续平底区间中选择中点。",
        "- 窗口一旦按速度选定，每条曲线只减去该窗口 raw vDeflection 的 median；不做线性斜率扣除，避免把可能的距离依赖物理信号拟合掉。",
        "",
        "| approach speed (µm/s) | selected window (nm) | time after motion start (s) | train score (mV) | held-out score (mV) | held-out vs fixed 100–220 nm | LOBO start range (nm) |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in speed_rows:
        improvement = 100.0 * float(row["validation_score_improvement_vs_fixed_fraction"])
        lines.append(
            f"| {row['speed_um_per_s']:g} | {row['selected_start_nm']:g}–{row['selected_stop_nm']:g} | "
            f"{row['selected_start_time_s']:.4f}–{row['selected_stop_time_s']:.4f} | "
            f"{1e3 * row['training_map_median_score_V']:.3f} | "
            f"{1e3 * row['validation_map_median_score_V']:.3f} | {improvement:+.1f}% | "
            f"{row['leave_one_block_out_start_nm_min']:g}–{row['leave_one_block_out_start_nm_max']:g} |"
        )
    lines.extend(
        [
            "",
            "这里的 detector score 只用于比较窗口稳定性，单位虽为 V，但不是误差条、likelihood 或置信区间。LOBO 是每次留出一个完整 acquisition block 后重新选择窗口；它检查窗口选择是否依赖某一 block。",
            "",
            "## 应用后的观测量",
            "",
            "下表以八张 map 的 map median 为统计单位。nN 只采用仓库中既有 D3 比例尺作辅助读数，raw V 仍是主结果。",
            "",
            "| speed (µm/s) | dynamic baseline − initial endpoint (V) | 320–360 nm corrected (V) | 380–420 nm corrected (V) | terminal 25 nm corrected (V) | within-window drift (V) |",
            "|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in speed_rows:
        lines.append(
            f"| {row['speed_um_per_s']:g} | {row['map_baseline_minus_initial_V_median_median']:.6f} | "
            f"{row['map_near_320_360nm_corrected_V_median_median']:.6f} | "
            f"{row['map_near_380_420nm_corrected_V_median_median']:.6f} | "
            f"{row['map_terminal_last_25nm_corrected_V_median_median']:.6f} | "
            f"{row['map_detector_drift_V_median_median']:.6f} |"
        )
    lines.extend(
        [
            "",
            "baseline window 自身的 corrected median 按定义为零，因此没有把它当作验证证据。真正的检查是：留出 blocks 的 detector score、窗口两侧 guard、完整 block 的 LOBO 选择，以及校正后曲线是否仍保留 300 nm 之后的变化。",
            "",
            "## 结论与边界",
            "",
            "1. **必须按速度延迟 baseline 检测。** 最高速 2.7 µm/s 的稳定窗口显著晚于低速；统一使用开头或固定 100–220 nm 会把启动瞬态混进 baseline。",
            "2. **0.1 µm/s 的精确起点约束较弱。** 该速度的瞬态幅度已接近 detector score 的噪声尺度，LOBO 给出 25–65 nm；35–155 nm 是可复现的操作选择，不应解释为精确的 relaxation length。",
            "3. **采用 constant median 而不是 linear detrend。** 这只移除速度/曲线特异的动态平台常数，不主动消去平台内残余斜率；表中的 within-window drift 就是仍保留的非平坦程度。",
            "4. **校正结果是 incremental response，不是 equilibrium force。** 速度相关平台可能含 bulk drag、cantilever relaxation、photodiode/control transient 及真正的远场物理响应。减去它是一种操作性比较基准，不证明这些成分是纯仪器背景。",
            "5. **retract 未使用此 baseline。** 所有 retract 均为 1 µm/s，且已有强负向、历史相关响应；把 approach 的平台规则套给 retract 会改变目标观测量。",
            "6. **没有做 surface-force 拟合。** 本包不做 contact/separation、hydrodynamic、PB、速度外推、回归或机制识别。",
            "",
            "## 文件",
            "",
            "- `candidate_window_scores.csv`：每个速度、每个 120 nm 候选窗口的 train/held-out map-median 分数及组成。",
            "- `window_selection_stability.csv`：完整 acquisition block 的 leave-one-block-out 选择。",
            "- `curve_dynamic_baselines.csv`：2048 条 approach 的选定窗口、baseline、guard、近表面增量和 detector diagnostics。",
            "- `corrected_curve_bins.npz`：逐曲线 5 nm binned raw-endpoint 与 dynamic-baseline corrected arrays；行号与 curve CSV 对齐。",
            "- `map_dynamic_baseline_summary.csv`、`speed_dynamic_baseline_summary.csv`：map 与 speed 层级汇总。",
            "- `binned_corrected_speed_summary.csv`、`line_position_residuals.csv`：曲线形状与 line-position diagnostics。",
            "- `figures/`：窗口选择、校正效果与留出验证。",
            "- `provenance.json`、`artifact_manifest.sha256`：参数、输入哈希和产物校验。",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    records, inventory = load_curves()
    candidate_rows, selected = build_candidate_table(records)
    stability_rows = build_window_stability(records)
    curve_rows, corrected_bins, endpoint_bins = apply_baseline(records, selected)
    map_rows = build_map_summary(curve_rows)
    speed_rows = build_speed_summary(map_rows, candidate_rows, stability_rows)
    binned_rows = build_binned_summary(records, corrected_bins, endpoint_bins)
    line_rows = build_line_position_summary(curve_rows)

    write_csv(OUT / "input_map_inventory.csv", inventory)
    write_csv(OUT / "candidate_window_scores.csv", candidate_rows)
    write_csv(OUT / "window_selection_stability.csv", stability_rows)
    write_csv(OUT / "curve_dynamic_baselines.csv", curve_rows)
    write_csv(OUT / "map_dynamic_baseline_summary.csv", map_rows)
    write_csv(OUT / "speed_dynamic_baseline_summary.csv", speed_rows)
    write_csv(OUT / "binned_corrected_speed_summary.csv", binned_rows)
    write_csv(OUT / "line_position_residuals.csv", line_rows)
    np.savez_compressed(
        OUT / "corrected_curve_bins.npz",
        scanner_travel_bin_centers_nm=BIN_CENTERS_NM,
        dynamic_baseline_corrected_V=corrected_bins,
        endpoint_0_10nm_referenced_V=endpoint_bins,
        curve_row_index=np.arange(len(records), dtype=np.int32),
        acquisition_order=np.asarray([record["acquisition_order"] for record in records], dtype=np.int16),
        speed_um_per_s=np.asarray([record["speed_um_per_s"] for record in records], dtype=np.float64),
        point_index=np.asarray([record["point_index"] for record in records], dtype=np.int16),
    )

    plot_candidate_optimization(candidate_rows)
    plot_detection_effect(binned_rows, speed_rows)
    plot_validation(binned_rows, speed_rows, map_rows, line_rows)
    (OUT / "REPORT.md").write_text(build_report(speed_rows, stability_rows), encoding="utf-8")

    provenance = {
        "analysis": "speed_conditioned_dynamic_baseline_no_surface_force_fit",
        "analysis_script": Path(__file__).relative_to(ROOT).as_posix(),
        "analysis_script_sha256": sha256_file(Path(__file__)),
        "input_download_manifest": DOWNLOAD_MANIFEST.relative_to(ROOT).as_posix(),
        "input_download_manifest_sha256": sha256_file(DOWNLOAD_MANIFEST),
        "source_observational_script": "analysis/analyze_12_09_26_997glycerol_observational.py",
        "source_observational_script_sha256": sha256_file(
            ROOT / "analysis" / "analyze_12_09_26_997glycerol_observational.py"
        ),
        "baseline_definition": "per-curve median raw vDeflection in a speed-specific settled in-motion window",
        "selection": {
            "candidate_start_nm": CANDIDATE_STARTS_NM.tolist(),
            "window_width_nm": WINDOW_WIDTH_NM,
            "guard_width_nm": GUARD_WIDTH_NM,
            "training_blocks": list(TRAIN_BLOCKS),
            "validation_blocks": list(VALIDATION_BLOCKS),
            "aggregation_unit": "pixel score -> map median -> median across maps",
            "score": "abs(first-quarter minus last-quarter median) + 0.5*(pre-guard mismatch + post-guard mismatch) + 0.25*median absolute adjacent-bin difference",
            "stable_selection_rule": "choose the centre of the contiguous candidate-start region containing the training minimum and within 1 percent of that minimum",
            "selected_windows_nm": {
                f"{speed:g}": [selected[speed], selected[speed] + WINDOW_WIDTH_NM]
                for speed in SPEEDS
            },
            "stability_check": "leave one complete acquisition block out and reselect",
        },
        "operations": {
            "binning": f"direct {BIN_WIDTH_NM:g} nm scanner-travel-bin medians; no interpolation",
            "baseline_value": "raw-sample median inside selected physical window",
            "subtraction": "constant only; no slope removal",
            "force_conversion": f"auxiliary only, prior D3 scale {D3_SCALE_NN_PER_V:.12g} nN/V",
            "retract": "not baseline-corrected",
            "fits_performed": [],
            "explicitly_not_performed": [
                "linear baseline slope fit",
                "contact or separation fit",
                "hydrodynamic fit",
                "PB or other surface-force fit",
                "speed regression or zero-speed extrapolation",
                "mechanism identification",
            ],
        },
        "interpretation_boundary": "corrected values are increments relative to a moving dynamic plateau, not absolute or equilibrium force",
        "integrity": {
            "input_map_count": len(inventory),
            "decoded_approach_curves": len(records),
            "input_files_sha256": {
                Path(row["source"]).name: row["sha256"] for row in inventory
            },
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "randomness": "none",
    }
    (OUT / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    artifacts = [
        ROOT / "analysis" / "download_12_09_26_keeper.py",
        ROOT / "analysis" / "analyze_12_09_26_997glycerol_observational.py",
        Path(__file__),
        *(path for path in OUT.rglob("*") if path.is_file() and path.name != "artifact_manifest.sha256"),
    ]
    manifest = "\n".join(
        f"{sha256_file(path)}  {path.relative_to(ROOT).as_posix()}" for path in sorted(artifacts)
    )
    (OUT / "artifact_manifest.sha256").write_text(manifest + "\n", encoding="utf-8")
    print(f"Wrote adaptive-baseline package to {OUT}", flush=True)


if __name__ == "__main__":
    main()

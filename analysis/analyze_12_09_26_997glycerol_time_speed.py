#!/usr/bin/env python3
"""Descriptive time/speed analysis for the 99.7 wt% glycerol D3 maps.

This consumes the speed-conditioned dynamic-baseline result and separates:

1. global map-level rank association with speed and elapsed time;
2. same-speed rank association with elapsed time;
3. exact later-minus-earlier differences within each palindrome pair; and
4. speed ranking of early/late pair means within each acquisition block.

No regression, partial correlation, contact/separation reconstruction, or
surface-force model is fitted.  Map, not pixel, is the statistical unit.
"""

from __future__ import annotations

import csv
from datetime import datetime
import hashlib
import json
from pathlib import Path
import platform
import sys
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
ADAPTIVE = ROOT / "analysis" / "glycerol_99p7_D3_adaptive_baseline_results"
OBSERVATIONAL = ROOT / "analysis" / "glycerol_99p7_D3_observational_results"
OUT = ROOT / "analysis" / "glycerol_99p7_D3_time_speed_results"
FIG = OUT / "figures"

SPEEDS = (0.1, 0.3, 0.9, 2.7)
COLORS = {0.1: "#2a9d8f", 0.3: "#457b9d", 0.9: "#e9c46a", 2.7: "#e76f51"}
BLOCK_COLORS = {1: "#264653", 2: "#2a9d8f", 3: "#e9c46a", 4: "#e76f51"}
D3_SCALE_NN_PER_V = 21.82557616516254

METRICS = {
    "baseline_transient_V": {
        "label": "Dynamic baseline − initial endpoint",
        "unit": "V",
        "source": "adaptive",
        "value": "baseline_minus_initial_V_median",
        "iqr": "baseline_minus_initial_V_spatial_IQR",
    },
    "dynamic_baseline_raw_V": {
        "label": "Absolute dynamic-baseline detector level",
        "unit": "V",
        "source": "adaptive",
        "value": "dynamic_baseline_raw_V_median",
        "iqr": "dynamic_baseline_raw_V_spatial_IQR",
    },
    "corrected_320_360nm_V": {
        "label": "Corrected increment at 320–360 nm travel",
        "unit": "V",
        "source": "adaptive",
        "value": "near_320_360nm_corrected_V_median",
        "iqr": "near_320_360nm_corrected_V_spatial_IQR",
    },
    "corrected_380_420nm_V": {
        "label": "Corrected increment at 380–420 nm travel",
        "unit": "V",
        "source": "adaptive",
        "value": "near_380_420nm_corrected_V_median",
        "iqr": "near_380_420nm_corrected_V_spatial_IQR",
    },
    "corrected_terminal_V": {
        "label": "Corrected terminal last 25 nm",
        "unit": "V",
        "source": "adaptive",
        "value": "terminal_last_25nm_corrected_V_median",
        "iqr": "terminal_last_25nm_corrected_V_spatial_IQR",
    },
    "baseline_window_drift_V": {
        "label": "Within-baseline-window drift magnitude",
        "unit": "V",
        "source": "adaptive",
        "value": "detector_drift_V_median",
        "iqr": "detector_drift_V_spatial_IQR",
    },
    "retract_minimum_delta_V": {
        "label": "Retract minimum, endpoint referenced",
        "unit": "V",
        "source": "observational",
        "value": "retract_minimum_delta_V_median",
        "iqr": "retract_minimum_delta_V_spatial_IQR",
    },
    "retract_minimum_travel_nm": {
        "label": "Travel position of retract minimum",
        "unit": "nm",
        "source": "observational",
        "value": "retract_minimum_scanner_travel_nm_median",
        "iqr": "retract_minimum_scanner_travel_nm_spatial_IQR",
    },
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_manifest(path: Path) -> int:
    count = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        target = ROOT / relative
        if not target.is_file() or sha256_file(target) != expected:
            raise RuntimeError(f"Manifest mismatch: {relative}")
        count += 1
    return count


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


def finite(values: Iterable[float]) -> np.ndarray:
    array = np.asarray(list(values), dtype=np.float64)
    return array[np.isfinite(array)]


def median(values: Iterable[float]) -> float:
    values = finite(values)
    return float(np.median(values)) if values.size else float("nan")


def quantiles(values: Iterable[float]) -> tuple[float, float, float]:
    values = finite(values)
    if not values.size:
        return float("nan"), float("nan"), float("nan")
    return tuple(float(value) for value in np.percentile(values, [25, 50, 75]))


def rank_rho(x: Iterable[float], y: Iterable[float]) -> tuple[float, int]:
    x_array = np.asarray(list(x), dtype=np.float64)
    y_array = np.asarray(list(y), dtype=np.float64)
    keep = np.isfinite(x_array) & np.isfinite(y_array)
    x_array = x_array[keep]
    y_array = y_array[keep]
    if x_array.size < 3 or np.unique(x_array).size < 2 or np.unique(y_array).size < 2:
        return float("nan"), int(x_array.size)
    return float(spearmanr(x_array, y_array).statistic), int(x_array.size)


def midpoint(start: str, stop: str) -> datetime:
    first = datetime.fromisoformat(start)
    last = datetime.fromisoformat(stop)
    return first + (last - first) / 2


def load_map_table() -> list[dict]:
    adaptive = {int(row["acquisition_order"]): row for row in read_csv(ADAPTIVE / "map_dynamic_baseline_summary.csv")}
    observational = {int(row["acquisition_order"]): row for row in read_csv(OBSERVATIONAL / "map_summary.csv")}
    inventory = {int(row["acquisition_order"]): row for row in read_csv(OBSERVATIONAL / "map_inventory.csv")}
    if set(adaptive) != set(range(1, 33)) or set(observational) != set(adaptive) or set(inventory) != set(adaptive):
        raise RuntimeError("Map-level input tables are not the same complete 32-map cohort")
    midpoints = [midpoint(inventory[order]["map_start_time"], inventory[order]["map_end_time"]) for order in range(1, 33)]
    first_midpoint = midpoints[0]
    output: list[dict] = []
    for order in range(1, 33):
        a = adaptive[order]
        o = observational[order]
        i = inventory[order]
        speed = float(a["speed_um_per_s"])
        if not np.isclose(speed, float(o["approach_speed_um_per_s"]), atol=1e-12) or not np.isclose(
            speed, float(i["approach_nominal_speed_um_per_s"]), atol=1e-12
        ):
            raise RuntimeError(f"Speed join mismatch for map {order}")
        if a["source"] != o["source"] or a["source"] != i["source"]:
            raise RuntimeError(f"Source join mismatch for map {order}")
        current_midpoint = midpoints[order - 1]
        row: dict = {
            "acquisition_order": order,
            "instrument_scan_number": int(a["instrument_scan_number"]),
            "block": int(a["block"]),
            "position_in_block": int(a["position_in_block"]),
            "palindrome_mate_order": int(o["palindrome_mate_order"]),
            "source": a["source"],
            "map_start_time": i["map_start_time"],
            "map_end_time": i["map_end_time"],
            "map_midpoint_time": current_midpoint.isoformat(),
            "elapsed_midpoint_min": (current_midpoint - first_midpoint).total_seconds() / 60.0,
            "map_duration_min": (datetime.fromisoformat(i["map_end_time"]) - datetime.fromisoformat(i["map_start_time"])).total_seconds() / 60.0,
            "speed_um_per_s": speed,
        }
        for key, spec in METRICS.items():
            source = a if spec["source"] == "adaptive" else o
            row[key] = float(source[spec["value"]])
            row[f"{key}_spatial_IQR"] = float(source[spec["iqr"]])
        output.append(row)
    if {speed: sum(row["speed_um_per_s"] == speed for row in output) for speed in SPEEDS} != {speed: 8 for speed in SPEEDS}:
        raise RuntimeError("Expected eight maps per approach speed")
    return output


def build_design_balance(map_rows: list[dict]) -> list[dict]:
    rows: list[dict] = []
    rho, n = rank_rho(
        (row["speed_um_per_s"] for row in map_rows),
        (row["elapsed_midpoint_min"] for row in map_rows),
    )
    rows.append(
        {
            "scope": "all_maps",
            "block": "",
            "speed_um_per_s": "",
            "n_maps": n,
            "spearman_rho_speed_vs_elapsed_time": rho,
            "elapsed_first_midpoint_min": min(row["elapsed_midpoint_min"] for row in map_rows),
            "elapsed_last_midpoint_min": max(row["elapsed_midpoint_min"] for row in map_rows),
        }
    )
    for block in range(1, 5):
        selected = [row for row in map_rows if row["block"] == block]
        rho, n = rank_rho(
            (row["speed_um_per_s"] for row in selected),
            (row["elapsed_midpoint_min"] for row in selected),
        )
        rows.append(
            {
                "scope": "within_block_raw_map_order",
                "block": block,
                "speed_um_per_s": "",
                "n_maps": n,
                "spearman_rho_speed_vs_elapsed_time": rho,
                "elapsed_first_midpoint_min": min(row["elapsed_midpoint_min"] for row in selected),
                "elapsed_last_midpoint_min": max(row["elapsed_midpoint_min"] for row in selected),
            }
        )
    for speed in SPEEDS:
        selected = [row for row in map_rows if row["speed_um_per_s"] == speed]
        rows.append(
            {
                "scope": "same_speed_time_coverage",
                "block": "",
                "speed_um_per_s": speed,
                "n_maps": len(selected),
                "spearman_rho_speed_vs_elapsed_time": "",
                "elapsed_first_midpoint_min": min(row["elapsed_midpoint_min"] for row in selected),
                "elapsed_last_midpoint_min": max(row["elapsed_midpoint_min"] for row in selected),
            }
        )
    return rows


def build_palindrome_pairs(map_rows: list[dict]) -> list[dict]:
    output: list[dict] = []
    for metric, spec in METRICS.items():
        for block in range(1, 5):
            block_rows = [row for row in map_rows if row["block"] == block]
            block_midpoint = median(row["elapsed_midpoint_min"] for row in block_rows)
            for speed in SPEEDS:
                pair = sorted(
                    [row for row in block_rows if row["speed_um_per_s"] == speed],
                    key=lambda row: row["position_in_block"],
                )
                if len(pair) != 2 or pair[0]["position_in_block"] + pair[1]["position_in_block"] != 9:
                    raise RuntimeError(f"Palindrome identity failed at block {block}, speed {speed}")
                early, late = pair
                output.append(
                    {
                        "metric": metric,
                        "metric_label": spec["label"],
                        "unit": spec["unit"],
                        "block": block,
                        "speed_um_per_s": speed,
                        "early_order": early["acquisition_order"],
                        "late_order": late["acquisition_order"],
                        "early_elapsed_midpoint_min": early["elapsed_midpoint_min"],
                        "late_elapsed_midpoint_min": late["elapsed_midpoint_min"],
                        "pair_time_gap_min": late["elapsed_midpoint_min"] - early["elapsed_midpoint_min"],
                        "pair_centre_elapsed_min": 0.5 * (early["elapsed_midpoint_min"] + late["elapsed_midpoint_min"]),
                        "pair_centre_offset_from_block_median_min": 0.5
                        * (early["elapsed_midpoint_min"] + late["elapsed_midpoint_min"])
                        - block_midpoint,
                        "early_value": early[metric],
                        "late_value": late[metric],
                        "later_minus_earlier": late[metric] - early[metric],
                        "palindrome_pair_mean": 0.5 * (early[metric] + late[metric]),
                    }
                )
    return output


def build_rank_tables(
    map_rows: list[dict], pair_rows: list[dict]
) -> tuple[list[dict], list[dict], list[dict]]:
    associations: list[dict] = []
    summaries: list[dict] = []
    history: list[dict] = []
    for metric, spec in METRICS.items():
        speed_rho, speed_n = rank_rho(
            (row["speed_um_per_s"] for row in map_rows), (row[metric] for row in map_rows)
        )
        time_rho, time_n = rank_rho(
            (row["elapsed_midpoint_min"] for row in map_rows), (row[metric] for row in map_rows)
        )
        associations.extend(
            [
                {
                    "metric": metric,
                    "metric_label": spec["label"],
                    "unit": spec["unit"],
                    "scope": "all_32_maps",
                    "predictor": "speed_um_per_s",
                    "speed_stratum_um_per_s": "",
                    "block": "",
                    "n_maps_or_pair_means": speed_n,
                    "spearman_rho": speed_rho,
                    "inference_boundary": "descriptive rank association; repeated sequential maps",
                },
                {
                    "metric": metric,
                    "metric_label": spec["label"],
                    "unit": spec["unit"],
                    "scope": "all_32_maps",
                    "predictor": "elapsed_midpoint_min",
                    "speed_stratum_um_per_s": "",
                    "block": "",
                    "n_maps_or_pair_means": time_n,
                    "spearman_rho": time_rho,
                    "inference_boundary": "descriptive serial association; speed protocol is balanced but not randomized",
                },
            ]
        )
        within_speed: list[float] = []
        for speed in SPEEDS:
            selected = [row for row in map_rows if row["speed_um_per_s"] == speed]
            rho, n = rank_rho(
                (row["elapsed_midpoint_min"] for row in selected), (row[metric] for row in selected)
            )
            within_speed.append(rho)
            associations.append(
                {
                    "metric": metric,
                    "metric_label": spec["label"],
                    "unit": spec["unit"],
                    "scope": "same_speed_across_four_blocks",
                    "predictor": "elapsed_midpoint_min",
                    "speed_stratum_um_per_s": speed,
                    "block": "",
                    "n_maps_or_pair_means": n,
                    "spearman_rho": rho,
                    "inference_boundary": "descriptive same-speed serial association; two maps per block",
                }
            )
        block_speed: list[float] = []
        high_low: list[float] = []
        metric_pairs = [row for row in pair_rows if row["metric"] == metric]
        for block in range(1, 5):
            selected = sorted(
                [row for row in metric_pairs if row["block"] == block],
                key=lambda row: row["speed_um_per_s"],
            )
            rho, n = rank_rho(
                (row["speed_um_per_s"] for row in selected),
                (row["palindrome_pair_mean"] for row in selected),
            )
            block_speed.append(rho)
            high_low.append(selected[-1]["palindrome_pair_mean"] - selected[0]["palindrome_pair_mean"])
            associations.append(
                {
                    "metric": metric,
                    "metric_label": spec["label"],
                    "unit": spec["unit"],
                    "scope": "within_block_palindrome_pair_means",
                    "predictor": "speed_um_per_s",
                    "speed_stratum_um_per_s": "",
                    "block": block,
                    "n_maps_or_pair_means": n,
                    "spearman_rho": rho,
                    "inference_boundary": "four pair means; arithmetic early/late mean cancels only first-order symmetric drift",
                }
            )
        summaries.append(
            {
                "metric": metric,
                "metric_label": spec["label"],
                "unit": spec["unit"],
                "global_rho_speed": speed_rho,
                "global_rho_elapsed_time": time_rho,
                "same_speed_time_rho_median": median(within_speed),
                "same_speed_time_rho_min": float(np.nanmin(within_speed)),
                "same_speed_time_rho_max": float(np.nanmax(within_speed)),
                "block_pair_mean_speed_rho_median": median(block_speed),
                "block_pair_mean_speed_rho_min": float(np.nanmin(block_speed)),
                "block_pair_mean_speed_rho_max": float(np.nanmax(block_speed)),
                "block_pair_mean_speed_rho_positive_fraction": float(np.mean(np.asarray(block_speed) > 0)),
                "pair_mean_2p7_minus_0p1_median": median(high_low),
                "pair_mean_2p7_minus_0p1_min": min(high_low),
                "pair_mean_2p7_minus_0p1_max": max(high_low),
            }
        )
        for speed in SPEEDS:
            selected = [row["later_minus_earlier"] for row in metric_pairs if row["speed_um_per_s"] == speed]
            q25, q50, q75 = quantiles(selected)
            history.append(
                {
                    "metric": metric,
                    "metric_label": spec["label"],
                    "unit": spec["unit"],
                    "speed_um_per_s": speed,
                    "block_pair_count": len(selected),
                    "later_minus_earlier_q25": q25,
                    "later_minus_earlier_median": q50,
                    "later_minus_earlier_q75": q75,
                    "later_minus_earlier_min": min(selected),
                    "later_minus_earlier_max": max(selected),
                    "positive_fraction": float(np.mean(np.asarray(selected) > 0)),
                }
            )
    return associations, summaries, history


def build_map_binned(map_rows: list[dict]) -> tuple[list[dict], np.ndarray, np.ndarray]:
    archive = np.load(ADAPTIVE / "corrected_curve_bins.npz", allow_pickle=False)
    corrected = archive["dynamic_baseline_corrected_V"]
    order = archive["acquisition_order"].astype(int)
    centres = archive["scanner_travel_bin_centers_nm"].astype(float)
    if corrected.shape != (2048, centres.size) or set(order) != set(range(1, 33)):
        raise RuntimeError("Unexpected corrected-curve archive shape or map identity")
    matrix = np.full((32, centres.size), np.nan)
    rows: list[dict] = []
    for map_row in map_rows:
        acquisition_order = int(map_row["acquisition_order"])
        selected = corrected[order == acquisition_order]
        if selected.shape[0] != 64:
            raise RuntimeError(f"Map {acquisition_order} does not have 64 corrected curves")
        matrix[acquisition_order - 1] = np.nanmedian(selected, axis=0)
        counts = np.sum(np.isfinite(selected), axis=0)
        for index, centre in enumerate(centres):
            rows.append(
                {
                    "acquisition_order": acquisition_order,
                    "block": map_row["block"],
                    "position_in_block": map_row["position_in_block"],
                    "speed_um_per_s": map_row["speed_um_per_s"],
                    "elapsed_midpoint_min": map_row["elapsed_midpoint_min"],
                    "scanner_travel_bin_center_nm": centre,
                    "finite_pixel_count": int(counts[index]),
                    "map_median_corrected_V": matrix[acquisition_order - 1, index],
                    "map_median_corrected_nN_prior_D3_scale": matrix[acquisition_order - 1, index]
                    * D3_SCALE_NN_PER_V,
                }
            )
    return rows, matrix, centres


def build_distance_tables(
    map_rows: list[dict], matrix: np.ndarray, centres: np.ndarray
) -> tuple[list[dict], list[dict]]:
    associations: list[dict] = []
    pairs: list[dict] = []
    speed = np.asarray([row["speed_um_per_s"] for row in map_rows])
    elapsed = np.asarray([row["elapsed_midpoint_min"] for row in map_rows])
    blocks = np.asarray([row["block"] for row in map_rows])
    positions = np.asarray([row["position_in_block"] for row in map_rows])
    orders = np.asarray([row["acquisition_order"] for row in map_rows])
    for index, centre in enumerate(centres):
        y = matrix[:, index]
        speed_rho, n = rank_rho(speed, y)
        time_rho, _ = rank_rho(elapsed, y)
        result: dict = {
            "scanner_travel_bin_center_nm": centre,
            "finite_map_count": int(np.count_nonzero(np.isfinite(y))),
            "global_rho_speed": speed_rho,
            "global_rho_elapsed_time": time_rho,
        }
        within_speed: list[float] = []
        for value in SPEEDS:
            selected = speed == value
            rho, _ = rank_rho(elapsed[selected], y[selected])
            within_speed.append(rho)
            result[f"same_speed_{value:g}_rho_time"] = rho
        result["same_speed_rho_time_median"] = median(within_speed)
        result["same_speed_rho_time_min"] = float(np.nanmin(within_speed))
        result["same_speed_rho_time_max"] = float(np.nanmax(within_speed))
        block_rhos: list[float] = []
        high_low: list[float] = []
        for block in range(1, 5):
            pair_means: list[float] = []
            for value in SPEEDS:
                selected = (blocks == block) & (speed == value)
                selected_orders = orders[selected]
                selected_positions = positions[selected]
                selected_values = y[selected]
                if selected_values.size != 2 or int(np.sum(selected_positions)) != 9:
                    raise RuntimeError("Distancewise palindrome pairing failed")
                sort = np.argsort(selected_positions)
                selected_values = selected_values[sort]
                selected_orders = selected_orders[sort]
                selected_positions = selected_positions[sort]
                pair_mean = float(np.mean(selected_values))
                pair_delta = float(selected_values[1] - selected_values[0])
                pair_means.append(pair_mean)
                pairs.append(
                    {
                        "scanner_travel_bin_center_nm": centre,
                        "block": block,
                        "speed_um_per_s": value,
                        "early_order": int(selected_orders[0]),
                        "late_order": int(selected_orders[1]),
                        "early_position_in_block": int(selected_positions[0]),
                        "late_position_in_block": int(selected_positions[1]),
                        "early_map_median_corrected_V": float(selected_values[0]),
                        "late_map_median_corrected_V": float(selected_values[1]),
                        "later_minus_earlier_V": pair_delta,
                        "palindrome_pair_mean_corrected_V": pair_mean,
                    }
                )
            rho, _ = rank_rho(SPEEDS, pair_means)
            block_rhos.append(rho)
            high_low.append(pair_means[-1] - pair_means[0])
            result[f"block_{block}_pair_mean_rho_speed"] = rho
            result[f"block_{block}_pair_mean_2p7_minus_0p1_V"] = high_low[-1]
        result["block_pair_mean_rho_speed_median"] = median(block_rhos)
        result["block_pair_mean_rho_speed_min"] = float(np.nanmin(block_rhos))
        result["block_pair_mean_rho_speed_max"] = float(np.nanmax(block_rhos))
        result["block_pair_mean_2p7_minus_0p1_V_median"] = median(high_low)
        result["block_pair_mean_2p7_minus_0p1_V_min"] = min(high_low)
        result["block_pair_mean_2p7_minus_0p1_V_max"] = max(high_low)
        associations.append(result)
    return associations, pairs


def block_boundaries(map_rows: list[dict]) -> list[float]:
    output: list[float] = []
    for order in (8, 16, 24):
        left = next(row for row in map_rows if row["acquisition_order"] == order)
        right = next(row for row in map_rows if row["acquisition_order"] == order + 1)
        output.append(0.5 * (left["elapsed_midpoint_min"] + right["elapsed_midpoint_min"]))
    return output


def plot_chronology(map_rows: list[dict]) -> None:
    panels = (
        "baseline_transient_V",
        "dynamic_baseline_raw_V",
        "corrected_320_360nm_V",
        "corrected_380_420nm_V",
        "corrected_terminal_V",
        "retract_minimum_delta_V",
    )
    fig, axes = plt.subplots(3, 2, figsize=(14.0, 11.0), sharex=True)
    x = np.asarray([row["elapsed_midpoint_min"] for row in map_rows])
    colors = [COLORS[row["speed_um_per_s"]] for row in map_rows]
    for ax, metric in zip(axes.flat, panels, strict=True):
        spec = METRICS[metric]
        y = np.asarray([row[metric] for row in map_rows])
        spread = np.asarray([row[f"{metric}_spatial_IQR"] for row in map_rows])
        ax.errorbar(x, y, yerr=spread / 2.0, fmt="none", color="0.70", linewidth=0.7, alpha=0.8)
        ax.plot(x, y, color="0.50", linewidth=0.8, alpha=0.9)
        ax.scatter(x, y, c=colors, s=35, edgecolor="white", linewidth=0.35, zorder=3)
        for boundary in block_boundaries(map_rows):
            ax.axvline(boundary, color="0.55", linestyle=":", linewidth=0.9)
        ax.set_title(spec["label"])
        ax.set_ylabel(spec["unit"])
        ax.grid(alpha=0.18)
    for ax in axes[-1, :]:
        ax.set_xlabel("elapsed time at map midpoint (min)")
    handles = [
        Line2D([], [], marker="o", linestyle="none", color=COLORS[speed], label=f"{speed:g} µm/s")
        for speed in SPEEDS
    ]
    fig.legend(handles=handles, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 0.965))
    fig.suptitle(
        "Map-level chronology after speed-conditioned dynamic baseline\n"
        "points are map medians; bars are half spatial IQR; dotted lines separate palindrome blocks",
        y=0.995,
        fontsize=14,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(FIG / "time_speed_chronology.png", dpi=200)
    plt.close(fig)


def plot_distance_associations(distance_rows: list[dict]) -> None:
    x = np.asarray([row["scanner_travel_bin_center_nm"] for row in distance_rows])
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.5), sharex=True)
    axes[0, 0].plot(x, [row["global_rho_speed"] for row in distance_rows], color="#d1495b", lw=2, label="speed")
    axes[0, 0].plot(x, [row["global_rho_elapsed_time"] for row in distance_rows], color="#30638e", lw=2, label="elapsed time")
    axes[0, 0].set_title("Global rank association across 32 maps")
    axes[0, 0].set_ylabel("Spearman ρ")
    axes[0, 0].legend(frameon=False)

    for speed in SPEEDS:
        axes[0, 1].plot(
            x,
            [row[f"same_speed_{speed:g}_rho_time"] for row in distance_rows],
            color=COLORS[speed],
            lw=1.6,
            label=f"{speed:g} µm/s",
        )
    axes[0, 1].plot(
        x,
        [row["same_speed_rho_time_median"] for row in distance_rows],
        color="black",
        lw=2.0,
        linestyle="--",
        label="median",
    )
    axes[0, 1].set_title("Time association within each speed (8 maps)")
    axes[0, 1].set_ylabel("Spearman ρ with elapsed time")
    axes[0, 1].legend(frameon=False, fontsize=8, ncol=2)

    for block in range(1, 5):
        axes[1, 0].plot(
            x,
            [row[f"block_{block}_pair_mean_rho_speed"] for row in distance_rows],
            color=BLOCK_COLORS[block],
            lw=1.5,
            label=f"block {block}",
        )
    axes[1, 0].plot(
        x,
        [row["block_pair_mean_rho_speed_median"] for row in distance_rows],
        color="black",
        lw=2.0,
        linestyle="--",
        label="median",
    )
    axes[1, 0].set_title("Speed rank after palindrome pair averaging")
    axes[1, 0].set_ylabel("Spearman ρ across four speeds")
    axes[1, 0].legend(frameon=False, fontsize=8, ncol=2)

    for block in range(1, 5):
        axes[1, 1].plot(
            x,
            [row[f"block_{block}_pair_mean_2p7_minus_0p1_V"] for row in distance_rows],
            color=BLOCK_COLORS[block],
            lw=1.5,
            label=f"block {block}",
        )
    axes[1, 1].plot(
        x,
        [row["block_pair_mean_2p7_minus_0p1_V_median"] for row in distance_rows],
        color="black",
        lw=2.0,
        linestyle="--",
        label="median",
    )
    axes[1, 1].set_title("Pair-mean high-minus-low speed contrast")
    axes[1, 1].set_ylabel("2.7 − 0.1 µm/s (V)")
    axes[1, 1].legend(frameon=False, fontsize=8, ncol=2)

    for ax in axes.flat:
        ax.axhline(0.0, color="0.4", linewidth=0.7)
        ax.set_xlim(0, 430)
        ax.grid(alpha=0.18)
    for ax in axes[-1, :]:
        ax.set_xlabel("scanner travel from approach start (nm)")
    fig.suptitle(
        "Distance-resolved time/speed associations\n"
        "dynamic-baseline-corrected map medians; rank statistics and exact palindrome contrasts only",
        fontsize=14,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIG / "distance_resolved_time_speed_associations.png", dpi=200)
    plt.close(fig)


def plot_scalar_rank_summary(summary_rows: list[dict]) -> None:
    columns = (
        ("global_rho_speed", "global\nspeed"),
        ("global_rho_elapsed_time", "global\ntime"),
        ("same_speed_time_rho_median", "same-speed\ntime median"),
        ("block_pair_mean_speed_rho_median", "pair-mean\nspeed median"),
    )
    matrix = np.asarray([[float(row[key]) for key, _ in columns] for row in summary_rows])
    labels = [row["metric_label"] for row in summary_rows]
    fig, ax = plt.subplots(figsize=(9.0, 6.2))
    image = ax.imshow(matrix, vmin=-1, vmax=1, cmap="coolwarm", aspect="auto")
    ax.set_xticks(range(len(columns)), [label for _, label in columns])
    ax.set_yticks(range(len(labels)), labels)
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = matrix[row_index, column_index]
            ax.text(
                column_index,
                row_index,
                f"{value:.2f}",
                ha="center",
                va="center",
                color="white" if abs(value) > 0.55 else "black",
                fontsize=9,
            )
    colorbar = fig.colorbar(image, ax=ax, shrink=0.85)
    colorbar.set_label("Spearman ρ")
    ax.set_title("Map-level rank-association summary\nno p-values: maps are sequential and not independent randomized trials")
    fig.tight_layout()
    fig.savefig(FIG / "scalar_rank_association_summary.png", dpi=200)
    plt.close(fig)


def plot_palindrome_separation(pair_rows: list[dict]) -> None:
    metrics = (
        "baseline_transient_V",
        "corrected_320_360nm_V",
        "corrected_380_420nm_V",
        "retract_minimum_delta_V",
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5))
    speed_x = np.arange(len(SPEEDS))
    for ax, metric in zip(axes.flat, metrics, strict=True):
        selected = [row for row in pair_rows if row["metric"] == metric]
        for block in range(1, 5):
            rows = sorted(
                [row for row in selected if row["block"] == block],
                key=lambda row: row["speed_um_per_s"],
            )
            ax.plot(
                speed_x,
                [row["palindrome_pair_mean"] for row in rows],
                color=BLOCK_COLORS[block],
                marker="o",
                linewidth=1.4,
                label=f"block {block}",
            )
        ax.set_xticks(speed_x, [f"{speed:g}" for speed in SPEEDS])
        ax.set_title(METRICS[metric]["label"])
        ax.set_ylabel(METRICS[metric]["unit"])
        ax.grid(alpha=0.18)
    for ax in axes[-1, :]:
        ax.set_xlabel("approach speed (µm/s)")
    axes[0, 0].legend(frameon=False, fontsize=8, ncol=2)
    fig.suptitle(
        "Speed ordering after early/late palindrome pair averaging\n"
        "each point is the arithmetic mean of two same-speed maps symmetric in block order",
        fontsize=14,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIG / "palindrome_pair_mean_speed_profiles.png", dpi=200)
    plt.close(fig)


def fmt(value: float) -> str:
    return f"{value:.3f}" if np.isfinite(value) else "n/a"


def build_report(
    design_rows: list[dict],
    summary_rows: list[dict],
    history_rows: list[dict],
    pair_rows: list[dict],
    distance_rows: list[dict],
) -> str:
    design_global = next(row for row in design_rows if row["scope"] == "all_maps")
    summary = {row["metric"]: row for row in summary_rows}
    pair_gaps = [
        float(row["pair_time_gap_min"])
        for row in pair_rows
        if row["metric"] == "corrected_380_420nm_V"
    ]
    pair_offsets = [
        abs(float(row["pair_centre_offset_from_block_median_min"]))
        for row in pair_rows
        if row["metric"] == "corrected_380_420nm_V"
    ]
    strong_speed_bins = [
        float(row["scanner_travel_bin_center_nm"])
        for row in distance_rows
        if float(row["global_rho_speed"]) >= 0.75
        and float(row["block_pair_mean_rho_speed_median"]) >= 0.75
    ]
    lines = [
        "# 12-09-26 · 99.7 wt% glycerol · D3：时间/速度相关性",
        "",
        "## 物理图像",
        "",
        "动态 baseline 校正后，剩余信号仍可能同时随 approach speed、实验经过时间、surface/probe conditioning 和 raster history 改变。四个 rotated palindrome blocks 让每个速度在每个 block 内出现两次，并位于对称 acquisition position；因此可以分别观察同速度早晚变化，以及先对称平均、再比较速度排序。这个设计能压低一阶时间漂移的混淆，但不能把顺序实验变成随机化速度因果实验。",
        "",
        "## 分析定义",
        "",
        "- 输入是 speed-conditioned dynamic-baseline 包；approach 每条曲线只做已验证的 constant-median baseline subtraction。",
        "- 统计单位是 8×8 map median，不把 64 pixels 当成独立重复。",
        "- 只计算 Spearman rank correlation、same-speed chronology、palindrome later−earlier exact difference 和 early/late arithmetic pair mean；不做 regression、partial correlation、time model 或 surface-force fit。",
        "- 不报告常规 p-value，因为 32 张 map 是串行、相关、非随机化观测；ρ 只描述当前序列的 rank association。",
        "- scanner travel 不是 tip–surface separation；320–360/380–420 nm 是固定 travel window，而不是绝对表面距离。",
        "",
        "## Protocol 的 speed/time 平衡",
        "",
        f"全部 32 maps 中 speed 与实际 map-midpoint elapsed time 的 Spearman ρ = **{design_global['spearman_rho_speed_vs_elapsed_time']:.3f}**。每个 palindrome block 内使用原始 8-map 顺序计算的 speed/time ρ 为："
        + ", ".join(
            f"block {row['block']} {row['spearman_rho_speed_vs_elapsed_time']:.3f}"
            for row in design_rows
            if row["scope"] == "within_block_raw_map_order"
        )
        + "。这说明速度标签与线性顺序基本正交，但非线性 drift、block 差异和速度依赖 map duration 仍然存在。",
        "",
        f"同速度 early/late pair 的实际时间间隔为 {min(pair_gaps):.2f}–{max(pair_gaps):.2f} min；各 pair 的平均时间偏离对应 block midpoint 最多 {max(pair_offsets):.3f} min。因此 arithmetic pair mean 很接近共同 block 时间中心，但 later−earlier difference 的幅值不能解释成统一时间间隔下的变化率。",
        "",
        "## Map-level rank association",
        "",
        "| Observable | global ρ(speed) | global ρ(time) | same-speed ρ(time), median [range] | block pair-mean ρ(speed), median [range] | pair-mean 2.7−0.1, median [range] |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for metric in METRICS:
        row = summary[metric]
        lines.append(
            f"| {row['metric_label']} ({row['unit']}) | {fmt(row['global_rho_speed'])} | "
            f"{fmt(row['global_rho_elapsed_time'])} | {fmt(row['same_speed_time_rho_median'])} "
            f"[{fmt(row['same_speed_time_rho_min'])}, {fmt(row['same_speed_time_rho_max'])}] | "
            f"{fmt(row['block_pair_mean_speed_rho_median'])} "
            f"[{fmt(row['block_pair_mean_speed_rho_min'])}, {fmt(row['block_pair_mean_speed_rho_max'])}] | "
            f"{row['pair_mean_2p7_minus_0p1_median']:.6g} "
            f"[{row['pair_mean_2p7_minus_0p1_min']:.6g}, {row['pair_mean_2p7_minus_0p1_max']:.6g}] |"
        )
    lines.extend(
        [
            "",
            "`global ρ(speed)` 利用全部 32 maps；`same-speed ρ(time)` 是四个速度各自八张 map 的结果再取中位数；`block pair-mean ρ(speed)` 在每个 block 内只剩四个速度点，因此 ρ 很离散，只用于方向/排序诊断。",
            "",
            "## 主要观测结果",
            "",
            f"1. **中距离增量保留一致的正速度排序。** 320–360 nm 的 global ρ(speed) = {summary['corrected_320_360nm_V']['global_rho_speed']:.3f}，四个 block 的 pair-mean ρ(speed) 全为 1；2.7−0.1 µm/s pair-mean 差的 block median 为 {summary['corrected_320_360nm_V']['pair_mean_2p7_minus_0p1_median']:.6f} V。380–420 nm 对应值为 ρ = {summary['corrected_380_420nm_V']['global_rho_speed']:.3f}、四个 block 均为 1、差值 {summary['corrected_380_420nm_V']['pair_mean_2p7_minus_0p1_median']:.6f} V。",
            f"2. **同一速度仍有明显 time/history relaxation。** 320–360 nm 的 same-speed ρ(time) 中位数为 {summary['corrected_320_360nm_V']['same_speed_time_rho_median']:.3f}，范围 {summary['corrected_320_360nm_V']['same_speed_time_rho_min']:.3f}–{summary['corrected_320_360nm_V']['same_speed_time_rho_max']:.3f}；其中 0.1 µm/s 是弱瞬态/噪声主导的例外。380–420 nm 的四个速度均为负，范围 {summary['corrected_380_420nm_V']['same_speed_time_rho_min']:.3f}–{summary['corrected_380_420nm_V']['same_speed_time_rho_max']:.3f}。",
            f"3. **速度排序只在一段 travel 区域内稳定。** global 和 block-pair-median speed ρ 同时不低于 0.75 的连续主区间是约 {min(strong_speed_bins):g}–{max(strong_speed_bins):g} nm；在早期区间看到的是相对动态 baseline 的 startup tail，约 427.5 nm 后则进入 terminal/contact-like 排序翻转。",
            f"4. **terminal response 与中距离相反。** terminal global ρ(speed) = {summary['corrected_terminal_V']['global_rho_speed']:.3f}，四个 block 的 pair-mean ρ 均为 −1；same-speed ρ(time) 中位数却为 {summary['corrected_terminal_V']['same_speed_time_rho_median']:.3f}。这表明 terminal load/contact alignment 与中距离增量不是同一个可直接互换的 observable。",
            f"5. **retract 更突出 history。** retract minimum 的 global ρ(time) = {summary['retract_minimum_delta_V']['global_rho_elapsed_time']:.3f}，same-speed ρ(time) 中位数 = {summary['retract_minimum_delta_V']['same_speed_time_rho_median']:.3f}；later−earlier 多数为正，即负向 minimum 随实验进行通常变得较浅。",
            "",
            "## Same-speed palindrome history",
            "",
            "下表给出每个 block 中同速度 later−earlier 的四个 exact differences 的中位数。若它明显不为零，说明即使速度固定，观测量仍随该 block 内的时间/history 改变。",
            "",
            "| Observable | 0.1 µm/s | 0.3 µm/s | 0.9 µm/s | 2.7 µm/s |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for metric in (
        "baseline_transient_V",
        "corrected_320_360nm_V",
        "corrected_380_420nm_V",
        "corrected_terminal_V",
        "retract_minimum_delta_V",
    ):
        values = {
            float(row["speed_um_per_s"]): row
            for row in history_rows
            if row["metric"] == metric
        }
        lines.append(
            f"| {METRICS[metric]['label']} ({METRICS[metric]['unit']}) | "
            + " | ".join(f"{values[speed]['later_minus_earlier_median']:.6g}" for speed in SPEEDS)
            + " |"
        )
    lines.extend(
        [
            "",
            "## 解释边界",
            "",
            "- 正的 speed rank association 只说明在本 protocol 的 pair-mean map signal 中较快 approach 通常对应较大增量；它不等于已识别的 hydrodynamic coefficient。",
            "- same-speed time association 与 later−earlier difference 是 history/systematic 的直接证据；它们限制任何把跨速度差异全部归因于速度的解释。",
            "- dynamic baseline 本身包含速度相关 motion-start response；绝对 baseline detector level 还会漂移和 reset，因此二者分别报告。",
            "- retract speed 始终为 1 µm/s。retract minimum 对 preceding approach speed 的相关性是接触历史关联，不是 retract velocity dependence。",
            "- terminal last-25-nm response 属于 contact-like 区域且随速度排序可与中距离相反，不能作为相同 separation 的 equilibrium force 比较。",
            "",
            "## 文件",
            "",
            "- `map_time_speed_observables.csv`：32-map chronology 与八个 map-level observables。",
            "- `rank_associations.csv`、`rank_association_summary.csv`：global、same-speed 与 block pair-mean rank associations。",
            "- `palindrome_history_pairs.csv`、`palindrome_history_summary.csv`：同速度 exact later−earlier differences 和 pair means。",
            "- `map_binned_corrected_curves.csv`：每张 map 的 5 nm corrected median curve。",
            "- `distancewise_rank_associations.csv`、`distancewise_palindrome_pairs.csv`：逐 scanner-travel bin 的时间/速度诊断。",
            "- `figures/`：chronology、distance-resolved correlation、scalar summary 与 palindrome speed profiles。",
            "- `provenance.json`、`artifact_manifest.sha256`：输入身份、算法边界和产物哈希。",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    adaptive_manifest_count = verify_manifest(ADAPTIVE / "artifact_manifest.sha256")
    observational_manifest_count = verify_manifest(OBSERVATIONAL / "artifact_manifest.sha256")

    map_rows = load_map_table()
    design_rows = build_design_balance(map_rows)
    pair_rows = build_palindrome_pairs(map_rows)
    association_rows, summary_rows, history_rows = build_rank_tables(map_rows, pair_rows)
    map_binned_rows, matrix, centres = build_map_binned(map_rows)
    distance_rows, distance_pair_rows = build_distance_tables(map_rows, matrix, centres)

    write_csv(OUT / "map_time_speed_observables.csv", map_rows)
    write_csv(OUT / "design_balance.csv", design_rows)
    write_csv(OUT / "rank_associations.csv", association_rows)
    write_csv(OUT / "rank_association_summary.csv", summary_rows)
    write_csv(OUT / "palindrome_history_pairs.csv", pair_rows)
    write_csv(OUT / "palindrome_history_summary.csv", history_rows)
    write_csv(OUT / "map_binned_corrected_curves.csv", map_binned_rows)
    write_csv(OUT / "distancewise_rank_associations.csv", distance_rows)
    write_csv(OUT / "distancewise_palindrome_pairs.csv", distance_pair_rows)

    plot_chronology(map_rows)
    plot_distance_associations(distance_rows)
    plot_scalar_rank_summary(summary_rows)
    plot_palindrome_separation(pair_rows)
    (OUT / "REPORT.md").write_text(
        build_report(
            design_rows,
            summary_rows,
            history_rows,
            pair_rows,
            distance_rows,
        ),
        encoding="utf-8",
    )

    provenance = {
        "analysis": "descriptive_map_level_time_speed_rank_and_palindrome_analysis",
        "analysis_script": Path(__file__).relative_to(ROOT).as_posix(),
        "analysis_script_sha256": sha256_file(Path(__file__)),
        "input_packages": {
            "adaptive_baseline": {
                "path": ADAPTIVE.relative_to(ROOT).as_posix(),
                "manifest_sha256": sha256_file(ADAPTIVE / "artifact_manifest.sha256"),
                "verified_manifest_entries": adaptive_manifest_count,
            },
            "observational": {
                "path": OBSERVATIONAL.relative_to(ROOT).as_posix(),
                "manifest_sha256": sha256_file(OBSERVATIONAL / "artifact_manifest.sha256"),
                "verified_manifest_entries": observational_manifest_count,
            },
        },
        "cohort": {
            "map_count": 32,
            "pixels_per_map": 64,
            "blocks": 4,
            "maps_per_block": 8,
            "speeds_um_per_s": list(SPEEDS),
            "maps_per_speed": 8,
            "palindrome_pairs_per_speed": 4,
        },
        "statistical_unit": "map median",
        "operations": {
            "chronology_time": "actual midpoint of first and last segment timestamps in each map",
            "global_association": "Spearman rank correlation across 32 map medians",
            "same_speed_time_association": "Spearman rank correlation across eight map medians at each speed",
            "palindrome_history": "exact later minus earlier map-median difference for each same-speed symmetric pair",
            "palindrome_speed": "arithmetic mean of each same-speed early/late pair, then Spearman rank across four speeds within a block",
            "distance_resolved_input": "median of 64 per-curve dynamic-baseline-corrected 5 nm bins within each map",
            "p_values": "not reported because maps are sequential correlated observations, not randomized independent trials",
            "fits_performed": [],
            "explicitly_not_performed": [
                "linear or nonlinear time regression",
                "partial correlation",
                "contact or separation reconstruction",
                "hydrodynamic fit",
                "PB or other surface-force fit",
                "zero-speed extrapolation",
            ],
        },
        "interpretation_boundary": "rank associations and exact paired differences are descriptive; they do not identify speed or time causally",
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "matplotlib": matplotlib.__version__,
        },
        "randomness": "none",
    }
    (OUT / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    artifacts = [
        Path(__file__),
        *(path for path in OUT.rglob("*") if path.is_file() and path.name != "artifact_manifest.sha256"),
    ]
    manifest = "\n".join(
        f"{sha256_file(path)}  {path.relative_to(ROOT).as_posix()}" for path in sorted(artifacts)
    )
    (OUT / "artifact_manifest.sha256").write_text(manifest + "\n", encoding="utf-8")
    print(f"Wrote time/speed package to {OUT}", flush=True)


if __name__ == "__main__":
    main()

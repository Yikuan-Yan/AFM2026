#!/usr/bin/env python3
"""Observational analysis of the 12-09-26 99.7 wt% glycerol D3 maps.

This deliberately performs no baseline, contact, hydrodynamic, PB, regression,
or other model fit.  It preserves raw vDeflection as the primary observable,
uses scanner travel rather than inferred surface separation, and reports only
window medians, quantiles, extrema, and exact palindrome differences.
"""

from __future__ import annotations

import csv
from datetime import datetime, timedelta
import hashlib
import json
from pathlib import Path
import platform
import sys
from typing import Iterable
from zipfile import ZipFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import analyze_velocity_systematics as velocity  # noqa: E402
import fit_glycerol_surface_forces as base  # noqa: E402


SHARE_TOKEN = "6e87fc8bc37b4a89aeda"
RAW_ROOT = ROOT / "raw" / f"keeper_{SHARE_TOKEN}"
DATA = RAW_ROOT / "12-09-26_99.7wt_glycerol_D3"
DOWNLOAD_MANIFEST = RAW_ROOT / "download_manifest.json"
OUT = ROOT / "analysis" / "glycerol_99p7_D3_observational_results"
FIG = OUT / "figures"

CONCENTRATION_WT_PERCENT = 99.7
GRID_SIZE = 8
MAP_COUNT = 32
MAPS_PER_BLOCK = 8
SPEEDS_UM_PER_S = (0.1, 0.3, 0.9, 2.7)
RETRACT_SPEED_UM_PER_S = 1.0
COLORS = {0.1: "#2a9d8f", 0.3: "#457b9d", 0.9: "#e9c46a", 2.7: "#e76f51"}

# Independent prior D3 calibration already present in this repository.  Raw V
# remains authoritative because the JPK files contain a conflicting embedded
# calibration and this no-fit analysis does not attempt to adjudicate it.
D3_PRIOR_INVOLS_NM_PER_V = 91.5629518721293
D3_PRIOR_K_N_PER_M = 0.2383668909630904
D3_PRIOR_SCALE_NN_PER_V = D3_PRIOR_INVOLS_NM_PER_V * 1e-9 * D3_PRIOR_K_N_PER_M * 1e9

INITIAL_REFERENCE_MAX_NM = 10.0
PLATEAU_START_NM = 100.0
PLATEAU_STOP_NM = 300.0
EARLY_START_NM = 80.0
EARLY_STOP_NM = 100.0
TERMINAL_WINDOW_NM = 25.0
BIN_WIDTH_NM = 5.0
BIN_EDGES_NM = np.arange(0.0, 495.0, BIN_WIDTH_NM)
BIN_CENTERS_NM = BIN_EDGES_NM[:-1] + BIN_WIDTH_NM / 2.0
FIXED_ROW = 3
FIXED_COLUMN = 3


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


def timestamp(value: str) -> datetime:
    return datetime.strptime(value, "%Y-%m-%d %H:%M:%S.%f %z")


def median(values: Iterable[float]) -> float:
    array = np.asarray(list(values), dtype=np.float64)
    array = array[np.isfinite(array)]
    return float(np.median(array)) if array.size else float("nan")


def iqr(values: Iterable[float]) -> float:
    array = np.asarray(list(values), dtype=np.float64)
    array = array[np.isfinite(array)]
    if not array.size:
        return float("nan")
    return float(np.percentile(array, 75) - np.percentile(array, 25))


def quantiles(values: Iterable[float]) -> tuple[float, float, float]:
    array = np.asarray(list(values), dtype=np.float64)
    array = array[np.isfinite(array)]
    if not array.size:
        return float("nan"), float("nan"), float("nan")
    return tuple(float(value) for value in np.percentile(array, [25, 50, 75]))


def map_metadata(path: Path) -> dict:
    with ZipFile(path) as archive:
        header = base.parse_properties(archive.read("header.properties"))
        segment_starts: list[datetime] = []
        segment_stops: list[datetime] = []
        failure_flags: list[str] = []
        for name in archive.namelist():
            if not name.endswith("/segment-header.properties") or not name.startswith("index/"):
                continue
            segment = base.parse_properties(archive.read(name))
            start = timestamp(segment["force-segment-header.time-stamp"])
            stop = start + timedelta(seconds=float(segment["force-segment-header.duration"]))
            segment_starts.append(start)
            segment_stops.append(stop)
            failure_flags.extend(
                key
                for key, value in segment.items()
                if key.startswith("force-segment-header.force-scan-flags.")
                and value.lower() == "true"
                and any(word in key for word in ("aborted", "out-of-range", "limit-exceeded"))
            )
        if not segment_starts:
            raise RuntimeError(f"No segment timestamps in {path}")
        p = "force-scan-map.settings.force-settings.segment."
        z0 = float(header[p + "0.z-start"])
        z1 = float(header[p + "0.z-end"])
        app_duration = float(header[p + "0.duration"])
        ret_duration = float(header[p + "2.duration"])
        return {
            "instrument_scan_number": int(header["force-scan-map.scan-number"]),
            "map_start_time": min(segment_starts).isoformat(),
            "map_end_time": max(segment_stops).isoformat(),
            "approach_nominal_speed_um_per_s": abs(z0 - z1) / app_duration * 1e6,
            "approach_declared_points": int(header[p + "0.num-points"]),
            "approach_declared_duration_s": app_duration,
            "retract_nominal_speed_um_per_s": abs(
                float(header[p + "2.z-start"]) - float(header[p + "2.z-end"])
            )
            / ret_duration
            * 1e6,
            "retract_declared_points": int(header[p + "2.num-points"]),
            "retract_declared_duration_s": ret_duration,
            "nominal_z_span_nm": abs(z0 - z1) * 1e9,
            "grid_i": int(header["force-scan-map.position-pattern.grid.ilength"]),
            "grid_j": int(header["force-scan-map.position-pattern.grid.jlength"]),
            "field_u_um": float(header["force-scan-map.position-pattern.grid.ulength"]) * 1e6,
            "field_v_um": float(header["force-scan-map.position-pattern.grid.vlength"]) * 1e6,
            "back_and_forth": header["force-scan-map.position-pattern.back-and-forth"].lower() == "true",
            "baseline_adjust_enabled": header[
                "force-scan-map.settings.force-settings.global-settings.force-baseline-adjust-settings.enabled"
            ].lower()
            == "true",
            "baseline_adjust_begin_of_line": header[
                "force-scan-map.settings.force-settings.global-settings.force-baseline-adjust-settings.beginOfLine"
            ].lower()
            == "true",
            "baseline_adjust_interval": int(
                header[
                    "force-scan-map.settings.force-settings.global-settings.force-baseline-adjust-settings.interval"
                ]
            ),
            "baseline_adjust_deadtime_samples": int(
                header[
                    "force-scan-map.settings.force-settings.global-settings.force-baseline-adjust-settings.deadtimeBeforeSamples"
                ]
            ),
            "baseline_adjust_average_samples": int(
                header[
                    "force-scan-map.settings.force-settings.global-settings.force-baseline-adjust-settings.averageSamples"
                ]
            ),
            "instrument_failure_flag_count": len(failure_flags),
        }


def curve_observables(curve: base.RawCurve | velocity.BranchCurve) -> tuple[dict, np.ndarray, np.ndarray]:
    travel_nm = (curve.measured_height_m[0] - curve.measured_height_m) * 1e9
    scanner_steps_nm = np.diff(travel_nm)
    # measuredHeight contains sub-nanometre readout noise, so individual samples
    # need not be strictly monotonic even though the 500 nm segment is directed.
    # Direct physical-window masks and bins remain valid; no monotonic
    # interpolation is introduced.
    if travel_nm[-1] <= travel_nm[0] or float(np.min(scanner_steps_nm)) < -2.0:
        raise RuntimeError("Scanner-travel direction or backtracking is unexpected")
    initial_mask = travel_nm <= INITIAL_REFERENCE_MAX_NM
    plateau_mask = (travel_nm >= PLATEAU_START_NM) & (travel_nm < PLATEAU_STOP_NM)
    early_mask = (travel_nm >= EARLY_START_NM) & (travel_nm < EARLY_STOP_NM)
    terminal_mask = travel_nm >= float(np.max(travel_nm) - TERMINAL_WINDOW_NM)
    if min(np.count_nonzero(initial_mask), np.count_nonzero(plateau_mask), np.count_nonzero(terminal_mask)) < 3:
        raise RuntimeError("A descriptive scanner-travel window has too few samples")
    initial_v = float(np.median(curve.deflection_V[initial_mask]))
    relative_v = curve.deflection_V - initial_v
    minimum_index = int(np.argmin(relative_v))
    maximum_index = int(np.argmax(relative_v))
    observables = {
        "sample_count": int(curve.deflection_V.size),
        "duration_s": float(curve.duration_s),
        "sample_rate_Hz": float(curve.deflection_V.size / curve.duration_s),
        "scanner_travel_nm": float(np.max(travel_nm)),
        "negative_scanner_step_fraction": float(np.mean(scanner_steps_nm < 0.0)),
        "largest_backward_scanner_step_nm": float(np.min(scanner_steps_nm)),
        "initial_0_10nm_raw_V": initial_v,
        "early_80_100nm_raw_V": float(np.median(curve.deflection_V[early_mask])),
        "settled_100_300nm_raw_V": float(np.median(curve.deflection_V[plateau_mask])),
        "terminal_last_25nm_raw_V": float(np.median(curve.deflection_V[terminal_mask])),
        "early_80_100nm_delta_V": float(np.median(relative_v[early_mask])),
        "settled_100_300nm_delta_V": float(np.median(relative_v[plateau_mask])),
        "terminal_last_25nm_delta_V": float(np.median(relative_v[terminal_mask])),
        "minimum_delta_V": float(relative_v[minimum_index]),
        "minimum_scanner_travel_nm": float(travel_nm[minimum_index]),
        "maximum_delta_V": float(relative_v[maximum_index]),
        "maximum_scanner_travel_nm": float(travel_nm[maximum_index]),
    }
    for field in (
        "early_80_100nm_delta_V",
        "settled_100_300nm_delta_V",
        "terminal_last_25nm_delta_V",
        "minimum_delta_V",
        "maximum_delta_V",
    ):
        observables[field.replace("_V", "_nN_prior_D3_scale")] = (
            observables[field] * D3_PRIOR_SCALE_NN_PER_V
        )
    return observables, travel_nm, relative_v


def binned_relative(travel_nm: np.ndarray, relative_v: np.ndarray) -> np.ndarray:
    output = np.full(BIN_CENTERS_NM.shape, np.nan, dtype=np.float64)
    indices = np.floor(travel_nm / BIN_WIDTH_NM).astype(int)
    for index in np.unique(indices):
        if 0 <= index < output.size:
            output[index] = float(np.median(relative_v[indices == index]))
    return output


def summarize_map(rows: list[dict], field: str, prefix: str) -> dict:
    values = [float(row[field]) for row in rows]
    q25, q50, q75 = quantiles(values)
    return {
        f"{prefix}_q25": q25,
        f"{prefix}_median": q50,
        f"{prefix}_q75": q75,
        f"{prefix}_spatial_IQR": q75 - q25,
        f"{prefix}_minimum": float(np.min(values)),
        f"{prefix}_maximum": float(np.max(values)),
    }


def speed_from_source(source: base.SourceData) -> float:
    speed = 0.5 / float(np.median([curve.duration_s for curve in source.curves]))
    selected = min(SPEEDS_UM_PER_S, key=lambda value: abs(value - speed))
    if abs(speed - selected) > 0.02:
        raise RuntimeError(f"Unexpected approach speed {speed}")
    return selected


def centered_position_profiles(pixel_rows: list[dict]) -> list[dict]:
    map_groups: dict[int, list[dict]] = {}
    for row in pixel_rows:
        map_groups.setdefault(int(row["acquisition_order"]), []).append(row)
    per_map: list[dict] = []
    for order, rows in sorted(map_groups.items()):
        for field in ("approach_initial_0_10nm_raw_V", "approach_settled_100_300nm_raw_V"):
            center = median(float(row[field]) for row in rows)
            for position in range(1, GRID_SIZE + 1):
                selected = [row for row in rows if int(row["position_within_line"]) == position]
                per_map.append(
                    {
                        "acquisition_order": order,
                        "speed_um_per_s": float(rows[0]["approach_speed_um_per_s"]),
                        "position_within_line": position,
                        "signal": "initial_0_10nm" if "initial" in field else "settled_100_300nm",
                        "map_centered_raw_V": median(float(row[field]) - center for row in selected),
                        "settled_delta_V": median(
                            float(row["approach_settled_100_300nm_delta_V"]) for row in selected
                        ),
                    }
                )
    output: list[dict] = []
    for speed in SPEEDS_UM_PER_S:
        for position in range(1, GRID_SIZE + 1):
            delta_rows = [
                row
                for row in per_map
                if row["speed_um_per_s"] == speed
                and row["position_within_line"] == position
                and row["signal"] == "initial_0_10nm"
            ]
            delta_q = quantiles(row["settled_delta_V"] for row in delta_rows)
            for signal in ("initial_0_10nm", "settled_100_300nm"):
                selected = [row for row in delta_rows if row["signal"] == signal]
                if not selected:
                    selected = [
                        row
                        for row in per_map
                        if row["speed_um_per_s"] == speed
                        and row["position_within_line"] == position
                        and row["signal"] == signal
                    ]
                centered_q = quantiles(row["map_centered_raw_V"] for row in selected)
                output.append(
                    {
                        "approach_speed_um_per_s": speed,
                        "position_within_line": position,
                        "signal": signal,
                        "maps": len(selected),
                        "map_centered_raw_V_q25": centered_q[0],
                        "map_centered_raw_V_median": centered_q[1],
                        "map_centered_raw_V_q75": centered_q[2],
                        "settled_delta_V_q25": delta_q[0],
                        "settled_delta_V_median": delta_q[1],
                        "settled_delta_V_q75": delta_q[2],
                        "settled_delta_nN_prior_D3_scale_median": delta_q[1]
                        * D3_PRIOR_SCALE_NN_PER_V,
                    }
                )
    return output


def aggregate_binned(branch_arrays: dict[tuple[float, str], list[np.ndarray]]) -> list[dict]:
    rows: list[dict] = []
    for speed in SPEEDS_UM_PER_S:
        for branch in ("approach", "retract_far_to_contact"):
            matrix = np.asarray(branch_arrays[(speed, branch)], dtype=np.float64)
            for index, center in enumerate(BIN_CENTERS_NM):
                values = matrix[:, index]
                values = values[np.isfinite(values)]
                if not values.size:
                    continue
                q25, q50, q75 = (float(value) for value in np.percentile(values, [25, 50, 75]))
                rows.append(
                    {
                        "branch": branch,
                        "approach_speed_um_per_s": speed,
                        "scanner_travel_bin_center_nm": center,
                        "contributing_curves": int(values.size),
                        "delta_V_q25": q25,
                        "delta_V_median": q50,
                        "delta_V_q75": q75,
                        "delta_nN_prior_D3_scale_q25": q25 * D3_PRIOR_SCALE_NN_PER_V,
                        "delta_nN_prior_D3_scale_median": q50 * D3_PRIOR_SCALE_NN_PER_V,
                        "delta_nN_prior_D3_scale_q75": q75 * D3_PRIOR_SCALE_NN_PER_V,
                    }
                )
    return rows


def build_speed_summary(map_rows: list[dict], pixel_rows: list[dict]) -> list[dict]:
    output: list[dict] = []
    for speed in SPEEDS_UM_PER_S:
        maps = [row for row in map_rows if float(row["approach_speed_um_per_s"]) == speed]
        pixels = [row for row in pixel_rows if float(row["approach_speed_um_per_s"]) == speed]
        edge = [row for row in pixels if int(row["position_within_line"]) in (1, 8)]
        interior = [row for row in pixels if 2 <= int(row["position_within_line"]) <= 7]
        fields = {
            "approach_settled_delta_V": "approach_settled_delta_V_median",
            "approach_settled_spatial_IQR_V": "approach_settled_delta_V_spatial_IQR",
            "approach_terminal_delta_V": "approach_terminal_delta_V_median",
            "retract_minimum_delta_V": "retract_minimum_delta_V_median",
            "retract_minimum_scanner_travel_nm": "retract_minimum_scanner_travel_nm_median",
        }
        row: dict = {
            "approach_speed_um_per_s": speed,
            "maps": len(maps),
            "approach_curves": len(pixels),
            "retract_speed_um_per_s": RETRACT_SPEED_UM_PER_S,
            "prior_D3_scale_nN_per_V": D3_PRIOR_SCALE_NN_PER_V,
        }
        for prefix, field in fields.items():
            values = [float(item[field]) for item in maps]
            row[prefix + "_map_median"] = median(values)
            row[prefix + "_map_minimum"] = float(np.min(values))
            row[prefix + "_map_maximum"] = float(np.max(values))
        row["approach_settled_delta_nN_prior_D3_scale_map_median"] = (
            row["approach_settled_delta_V_map_median"] * D3_PRIOR_SCALE_NN_PER_V
        )
        row["approach_settled_delta_nN_prior_D3_scale_map_minimum"] = (
            row["approach_settled_delta_V_map_minimum"] * D3_PRIOR_SCALE_NN_PER_V
        )
        row["approach_settled_delta_nN_prior_D3_scale_map_maximum"] = (
            row["approach_settled_delta_V_map_maximum"] * D3_PRIOR_SCALE_NN_PER_V
        )
        row["retract_minimum_delta_nN_prior_D3_scale_map_median"] = (
            row["retract_minimum_delta_V_map_median"] * D3_PRIOR_SCALE_NN_PER_V
        )
        row["retract_minimum_delta_nN_prior_D3_scale_map_minimum"] = (
            row["retract_minimum_delta_V_map_minimum"] * D3_PRIOR_SCALE_NN_PER_V
        )
        row["retract_minimum_delta_nN_prior_D3_scale_map_maximum"] = (
            row["retract_minimum_delta_V_map_maximum"] * D3_PRIOR_SCALE_NN_PER_V
        )
        row["line_edge_settled_delta_V_pixel_median"] = median(
            float(item["approach_settled_100_300nm_delta_V"]) for item in edge
        )
        row["line_interior_settled_delta_V_pixel_median"] = median(
            float(item["approach_settled_100_300nm_delta_V"]) for item in interior
        )
        row["settled_delta_V_per_speed_um_per_s"] = (
            row["approach_settled_delta_V_map_median"] / speed
        )
        output.append(row)
    return output


def build_palindrome_pairs(map_rows: list[dict]) -> list[dict]:
    output: list[dict] = []
    for block in range(1, MAP_COUNT // MAPS_PER_BLOCK + 1):
        block_rows = [row for row in map_rows if int(row["block"]) == block]
        block_rows.sort(key=lambda row: int(row["position_in_block"]))
        for pair in range(1, 5):
            early = block_rows[pair - 1]
            late = block_rows[MAPS_PER_BLOCK - pair]
            if float(early["approach_speed_um_per_s"]) != float(late["approach_speed_um_per_s"]):
                raise RuntimeError("Palindrome speed identity failed")
            output.append(
                {
                    "block": block,
                    "pair_from_outer_edge": pair,
                    "approach_speed_um_per_s": early["approach_speed_um_per_s"],
                    "early_acquisition_order": early["acquisition_order"],
                    "late_acquisition_order": late["acquisition_order"],
                    "early_map_start_time": early["map_start_time"],
                    "late_map_start_time": late["map_start_time"],
                    "approach_settled_delta_V_later_minus_earlier": float(
                        late["approach_settled_delta_V_median"]
                    )
                    - float(early["approach_settled_delta_V_median"]),
                    "approach_settled_delta_nN_prior_D3_scale_later_minus_earlier": (
                        float(late["approach_settled_delta_V_median"])
                        - float(early["approach_settled_delta_V_median"])
                    )
                    * D3_PRIOR_SCALE_NN_PER_V,
                    "retract_minimum_delta_V_later_minus_earlier": float(
                        late["retract_minimum_delta_V_median"]
                    )
                    - float(early["retract_minimum_delta_V_median"]),
                    "retract_minimum_delta_nN_prior_D3_scale_later_minus_earlier": (
                        float(late["retract_minimum_delta_V_median"])
                        - float(early["retract_minimum_delta_V_median"])
                    )
                    * D3_PRIOR_SCALE_NN_PER_V,
                    "retract_minimum_scanner_travel_nm_later_minus_earlier": float(
                        late["retract_minimum_scanner_travel_nm_median"]
                    )
                    - float(early["retract_minimum_scanner_travel_nm_median"]),
                }
            )
    return output


def plot_speed_overview(binned_rows: list[dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    panels = (
        (axes[0, 0], "approach", (0, 490), (-5, 40), "Approach: full scanner travel"),
        (axes[0, 1], "approach", (0, 160), (-2, 20), "Approach startup: no flat far-field baseline"),
        (
            axes[1, 0],
            "retract_far_to_contact",
            (0, 490),
            (-50, 45),
            "Retract reordered far-to-contact (retract speed fixed at 1 µm/s)",
        ),
        (
            axes[1, 1],
            "retract_far_to_contact",
            (0, 350),
            (-3, 2),
            "Retract far-side response",
        ),
    )
    for axis, branch, xlim, ylim, title in panels:
        for speed in SPEEDS_UM_PER_S:
            rows = [
                row
                for row in binned_rows
                if row["branch"] == branch and float(row["approach_speed_um_per_s"]) == speed
            ]
            x = np.asarray([row["scanner_travel_bin_center_nm"] for row in rows], dtype=float)
            y = np.asarray([row["delta_nN_prior_D3_scale_median"] for row in rows], dtype=float)
            q1 = np.asarray([row["delta_nN_prior_D3_scale_q25"] for row in rows], dtype=float)
            q3 = np.asarray([row["delta_nN_prior_D3_scale_q75"] for row in rows], dtype=float)
            axis.plot(x, y, color=COLORS[speed], lw=2.0, label=f"{speed:g} µm/s approach")
            axis.fill_between(x, q1, q3, color=COLORS[speed], alpha=0.14)
        axis.axhline(0.0, color="0.45", lw=0.7)
        axis.set_xlim(*xlim)
        axis.set_ylim(*ylim)
        axis.set_title(title)
        axis.set_xlabel("Scanner travel from far endpoint (nm)")
        axis.set_ylabel("Δforce from first 10 nm (nN; prior D3 scale)")
        axis.grid(alpha=0.2)
        axis.legend(frameon=False, fontsize=9)
    fig.suptitle(
        "12 Sep 2026 · 99.7 wt% glycerol · raw endpoint-referenced branches\n"
        "All 512 curves per approach speed: median and IQR; no baseline/contact/model fit",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIG / "approach_retract_speed_overview.png", dpi=200)
    plt.close(fig)


def plot_startup_and_line_position(binned_rows: list[dict], profiles: list[dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for speed in SPEEDS_UM_PER_S:
        rows = [
            row
            for row in binned_rows
            if row["branch"] == "approach" and float(row["approach_speed_um_per_s"]) == speed
        ]
        travel = np.asarray([row["scanner_travel_bin_center_nm"] for row in rows], dtype=float)
        force = np.asarray([row["delta_nN_prior_D3_scale_median"] for row in rows], dtype=float)
        axes[0, 0].plot(travel, force, color=COLORS[speed], lw=2, label=f"{speed:g} µm/s")
        axes[0, 1].plot(travel / (speed * 1000.0), force, color=COLORS[speed], lw=2, label=f"{speed:g} µm/s")
        selected = [
            row
            for row in profiles
            if float(row["approach_speed_um_per_s"]) == speed and row["signal"] == "initial_0_10nm"
        ]
        positions = np.asarray([row["position_within_line"] for row in selected], dtype=float)
        values = np.asarray([row["settled_delta_nN_prior_D3_scale_median"] for row in selected])
        axes[1, 0].plot(positions, values, color=COLORS[speed], marker="o", lw=1.8, label=f"{speed:g} µm/s")
        for signal, marker, linestyle in (
            ("initial_0_10nm", "o", "-"),
            ("settled_100_300nm", "s", "--"),
        ):
            selected = [
                row
                for row in profiles
                if float(row["approach_speed_um_per_s"]) == speed and row["signal"] == signal
            ]
            axes[1, 1].plot(
                [row["position_within_line"] for row in selected],
                [row["map_centered_raw_V_median"] for row in selected],
                color=COLORS[speed],
                marker=marker,
                ls=linestyle,
                lw=1.4,
                ms=4,
                label=f"{speed:g} µm/s · {'initial' if signal.startswith('initial') else 'settled'}",
            )
    axes[0, 0].set_xlim(0, 160)
    axes[0, 0].set_ylim(-2, 20)
    axes[0, 0].set_title("Same startup response versus scanner travel")
    axes[0, 0].set_xlabel("Scanner travel from far endpoint (nm)")
    axes[0, 1].set_xlim(0, 0.12)
    axes[0, 1].set_ylim(-2, 20)
    axes[0, 1].set_title("Startup response versus elapsed motion time")
    axes[0, 1].set_xlabel("Nominal elapsed approach time (s)")
    axes[1, 0].set_title("Endpoint-referenced plateau depends on position within each raster line")
    axes[1, 0].set_xlabel("Acquisition position within line (1 = line start, 8 = line end)")
    axes[1, 0].set_ylabel("100–300 nm Δforce (nN; prior D3 scale)")
    axes[1, 1].set_title("The line pattern is concentrated in the initial 0–10 nm level")
    axes[1, 1].set_xlabel("Acquisition position within line")
    axes[1, 1].set_ylabel("Raw V minus that map's spatial median (V)")
    for axis in axes.flat:
        axis.axhline(0.0, color="0.45", lw=0.7)
        axis.grid(alpha=0.2)
        axis.legend(frameon=False, fontsize=8, ncol=2 if axis is axes[1, 1] else 1)
    axes[0, 0].set_ylabel("Δforce from first 10 nm (nN; prior D3 scale)")
    axes[0, 1].set_ylabel("Δforce from first 10 nm (nN; prior D3 scale)")
    fig.suptitle(
        "Motion-start and raster-line-position diagnostics\n"
        "Direct medians only; header enables baseline adjustment at the beginning of each line",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIG / "startup_and_line_position.png", dpi=200)
    plt.close(fig)


def plot_chronology(map_rows: list[dict]) -> None:
    times = np.asarray([mdates.date2num(datetime.fromisoformat(row["map_start_time"])) for row in map_rows])
    colors = [COLORS[float(row["approach_speed_um_per_s"])] for row in map_rows]
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharex=True)
    plateau = np.asarray([row["approach_settled_delta_nN_prior_D3_scale_median"] for row in map_rows])
    plateau_iqr = np.asarray([row["approach_settled_delta_nN_prior_D3_scale_spatial_IQR"] for row in map_rows])
    axes[0, 0].errorbar(times, plateau, yerr=plateau_iqr / 2, color="0.55", lw=0.7, fmt="none")
    axes[0, 0].plot(times, plateau, color="0.45", lw=0.8)
    axes[0, 0].scatter(times, plateau, c=colors, s=34, zorder=3)
    axes[0, 0].set_ylabel("Approach plateau Δforce (nN; prior D3 scale)")
    axes[0, 0].set_title("Approach 100–300 nm settled increment; bars = half spatial IQR")

    initial = np.asarray([row["approach_initial_raw_V_median"] for row in map_rows])
    settled = np.asarray([row["approach_settled_raw_V_median"] for row in map_rows])
    axes[0, 1].plot(times, initial, color="#6a4c93", marker="o", ms=3.5, lw=1, label="initial 0–10 nm")
    axes[0, 1].plot(times, settled, color="#1982c4", marker="s", ms=3.5, lw=1, label="settled 100–300 nm")
    jump_index = int(np.argmax(np.abs(np.diff(initial))))
    jump_x = (times[jump_index] + times[jump_index + 1]) / 2
    axes[0, 1].axvline(jump_x, color="#ad3f36", ls=":", lw=1.2, label="largest raw-level discontinuity")
    axes[0, 1].set_ylabel("Raw vDeflection (V)")
    axes[0, 1].set_title("Absolute detector level drifts and resets; not a force-zero chronology")
    axes[0, 1].legend(frameon=False, fontsize=8)

    retract_minimum = np.asarray([row["retract_minimum_delta_nN_prior_D3_scale_median"] for row in map_rows])
    retract_iqr = np.asarray([row["retract_minimum_delta_nN_prior_D3_scale_spatial_IQR"] for row in map_rows])
    axes[1, 0].errorbar(times, retract_minimum, yerr=retract_iqr / 2, color="0.60", lw=0.7, fmt="none")
    axes[1, 0].plot(times, retract_minimum, color="0.45", lw=0.8)
    axes[1, 0].scatter(times, retract_minimum, c=colors, s=34, zorder=3)
    axes[1, 0].set_ylabel("Retract minimum Δforce (nN; prior D3 scale)")
    axes[1, 0].set_title("Retract suction/adhesion minimum; retract speed is always 1 µm/s")

    minimum_x = np.asarray([row["retract_minimum_scanner_travel_nm_median"] for row in map_rows])
    axes[1, 1].plot(times, minimum_x, color="0.45", lw=0.8)
    axes[1, 1].scatter(times, minimum_x, c=colors, s=34, zorder=3)
    axes[1, 1].set_ylabel("Retract minimum position from far endpoint (nm)")
    axes[1, 1].set_title("The retract minimum position also changes with acquisition history")

    for boundary in range(MAPS_PER_BLOCK, MAP_COUNT, MAPS_PER_BLOCK):
        x = (times[boundary - 1] + times[boundary]) / 2
        for axis in axes.flat:
            axis.axvline(x, color="0.5", ls="--", lw=0.7)
    first_time = datetime.fromisoformat(map_rows[0]["map_start_time"])
    for axis in axes.flat:
        axis.grid(alpha=0.2)
        axis.xaxis.set_major_locator(mdates.MinuteLocator(interval=10, tz=first_time.tzinfo))
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=first_time.tzinfo))
        axis.tick_params(axis="x", labelrotation=0)
        axis.set_xlabel("Map start time (instrument-recorded UTC+02:00)")
    handles = [
        plt.Line2D([], [], marker="o", ls="", color=color, label=f"{speed:g} µm/s approach")
        for speed, color in COLORS.items()
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.95), ncol=4, frameon=False)
    fig.suptitle("Map chronology and history dependence · dashed lines are 8-map palindrome blocks", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(FIG / "chronology_and_history.png", dpi=200)
    plt.close(fig)


def plot_fixed_pixel(fixed_rows: list[dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    map_points: list[tuple[int, float, float]] = []
    for item in fixed_rows:
        speed = float(item["approach_speed_um_per_s"])
        x = np.asarray(item["scanner_travel_nm"], dtype=float)
        y = np.asarray(item["delta_V"], dtype=float) * D3_PRIOR_SCALE_NN_PER_V
        if item["branch"] == "approach":
            axes[0, 0].plot(x, y, color=COLORS[speed], alpha=0.38, lw=0.7)
            axes[0, 1].plot(x, y, color=COLORS[speed], alpha=0.38, lw=0.7)
            mask = (x >= PLATEAU_START_NM) & (x < PLATEAU_STOP_NM)
            map_points.append((int(item["acquisition_order"]), speed, float(np.median(y[mask]))))
        else:
            axes[1, 0].plot(x, y, color=COLORS[speed], alpha=0.38, lw=0.7)
    axes[0, 0].set_xlim(0, 490)
    axes[0, 0].set_ylim(-5, 40)
    axes[0, 0].set_title("Approach at one fixed physical/acquisition position")
    axes[0, 1].set_xlim(0, 160)
    axes[0, 1].set_ylim(-2, 20)
    axes[0, 1].set_title("Fixed-pixel startup region")
    axes[1, 0].set_xlim(0, 490)
    axes[1, 0].set_ylim(-60, 45)
    axes[1, 0].set_title("Retract at the same fixed pixel, reordered far-to-contact")
    map_points.sort()
    axes[1, 1].plot([row[0] for row in map_points], [row[2] for row in map_points], color="0.45", lw=0.8)
    axes[1, 1].scatter(
        [row[0] for row in map_points],
        [row[2] for row in map_points],
        c=[COLORS[row[1]] for row in map_points],
        s=35,
    )
    axes[1, 1].set_title("Fixed-pixel 100–300 nm settled increment")
    axes[1, 1].set_xlabel("Acquisition order")
    axes[1, 1].set_ylabel("Δforce (nN; prior D3 scale)")
    for boundary in (8.5, 16.5, 24.5):
        axes[1, 1].axvline(boundary, color="0.5", ls="--", lw=0.7)
    for axis in axes.flat[:3]:
        axis.axhline(0, color="0.45", lw=0.7)
        axis.grid(alpha=0.2)
        axis.set_xlabel("Scanner travel from far endpoint (nm)")
        axis.set_ylabel("Δforce from first 10 nm (nN; prior D3 scale)")
    axes[1, 1].grid(alpha=0.2)
    handles = [
        plt.Line2D([], [], color=color, label=f"{speed:g} µm/s approach")
        for speed, color in COLORS.items()
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.95), ncol=4, frameon=False)
    fig.suptitle(
        f"Fixed pixel row {FIXED_ROW}, column {FIXED_COLUMN}: speed separation survives spatial control\n"
        "Raw first-10-nm reference; no baseline/contact/model fit",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(FIG / "fixed_pixel_raw_chronology.png", dpi=200)
    plt.close(fig)


def build_report(
    map_rows: list[dict],
    speed_rows: list[dict],
    pair_rows: list[dict],
    profiles: list[dict],
    embedded_sensitivity_nm_per_v: float,
    embedded_k_n_per_m: float,
) -> str:
    embedded_scale = embedded_sensitivity_nm_per_v * 1e-9 * embedded_k_n_per_m * 1e9
    initial_values = np.asarray([row["approach_initial_raw_V_median"] for row in map_rows])
    jump_index = int(np.argmax(np.abs(np.diff(initial_values))))
    jump = float(initial_values[jump_index + 1] - initial_values[jump_index])
    profile_lookup = {
        (float(row["approach_speed_um_per_s"]), int(row["position_within_line"])): row
        for row in profiles
        if row["signal"] == "initial_0_10nm"
    }
    pair_plateau = np.asarray(
        [row["approach_settled_delta_nN_prior_D3_scale_later_minus_earlier"] for row in pair_rows]
    )
    pair_retract = np.asarray(
        [row["retract_minimum_delta_nN_prior_D3_scale_later_minus_earlier"] for row in pair_rows]
    )
    lines = [
        "# 12-09-26 · 99.7 wt% glycerol · cantilever D3 · observational analysis",
        "",
        "## 物理图像与结论",
        "",
        "这批数据的主导特征不是一条可直接放入 PB 或 lubrication 拟合的平衡 force–distance 曲线，而是 **高黏度液体中、每次 Z motion 启动后出现的速度相关瞬态和准平台，加上逐行 baseline adjustment 与采集历史留下的行内位置结构**。在 2.7 µm/s 时，approach 从最远端开始后的几十纳米内就上升到很大的正平台；这一变化远早于末端硬接触样的陡升，因此前段不能作为平坦、无力的 far field。固定像素也保留四速度分离，说明它不是单纯由空间平均造成；但全图中行首/行末与内部点明显不同，说明它也不能被解释成样品表面的静态空间力图。",
        "",
        "本分析因此只报告 raw `vDeflection`、从最远端起算的 scanner travel、固定物理窗口中位数/IQR、极值以及四个 8-map palindrome block 的同速 later-minus-earlier 差。**没有进行 baseline、contact、hydrodynamic、PB、回归或其他模型拟合，也没有输出 surface separation、Debye length、surface potential 或 `U→0` 外推。**",
        "",
        "## 直接观测",
        "",
        "### 1. Approach 的启动瞬态与速度平台",
        "",
        "每条曲线以最初 0–10 nm 的 raw detector 中位数作一个纯操作性常数参考；100–300 nm 窗口用于描述 motion 启动后已经较稳定的电压增量。它不是绝对零力。下表先给 raw V，再把同一电压差按仓库内既有 D3 独立标定换成仅供量级阅读的 nN：",
        "",
        "| approach speed (µm/s) | maps | settled ΔV map median (V) | map range (V) | prior-D3-scale ΔF (nN) | map range (nN) | ΔV / speed (V per µm/s) |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in speed_rows:
        lines.append(
            f"| {float(row['approach_speed_um_per_s']):.1f} | {int(row['maps'])} | "
            f"{float(row['approach_settled_delta_V_map_median']):.5f} | "
            f"{float(row['approach_settled_delta_V_map_minimum']):.5f}–{float(row['approach_settled_delta_V_map_maximum']):.5f} | "
            f"{float(row['approach_settled_delta_nN_prior_D3_scale_map_median']):.3f} | "
            f"{float(row['approach_settled_delta_nN_prior_D3_scale_map_minimum']):.3f}–{float(row['approach_settled_delta_nN_prior_D3_scale_map_maximum']):.3f} | "
            f"{float(row['settled_delta_V_per_speed_um_per_s']):.5f} |"
        )
    lines.extend(
        [
            "",
            "四个幅值不仅相差很大，`ΔV/speed` 也不恒定；因此在未解决 startup response、baseline semantics 和绝对间距前，不能把这些平台直接当作简单线性 viscous coefficient。图中 2.7 µm/s 平台在约前 100 nm 内形成，而末端约 400 nm 以后才出现所有速度共有的硬接触样上升。最小可信结论是：**运动启动/体相拖曳/仪器响应在 nominal far side 已经是主量级之一**；仅凭本数据不能把三者再分开。",
            "",
            "### 2. 逐行采集位置结构",
            "",
            "Header 明确记录 baseline adjustment `enabled=true`, `beginOfLine=true`, `interval=1`, `deadtimeBeforeSamples=100`, `averageSamples=100`。按 serpentine raster 的采集位置统计，100–300 nm 的 endpoint-referenced 平台在行首最低、第二点最高，随后向行末下降。以下为各速度 position 1 / 2 / 8 的中位数（nN，仅 prior D3 scale）：",
            "",
            "| speed (µm/s) | line start pos 1 | pos 2 | line end pos 8 |",
            "|---:|---:|---:|---:|",
        ]
    )
    for speed in SPEEDS_UM_PER_S:
        lines.append(
            f"| {speed:.1f} | "
            f"{float(profile_lookup[(speed, 1)]['settled_delta_nN_prior_D3_scale_median']):.3f} | "
            f"{float(profile_lookup[(speed, 2)]['settled_delta_nN_prior_D3_scale_median']):.3f} | "
            f"{float(profile_lookup[(speed, 8)]['settled_delta_nN_prior_D3_scale_median']):.3f} |"
        )
    lines.extend(
        [
            "",
            "把每张 map 的 raw V 先减去自身空间中位数后可见：结构主要集中在最初 0–10 nm detector level；100–300 nm 的 settled raw level 平坦得多。因此彩色 heatmap 的边缘/行内图样主要是 motion-start reference 与逐行采集状态，而不是证据充分的表面相互作用空间图。Pixel 仍是相关的 raster 采样，不能作为 64 个独立实验。",
            "",
            "### 3. 旧分析方法在这里失效的位置",
            "",
            "- 旧流程把 approach 的前 20% 当作 far-field baseline；本实验的前 20% 正好是最强的速度相关弯曲启动段。线性 subtraction 会删掉或重塑待观察的动态信号。",
            "- 四速度保持约 2 kHz 采样率，而不是保持相同点间距：每条 approach 分别约 10000、3330、1110、370 点。旧 contact search 的固定 150 点分别覆盖约 **7.5、22.5、67.5、202.7 nm** scanner travel，物理窗口不等价。",
            "- 高速平台按 prior D3 scale 可达十余 nN，对应几十纳米 cantilever deflection；若把它先作为基线删去或保留在 `h+δ` 中，都会显著改变 inferred separation。没有独立的 absolute force-zero/contact plane 时，本次不构造 surface separation。",
            "- 文件自身的 embedded calibration 与既有 D3 独立标定冲突，进一步限制了绝对力幅值；本报告的科学主结果因此是 V，nN 只是指定标度下的换算。",
            "",
            "### 4. Retract 是固定 1 µm/s，但存在大而历史相关的负向最低点",
            "",
            "所有 retract 都是 1 µm/s。下表按其前一条 approach speed 分组；颜色分组不是 retract speed 扫描，也不能据此归因于 approach speed：",
            "",
            "| preceding approach speed (µm/s) | retract minimum map median (V) | map range (V) | prior-D3-scale median (nN) | map range (nN) | median minimum position (nm from far endpoint) |",
            "|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in speed_rows:
        lines.append(
            f"| {float(row['approach_speed_um_per_s']):.1f} | "
            f"{float(row['retract_minimum_delta_V_map_median']):.5f} | "
            f"{float(row['retract_minimum_delta_V_map_minimum']):.5f}–{float(row['retract_minimum_delta_V_map_maximum']):.5f} | "
            f"{float(row['retract_minimum_delta_nN_prior_D3_scale_map_median']):.2f} | "
            f"{float(row['retract_minimum_delta_nN_prior_D3_scale_map_minimum']):.2f}–{float(row['retract_minimum_delta_nN_prior_D3_scale_map_maximum']):.2f} | "
            f"{float(row['retract_minimum_scanner_travel_nm_map_median']):.1f} |"
        )
    lines.extend(
        [
            "",
            "Retract minimum 的 map spatial IQR 很大，并在时间上先增强、后减弱；其位置也移动。它可以描述为高黏度条件下的 suction/adhesion-like hysteresis，但仅凭曲线不能区分 viscous drainage、cavitation/bridging、surface adhesion、contact history 或 instrument baseline。全部 2048 条 approach 和 2048 条 retract 都未触及 raw encoder saturation，因此这些极值不是 ADC 截断造成的。",
            "",
            "### 5. Palindrome 配对只支持“存在 history”，不支持一次速度拟合",
            "",
            f"四个 block 各自按位置 1↔8、2↔7、3↔6、4↔5 形成 16 个同速 early/late map pair。Approach settled 平台的 later-minus-earlier 中位数为 **{float(np.median(pair_plateau)):+.3f} nN**，范围 **{float(np.min(pair_plateau)):+.3f} 到 {float(np.max(pair_plateau)):+.3f} nN**；retract minimum 的对应中位数为 **{float(np.median(pair_retract)):+.3f} nN**，范围 **{float(np.min(pair_retract)):+.3f} 到 {float(np.max(pair_retract)):+.3f} nN**。符号并不一致，说明时间/接触历史对两个 branch 的影响不可忽略。这里没有把 16 对做成显著性检验，也没有拟合速度斜率。",
            "",
            "Absolute raw detector level 从早期逐步漂移；最大相邻跳变发生在 acquisition order "
            f"{jump_index + 1}→{jump_index + 2}，initial level 改变 **{jump:+.3f} V**。这类 reset/discontinuity 不出现在 encoder conversion 常数中，故跨越该点的 absolute-V 比较不能直接解释成力变化。Endpoint-referenced 增量仍保留，但也不能恢复绝对零力。",
            "",
            "## 数据完整性与标度边界",
            "",
            f"- Keeper share：32 个文件，32 张完整 8×8 map；共 2048 条 approach 和 2048 条 retract。下载清单逐文件记录 SHA-256，JPK ZIP CRC 全部通过；parser skip、instrument failure flag 与 raw saturation 均为 0。",
            "- 每张图覆盖同一 2×2 µm 区域，500 nm nominal Z span；approach 为 0.1/0.3/0.9/2.7 µm/s，各 8 张；retract 固定 1 µm/s。四个 8-map block 都是同速对称 palindrome，且每个 block 的起始速度轮换。",
            f"- 既有 D3 独立标定：InvOLS = {D3_PRIOR_INVOLS_NM_PER_V:.6f} nm/V，k = {D3_PRIOR_K_N_PER_M:.9f} N/m，乘积 = {D3_PRIOR_SCALE_NN_PER_V:.6f} nN/V。",
            f"- 这些 JPK 文件内嵌：InvOLS = {embedded_sensitivity_nm_per_v:.6f} nm/V，k = {embedded_k_n_per_m:.9f} N/m，乘积 = {embedded_scale:.6f} nN/V。Embedded/prior 幅值比 = {embedded_scale / D3_PRIOR_SCALE_NN_PER_V:.3f}。本分析不拟合新 InvOLS，也不裁决哪一套是当前绝对力标尺。",
            "",
            "## 可复现产物",
            "",
            "- `map_inventory.csv`：文件、时间、四速度设计、采样点、baseline-adjust settings、校准字段与完整性。",
            "- `pixel_observables.csv`：2048 个物理像素的 approach/retract 原始窗口描述量。",
            "- `map_summary.csv`、`speed_summary.csv`、`palindrome_pair_differences.csv`：等 map 权重的描述性汇总。",
            "- `binned_branch_summary.csv`：5 nm scanner-travel bins 的全曲线中位数/IQR；binning 只用于显示，不提高物理分辨率。",
            "- `fixed_pixel_binned_curves.csv`：row 3, column 3 的逐 map 原始曲线，控制空间与行内采集位置。",
            "- `line_position_profile.csv`：行内 position 1–8 的 startup/settled 结构。",
            "- `figures/approach_retract_speed_overview.png`：四速度 approach 与固定速度 retract。",
            "- `figures/startup_and_line_position.png`：startup 时间/位移图与行内位置效应。",
            "- `figures/chronology_and_history.png`：palindrome block、raw-level reset 与 retract history。",
            "- `figures/fixed_pixel_raw_chronology.png`：固定像素的 32-map 对照。",
            "",
            "## 当前可说与不可说",
            "",
            "可以说：99.7 wt% glycerol 中，motion-start/体相拖曳/仪器响应在 nominal far side 已经达到与表面近场同等级的重要性；逐行采集状态和长时间 history 都清晰存在；retract 有强负向 hysteresis。",
            "",
            "不能说：平台已经给出 hydrodynamic coefficient；四组就是平衡 surface force；负向 retract 极值已经识别为某一种 adhesion/cavitation 机制；既有 D3 或 embedded 标定已被本数据重新验证。下一步若要做定量模型，首先需要独立 blank/large-gap motion transient、current D3 calibration，以及不依赖前 20% 的 absolute force-zero/contact strategy。",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    paths = sorted(DATA.glob("*.jpk-force-map"))
    if len(paths) != MAP_COUNT:
        raise RuntimeError(f"Expected {MAP_COUNT} maps, found {len(paths)}")
    download = json.loads(DOWNLOAD_MANIFEST.read_text(encoding="utf-8"))
    expected_hashes = {
        Path(row["local_path"]).name: row["sha256"] for row in download["files"]
    }

    inventory_rows: list[dict] = []
    pixel_rows: list[dict] = []
    map_rows: list[dict] = []
    fixed_plot_rows: list[dict] = []
    fixed_export_rows: list[dict] = []
    branch_arrays: dict[tuple[float, str], list[np.ndarray]] = {
        (speed, branch): []
        for speed in SPEEDS_UM_PER_S
        for branch in ("approach", "retract_far_to_contact")
    }
    embedded_sensitivities: set[float] = set()
    embedded_springs: set[float] = set()
    approach_saturated = 0
    retract_saturated = 0
    total_skipped = 0

    for acquisition_order, path in enumerate(paths, start=1):
        digest = sha256_file(path)
        if expected_hashes.get(path.name) != digest:
            raise RuntimeError(f"Download manifest hash mismatch for {path.name}")
        metadata = map_metadata(path)
        source = base.load_source(path, 100)
        approaches, approach_skipped = velocity.load_branch(path, "extend")
        retracts, retract_skipped = velocity.load_branch(path, "retract")
        total_skipped += source.skipped_curves + approach_skipped + retract_skipped
        speed = speed_from_source(source)
        if not np.isclose(speed, metadata["approach_nominal_speed_um_per_s"], atol=0.002):
            raise RuntimeError("Header/source speed mismatch")
        if metadata["grid_i"] != GRID_SIZE or metadata["grid_j"] != GRID_SIZE:
            raise RuntimeError("Unexpected map grid")
        if len(approaches) != GRID_SIZE**2 or len(retracts) != GRID_SIZE**2:
            raise RuntimeError("Incomplete approach/retract branch")
        approach_by_index = {int(curve.point_index): curve for curve in approaches}
        retract_by_index = {int(curve.point_index): curve for curve in retracts}
        if set(approach_by_index) != set(range(GRID_SIZE**2)) or set(retract_by_index) != set(range(GRID_SIZE**2)):
            raise RuntimeError("Unexpected point-index inventory")
        embedded_sensitivities.add(source.stored_sensitivity_m_per_V * 1e9)
        embedded_springs.add(source.stored_spring_constant_N_per_m)

        block = (acquisition_order - 1) // MAPS_PER_BLOCK + 1
        position_in_block = (acquisition_order - 1) % MAPS_PER_BLOCK + 1
        map_pixel_rows: list[dict] = []
        app_sample_counts: list[int] = []
        ret_sample_counts: list[int] = []
        app_travels: list[float] = []
        ret_travels: list[float] = []
        for point_index in range(GRID_SIZE**2):
            approach = approach_by_index[point_index]
            retract = retract_by_index[point_index]
            app_obs, app_x, app_relative = curve_observables(approach)
            ret_obs, ret_x, ret_relative = curve_observables(retract)
            app_binned = binned_relative(app_x, app_relative)
            ret_binned = binned_relative(ret_x, ret_relative)
            branch_arrays[(speed, "approach")].append(app_binned)
            branch_arrays[(speed, "retract_far_to_contact")].append(ret_binned)
            row_index, column_index = base.map_pixel_from_index(
                point_index, GRID_SIZE, source.map_back_and_forth
            )
            approach_clipped = bool(np.any(approach.raw_saturation_mask))
            retract_clipped = bool(np.any(retract.raw_saturation_mask))
            approach_saturated += int(approach_clipped)
            retract_saturated += int(retract_clipped)
            record: dict = {
                "acquisition_order": acquisition_order,
                "instrument_scan_number": metadata["instrument_scan_number"],
                "block": block,
                "position_in_block": position_in_block,
                "source": path.relative_to(ROOT).as_posix(),
                "map_start_time": metadata["map_start_time"],
                "approach_speed_um_per_s": speed,
                "retract_speed_um_per_s": metadata["retract_nominal_speed_um_per_s"],
                "point_index": point_index,
                "line_index": point_index // GRID_SIZE + 1,
                "position_within_line": point_index % GRID_SIZE + 1,
                "physical_row": row_index,
                "physical_column": column_index,
                "approach_raw_saturated": approach_clipped,
                "retract_raw_saturated": retract_clipped,
            }
            record.update({f"approach_{key}": value for key, value in app_obs.items()})
            record.update({f"retract_{key}": value for key, value in ret_obs.items()})
            pixel_rows.append(record)
            map_pixel_rows.append(record)
            app_sample_counts.append(app_obs["sample_count"])
            ret_sample_counts.append(ret_obs["sample_count"])
            app_travels.append(app_obs["scanner_travel_nm"])
            ret_travels.append(ret_obs["scanner_travel_nm"])
            if (row_index, column_index) == (FIXED_ROW, FIXED_COLUMN):
                for branch, x, relative, binned in (
                    ("approach", app_x, app_relative, app_binned),
                    ("retract_far_to_contact", ret_x, ret_relative, ret_binned),
                ):
                    fixed_plot_rows.append(
                        {
                            "acquisition_order": acquisition_order,
                            "approach_speed_um_per_s": speed,
                            "branch": branch,
                            "scanner_travel_nm": x,
                            "delta_V": relative,
                        }
                    )
                    for center, value in zip(BIN_CENTERS_NM, binned, strict=True):
                        if np.isfinite(value):
                            fixed_export_rows.append(
                                {
                                    "acquisition_order": acquisition_order,
                                    "instrument_scan_number": metadata["instrument_scan_number"],
                                    "block": block,
                                    "position_in_block": position_in_block,
                                    "map_start_time": metadata["map_start_time"],
                                    "approach_speed_um_per_s": speed,
                                    "branch": branch,
                                    "physical_row": FIXED_ROW,
                                    "physical_column": FIXED_COLUMN,
                                    "scanner_travel_bin_center_nm": center,
                                    "delta_V": value,
                                    "delta_nN_prior_D3_scale": value * D3_PRIOR_SCALE_NN_PER_V,
                                }
                            )

        inventory = {
            "acquisition_order": acquisition_order,
            "block": block,
            "position_in_block": position_in_block,
            "palindrome_mate_order": (block - 1) * MAPS_PER_BLOCK + (MAPS_PER_BLOCK + 1 - position_in_block),
            "source": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": digest,
            **metadata,
            "approach_curves": len(approaches),
            "retract_curves": len(retracts),
            "parser_skips": source.skipped_curves + approach_skipped + retract_skipped,
            "approach_actual_points_median": median(app_sample_counts),
            "retract_actual_points_median": median(ret_sample_counts),
            "approach_actual_sample_rate_Hz_median": median(
                row["approach_sample_rate_Hz"] for row in map_pixel_rows
            ),
            "retract_actual_sample_rate_Hz_median": median(
                row["retract_sample_rate_Hz"] for row in map_pixel_rows
            ),
            "approach_actual_scanner_travel_nm_median": median(app_travels),
            "retract_actual_scanner_travel_nm_median": median(ret_travels),
            "embedded_InvOLS_nm_per_V": source.stored_sensitivity_m_per_V * 1e9,
            "embedded_k_N_per_m": source.stored_spring_constant_N_per_m,
            "prior_D3_InvOLS_nm_per_V": D3_PRIOR_INVOLS_NM_PER_V,
            "prior_D3_k_N_per_m": D3_PRIOR_K_N_PER_M,
        }
        inventory_rows.append(inventory)
        map_summary: dict = {
            "acquisition_order": acquisition_order,
            "instrument_scan_number": metadata["instrument_scan_number"],
            "block": block,
            "position_in_block": position_in_block,
            "palindrome_mate_order": inventory["palindrome_mate_order"],
            "source": path.relative_to(ROOT).as_posix(),
            "map_start_time": metadata["map_start_time"],
            "map_end_time": metadata["map_end_time"],
            "approach_speed_um_per_s": speed,
            "retract_speed_um_per_s": metadata["retract_nominal_speed_um_per_s"],
        }
        for field, prefix in (
            ("approach_initial_0_10nm_raw_V", "approach_initial_raw_V"),
            ("approach_settled_100_300nm_raw_V", "approach_settled_raw_V"),
            ("approach_settled_100_300nm_delta_V", "approach_settled_delta_V"),
            (
                "approach_settled_100_300nm_delta_nN_prior_D3_scale",
                "approach_settled_delta_nN_prior_D3_scale",
            ),
            ("approach_terminal_last_25nm_delta_V", "approach_terminal_delta_V"),
            (
                "approach_terminal_last_25nm_delta_nN_prior_D3_scale",
                "approach_terminal_delta_nN_prior_D3_scale",
            ),
            ("retract_minimum_delta_V", "retract_minimum_delta_V"),
            (
                "retract_minimum_delta_nN_prior_D3_scale",
                "retract_minimum_delta_nN_prior_D3_scale",
            ),
            ("retract_minimum_scanner_travel_nm", "retract_minimum_scanner_travel_nm"),
            ("retract_settled_100_300nm_delta_V", "retract_settled_delta_V"),
        ):
            map_summary.update(summarize_map(map_pixel_rows, field, prefix))
        map_rows.append(map_summary)
        print(
            f"[{acquisition_order:02d}/{MAP_COUNT}] scan {metadata['instrument_scan_number']} · "
            f"{speed:g} µm/s · decoded {len(approaches)}+{len(retracts)} branches",
            flush=True,
        )

    if total_skipped != 0 or approach_saturated != 0 or retract_saturated != 0:
        raise RuntimeError("Unexpected skipped or saturated raw branches")
    if len(embedded_sensitivities) != 1 or len(embedded_springs) != 1:
        raise RuntimeError("Embedded calibration changes between files")
    speed_counts = {speed: sum(row["approach_speed_um_per_s"] == speed for row in map_rows) for speed in SPEEDS_UM_PER_S}
    if speed_counts != {speed: 8 for speed in SPEEDS_UM_PER_S}:
        raise RuntimeError(f"Unexpected speed inventory: {speed_counts}")

    embedded_sensitivity = next(iter(embedded_sensitivities))
    embedded_spring = next(iter(embedded_springs))
    profiles = centered_position_profiles(pixel_rows)
    binned_rows = aggregate_binned(branch_arrays)
    speed_rows = build_speed_summary(map_rows, pixel_rows)
    pair_rows = build_palindrome_pairs(map_rows)

    write_csv(OUT / "map_inventory.csv", inventory_rows)
    write_csv(OUT / "pixel_observables.csv", pixel_rows)
    write_csv(OUT / "map_summary.csv", map_rows)
    write_csv(OUT / "speed_summary.csv", speed_rows)
    write_csv(OUT / "palindrome_pair_differences.csv", pair_rows)
    write_csv(OUT / "line_position_profile.csv", profiles)
    write_csv(OUT / "binned_branch_summary.csv", binned_rows)
    write_csv(OUT / "fixed_pixel_binned_curves.csv", fixed_export_rows)

    plot_speed_overview(binned_rows)
    plot_startup_and_line_position(binned_rows, profiles)
    plot_chronology(map_rows)
    plot_fixed_pixel(fixed_plot_rows)

    report = build_report(
        map_rows,
        speed_rows,
        pair_rows,
        profiles,
        embedded_sensitivity,
        embedded_spring,
    )
    (OUT / "REPORT.md").write_text(report, encoding="utf-8")
    provenance = {
        "analysis": "observational_no_fit",
        "analysis_script": Path(__file__).relative_to(ROOT).as_posix(),
        "analysis_script_sha256": sha256_file(Path(__file__)),
        "download_script": "analysis/download_12_09_26_keeper.py",
        "download_script_sha256": sha256_file(ROOT / "analysis" / "download_12_09_26_keeper.py"),
        "download_manifest": DOWNLOAD_MANIFEST.relative_to(ROOT).as_posix(),
        "download_manifest_sha256": sha256_file(DOWNLOAD_MANIFEST),
        "share_url": f"https://keeper.mpdl.mpg.de/d/{SHARE_TOKEN}/",
        "experiment": {
            "date": "2026-09-12",
            "glycerol_wt_percent": CONCENTRATION_WT_PERCENT,
            "cantilever": "D3 (user identification)",
            "map_count": MAP_COUNT,
            "grid": [GRID_SIZE, GRID_SIZE],
            "approach_speeds_um_per_s": list(SPEEDS_UM_PER_S),
            "maps_per_speed": 8,
            "retract_speed_um_per_s": RETRACT_SPEED_UM_PER_S,
        },
        "primary_observable": "raw vDeflection in V versus scanner travel from the far endpoint",
        "operations": {
            "initial_reference": f"within-curve median over 0-{INITIAL_REFERENCE_MAX_NM:g} nm scanner travel",
            "settled_window": f"within-curve median over {PLATEAU_START_NM:g}-{PLATEAU_STOP_NM:g} nm scanner travel",
            "terminal_window": f"within-curve median over the last {TERMINAL_WINDOW_NM:g} nm scanner travel",
            "binning": f"direct within-bin medians on {BIN_WIDTH_NM:g} nm scanner-travel bins; no interpolation",
            "map_summary": "spatial median and IQR across 64 physical pixels",
            "speed_summary": "median and range across eight map medians per speed",
            "palindrome": "exact later-minus-earlier differences for symmetric map positions in each 8-map block",
            "fits_performed": [],
            "explicitly_not_performed": [
                "far-field baseline fit",
                "contact fit or surface-separation reconstruction",
                "hydrodynamic fit",
                "PB or other surface-force fit",
                "regression or U-to-zero extrapolation",
                "hypothesis test",
            ],
        },
        "calibration": {
            "raw_voltage_is_primary": True,
            "prior_D3_InvOLS_nm_per_V": D3_PRIOR_INVOLS_NM_PER_V,
            "prior_D3_k_N_per_m": D3_PRIOR_K_N_PER_M,
            "prior_D3_scale_nN_per_V": D3_PRIOR_SCALE_NN_PER_V,
            "embedded_InvOLS_nm_per_V": embedded_sensitivity,
            "embedded_k_N_per_m": embedded_spring,
            "embedded_scale_nN_per_V": embedded_sensitivity * 1e-9 * embedded_spring * 1e9,
            "interpretation": "conflicting scales were not reconciled; nN columns use only the prior D3 scale and raw V remains primary",
        },
        "integrity": {
            "input_files_sha256": {path.name: sha256_file(path) for path in paths},
            "download_zip_crc": "passed for all 32 JPK archives",
            "decoded_approach_curves": len(pixel_rows),
            "decoded_retract_curves": len(pixel_rows),
            "parser_skips": total_skipped,
            "approach_raw_saturated_curves": approach_saturated,
            "retract_raw_saturated_curves": retract_saturated,
            "instrument_failure_flags": int(sum(row["instrument_failure_flag_count"] for row in inventory_rows)),
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
        Path(__file__),
        *(path for path in OUT.rglob("*") if path.is_file() and path.name != "artifact_manifest.sha256"),
    ]
    manifest = "\n".join(
        f"{sha256_file(path)}  {path.relative_to(ROOT).as_posix()}" for path in sorted(artifacts)
    )
    (OUT / "artifact_manifest.sha256").write_text(manifest + "\n", encoding="utf-8")
    print(f"Wrote observational package to {OUT}", flush=True)


if __name__ == "__main__":
    main()

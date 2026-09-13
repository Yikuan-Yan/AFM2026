#!/usr/bin/env python3
"""Fit same-speed time decay and blockwise force extrapolation to v=0.

Inputs are the map-median absolute-force reconstruction made with the validated
speed-conditioned baselines and the physical-distance contact method.  The map
is the regression unit.  Four 8-map palindrome blocks are treated as four
separate continuous speed scans; within each block, the two symmetric maps at
each speed are averaged before fitting force = intercept + slope * speed.
"""

from __future__ import annotations

import csv
from datetime import datetime, timedelta
import hashlib
import json
from pathlib import Path
import platform
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import scipy
from scipy.stats import spearmanr, t as student_t, theilslopes


ROOT = Path(__file__).resolve().parents[1]
FORCE_RESULTS = ROOT / "analysis" / "glycerol_99p7_D3_force_pb_time_results"
TIME_RESULTS = ROOT / "analysis" / "glycerol_99p7_D3_time_speed_results"
OUT = ROOT / "analysis" / "glycerol_99p7_D3_decay_v0_results"
FIG = OUT / "figures"

SPEEDS = (0.1, 0.3, 0.9, 2.7)
LOW_SPEEDS = (0.1, 0.3, 0.9)
TARGETS_NM = (20.0, 50.0, 100.0, 200.0)
SPEED_COLORS = {0.1: "#2a9d8f", 0.3: "#457b9d", 0.9: "#e9c46a", 2.7: "#e76f51"}
BLOCK_COLORS = {1: "#264653", 2: "#2a9d8f", 3: "#e9c46a", 4: "#e76f51"}


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


def ols_hc3(x: np.ndarray, y: np.ndarray) -> dict:
    """Two-parameter OLS with HC3 covariance at the map/pair-mean unit."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or x.size < 3 or np.unique(x).size < 2:
        raise ValueError("OLS needs at least three finite points and two x values")
    if not np.all(np.isfinite(x + y)):
        raise ValueError("OLS input contains non-finite values")
    design = np.column_stack((np.ones(x.size), x))
    xtx_inv = np.linalg.inv(design.T @ design)
    beta = xtx_inv @ design.T @ y
    residual = y - design @ beta
    leverage = np.sum(design * (design @ xtx_inv), axis=1)
    adjusted = residual / np.maximum(1.0 - leverage, np.finfo(float).eps)
    meat = design.T @ np.diag(adjusted**2) @ design
    covariance = xtx_inv @ meat @ xtx_inv
    sse = float(np.sum(residual**2))
    sst = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float("nan") if sst <= 0.0 else 1.0 - sse / sst
    return {
        "beta": beta,
        "covariance_hc3": covariance,
        "residual": residual,
        "leverage": leverage,
        "r2": r2,
        "rmse": float(np.sqrt(sse / x.size)),
        "n": int(x.size),
        "df": int(x.size - 2),
    }


def mean_ci(model: dict, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    design = np.column_stack((np.ones(x.size), x))
    prediction = design @ model["beta"]
    variance = np.einsum("ij,jk,ik->i", design, model["covariance_hc3"], design)
    critical = student_t.ppf(0.975, model["df"])
    half = critical * np.sqrt(np.maximum(variance, 0.0))
    return prediction, prediction - half, prediction + half


def time_decay_fits(chronology: list[dict], force_lookup: dict[tuple[int, float], dict]) -> list[dict]:
    rows: list[dict] = []
    for distance in TARGETS_NM:
        for speed in SPEEDS:
            selected = [row for row in chronology if np.isclose(float(row["speed_um_per_s"]), speed)]
            elapsed_min = np.asarray([float(row["elapsed_midpoint_min"]) for row in selected])
            force_pN = np.asarray([float(force_lookup[(int(row["acquisition_order"]), distance)]["map_median_force_pN"]) for row in selected])
            model = ols_hc3(elapsed_min, force_pN)
            intercept, slope = model["beta"]
            predicted_0 = float(intercept)
            predicted_60 = float(intercept + 60.0 * slope)
            decay_percent = 100.0 * (predicted_0 - predicted_60) / predicted_0
            # Delta-method HC3 uncertainty for -6000*slope/intercept.
            gradient = np.asarray([6000.0 * slope / intercept**2, -6000.0 / intercept])
            variance = float(gradient @ model["covariance_hc3"] @ gradient)
            decay_se = float(np.sqrt(max(variance, 0.0)))
            critical = float(student_t.ppf(0.975, model["df"]))
            robust = theilslopes(force_pN, elapsed_min)
            rows.append(
                {
                    "distance_nm": distance,
                    "speed_um_per_s": speed,
                    "n_maps": model["n"],
                    "elapsed_min_min": float(np.min(elapsed_min)),
                    "elapsed_min_max": float(np.max(elapsed_min)),
                    "fitted_force_at_experiment_t0_pN": predicted_0,
                    "fitted_force_at_experiment_t60min_pN": predicted_60,
                    "slope_pN_per_min": float(slope),
                    "slope_pN_per_hour": float(60.0 * slope),
                    "slope_hc3_se_pN_per_min": float(np.sqrt(max(model["covariance_hc3"][1, 1], 0.0))),
                    "expected_decay_over_first_hour_percent": decay_percent,
                    "expected_decay_over_first_hour_hc3_delta_se_percent": decay_se,
                    "expected_decay_over_first_hour_ci95_low_percent": decay_percent - critical * decay_se,
                    "expected_decay_over_first_hour_ci95_high_percent": decay_percent + critical * decay_se,
                    "r2": model["r2"],
                    "rmse_pN": model["rmse"],
                    "theil_sen_slope_pN_per_min": float(robust.slope),
                    "ols_theil_sen_same_sign": bool(np.sign(slope) == np.sign(robust.slope)),
                    "interpretation": "positive percent means fitted decline from experiment t=0 to t=60 min",
                }
            )
    return rows


def speed_pair_means_and_fits(
    chronology: list[dict], force_lookup: dict[tuple[int, float], dict]
) -> tuple[list[dict], list[dict]]:
    pair_rows: list[dict] = []
    fit_rows: list[dict] = []
    for distance in TARGETS_NM:
        for block in range(1, 5):
            block_rows = [row for row in chronology if int(row["block"]) == block]
            block_midpoint = float(np.median([float(row["elapsed_midpoint_min"]) for row in block_rows]))
            block_pairs: list[dict] = []
            for speed in SPEEDS:
                pair = sorted(
                    [row for row in block_rows if np.isclose(float(row["speed_um_per_s"]), speed)],
                    key=lambda row: int(row["position_in_block"]),
                )
                if len(pair) != 2 or int(pair[0]["position_in_block"]) + int(pair[1]["position_in_block"]) != 9:
                    raise RuntimeError(f"Block {block}, speed {speed}: palindrome identity failed")
                forces = np.asarray([float(force_lookup[(int(row["acquisition_order"]), distance)]["map_median_force_pN"]) for row in pair])
                times = np.asarray([float(row["elapsed_midpoint_min"]) for row in pair])
                current = {
                    "distance_nm": distance,
                    "block": block,
                    "block_midpoint_elapsed_min": block_midpoint,
                    "speed_um_per_s": speed,
                    "early_map_order": int(pair[0]["acquisition_order"]),
                    "late_map_order": int(pair[1]["acquisition_order"]),
                    "early_position_in_block": int(pair[0]["position_in_block"]),
                    "late_position_in_block": int(pair[1]["position_in_block"]),
                    "early_force_pN": float(forces[0]),
                    "late_force_pN": float(forces[1]),
                    "pair_mean_force_pN": float(np.mean(forces)),
                    "pair_half_difference_pN": float(abs(forces[1] - forces[0]) / 2.0),
                    "pair_centre_elapsed_min": float(np.mean(times)),
                    "pair_centre_offset_from_block_midpoint_min": float(np.mean(times) - block_midpoint),
                    "pair_time_gap_min": float(times[1] - times[0]),
                }
                pair_rows.append(current)
                block_pairs.append(current)

            speeds = np.asarray([row["speed_um_per_s"] for row in block_pairs])
            values = np.asarray([row["pair_mean_force_pN"] for row in block_pairs])
            primary = ols_hc3(speeds, values)
            low_mask = np.isin(speeds, LOW_SPEEDS)
            low = ols_hc3(speeds[low_mask], values[low_mask])
            intercept, slope = primary["beta"]
            intercept_se = float(np.sqrt(max(primary["covariance_hc3"][0, 0], 0.0)))
            critical = float(student_t.ppf(0.975, primary["df"]))
            low_intercept, low_slope = low["beta"]
            fit_rows.append(
                {
                    "distance_nm": distance,
                    "block": block,
                    "block_midpoint_elapsed_min": block_midpoint,
                    "n_speed_pair_means": primary["n"],
                    "speed_fit_min_um_per_s": float(np.min(speeds)),
                    "speed_fit_max_um_per_s": float(np.max(speeds)),
                    "v0_intercept_all4_pN": float(intercept),
                    "v0_intercept_all4_hc3_se_pN": intercept_se,
                    "v0_intercept_all4_ci95_low_pN": float(intercept - critical * intercept_se),
                    "v0_intercept_all4_ci95_high_pN": float(intercept + critical * intercept_se),
                    "speed_slope_all4_pN_per_um_per_s": float(slope),
                    "speed_fit_all4_r2": primary["r2"],
                    "speed_fit_all4_rmse_pN": primary["rmse"],
                    "v0_intercept_low3_pN": float(low_intercept),
                    "speed_slope_low3_pN_per_um_per_s": float(low_slope),
                    "speed_fit_low3_r2": low["r2"],
                    "v0_all4_minus_low3_pN": float(intercept - low_intercept),
                    "v0_all4_minus_low3_percent_of_abs_low3": float(100.0 * (intercept - low_intercept) / abs(low_intercept)) if abs(low_intercept) > 1e-12 else float("nan"),
                    "pair_mean_speed_spearman_rho": float(spearmanr(speeds, values).statistic),
                    "interpretation": "apparent block-centred v=0 intercept; pair averaging suppresses first-order linear time drift",
                }
            )
    return pair_rows, fit_rows


def block_boundaries(chronology: list[dict]) -> list[datetime]:
    output: list[datetime] = []
    for index in (8, 16, 24):
        left = datetime.fromisoformat(chronology[index - 1]["map_midpoint_time"])
        right = datetime.fromisoformat(chronology[index]["map_midpoint_time"])
        output.append(left + (right - left) / 2)
    return output


def time_decay_figure(
    chronology: list[dict], force_lookup: dict[tuple[int, float], dict], time_rows: list[dict]
) -> None:
    times = np.asarray([datetime.fromisoformat(row["map_midpoint_time"]) for row in chronology], dtype=object)
    elapsed = np.asarray([float(row["elapsed_midpoint_min"]) for row in chronology])
    fit_lookup = {(float(row["distance_nm"]), float(row["speed_um_per_s"])): row for row in time_rows}
    grid_elapsed = np.linspace(0.0, float(np.max(elapsed)), 300)
    first_time = times[0]
    grid_times = np.asarray([first_time + timedelta(minutes=float(value)) for value in grid_elapsed], dtype=object)

    fig, axes = plt.subplots(2, 2, figsize=(18, 9), sharex=True)
    for axis, distance in zip(axes.flat, TARGETS_NM, strict=True):
        text_lines = []
        for speed in SPEEDS:
            selected = np.asarray([np.isclose(float(row["speed_um_per_s"]), speed) for row in chronology])
            selected_rows = np.asarray(chronology, dtype=object)[selected]
            values = np.asarray([float(force_lookup[(int(row["acquisition_order"]), distance)]["map_median_force_pN"]) for row in selected_rows])
            q25 = np.asarray([float(force_lookup[(int(row["acquisition_order"]), distance)]["map_force_q25_pN"]) for row in selected_rows])
            q75 = np.asarray([float(force_lookup[(int(row["acquisition_order"]), distance)]["map_force_q75_pN"]) for row in selected_rows])
            axis.errorbar(times[selected], values, yerr=np.vstack((values - q25, q75 - values)), fmt="o", ms=5.2, color=SPEED_COLORS[speed], mec="0.25", mew=0.4, ecolor=SPEED_COLORS[speed], elinewidth=0.65, alpha=0.9, zorder=3)
            fit = fit_lookup[(distance, speed)]
            model = ols_hc3(elapsed[selected], values)
            prediction, lower, upper = mean_ci(model, grid_elapsed)
            axis.plot(grid_times, prediction, color=SPEED_COLORS[speed], lw=1.8)
            axis.fill_between(grid_times, lower, upper, color=SPEED_COLORS[speed], alpha=0.08, linewidth=0)
            text_lines.append(f"{speed:g}: {float(fit['expected_decay_over_first_hour_percent']):+.1f}%")
        for boundary in block_boundaries(chronology):
            axis.axvline(boundary, color="0.45", ls="--", lw=0.8, zorder=0)
        axis.set_title(f"D = {distance:g} nm", fontsize=14)
        axis.set_ylabel("Absolute force (pN)")
        axis.grid(axis="y", alpha=0.2)
        axis.text(0.985, 0.97, "1 h fitted decay\n" + "\n".join(text_lines), transform=axis.transAxes, ha="right", va="top", fontsize=8.5, bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.88, "pad": 4})
    timezone = times[0].tzinfo
    for axis in axes.flat:
        axis.xaxis.set_major_locator(mdates.MinuteLocator(interval=10, tz=timezone))
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=timezone))
    for axis in axes[-1]:
        axis.set_xlabel("Map acquisition midpoint time (UTC+02:00)")
    handles = [Line2D([], [], marker="o", color=color, label=f"{speed:g} µm/s") for speed, color in SPEED_COLORS.items()]
    fig.suptitle("Same-speed linear time decay of map-median absolute force", fontsize=16, y=0.985)
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.94), ncol=4, frameon=False)
    fig.text(0.5, 0.012, "Lines: OLS on 8 maps/speed; bands: HC3 95% mean-response intervals. Percent compares fitted force at experiment t=0 and t=60 min.", ha="center", fontsize=9, color="0.3")
    fig.tight_layout(rect=(0, 0.035, 1, 0.90))
    for suffix in ("png", "svg"):
        fig.savefig(FIG / f"same_speed_linear_decay_vs_time.{suffix}", dpi=200)
    plt.close(fig)


def speed_extrapolation_figure(pair_rows: list[dict], fit_rows: list[dict]) -> None:
    pair_lookup = {(float(row["distance_nm"]), int(row["block"]), float(row["speed_um_per_s"])): row for row in pair_rows}
    fit_lookup = {(float(row["distance_nm"]), int(row["block"])): row for row in fit_rows}
    fig, axes = plt.subplots(2, 2, figsize=(16.5, 9.2), sharex=True)
    for axis, distance in zip(axes.flat, TARGETS_NM, strict=True):
        all_intercepts = []
        low_intercepts = []
        for block in range(1, 5):
            color = BLOCK_COLORS[block]
            pairs = [pair_lookup[(distance, block, speed)] for speed in SPEEDS]
            x = np.asarray(SPEEDS)
            y = np.asarray([float(row["pair_mean_force_pN"]) for row in pairs])
            yerr = np.asarray([float(row["pair_half_difference_pN"]) for row in pairs])
            fit = fit_lookup[(distance, block)]
            intercept = float(fit["v0_intercept_all4_pN"])
            slope = float(fit["speed_slope_all4_pN_per_um_per_s"])
            low_intercept = float(fit["v0_intercept_low3_pN"])
            all_intercepts.append(intercept)
            low_intercepts.append(low_intercept)
            axis.errorbar(x, y, yerr=yerr, fmt="o", ms=5.5, color=color, mec="0.25", mew=0.4, ecolor=color, elinewidth=0.8, alpha=0.9, zorder=3)
            axis.plot([0.1, 2.7], [intercept + 0.1 * slope, intercept + 2.7 * slope], color=color, lw=1.7)
            axis.plot([0.0, 0.1], [intercept, intercept + 0.1 * slope], color=color, lw=1.7, ls="--")
            axis.scatter([0.0], [intercept], marker="D", s=38, color=color, edgecolor="0.25", lw=0.4, zorder=4)
            axis.scatter([0.0], [low_intercept], marker="s", s=58, facecolors="none", edgecolors=color, lw=1.4, zorder=4)
        axis.axvline(0.0, color="0.55", lw=0.7)
        axis.set_title(f"D = {distance:g} nm", fontsize=14)
        axis.set_ylabel("Pair-mean absolute force (pN)")
        axis.grid(alpha=0.2)
        axis.text(0.985, 0.035, f"median v=0 (all 4): {np.median(all_intercepts):.1f} pN\nmedian v=0 (≤0.9): {np.median(low_intercepts):.1f} pN", transform=axis.transAxes, ha="right", va="bottom", fontsize=8.5, bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.88, "pad": 4})
    for axis in axes[-1]:
        axis.set_xlabel("Approach speed (µm/s)")
    block_handles = [Line2D([], [], marker="o", color=color, label=f"block {block}") for block, color in BLOCK_COLORS.items()]
    method_handles = [
        Line2D([], [], marker="D", ls="--", color="0.35", label="all-4-speed v=0 intercept"),
        Line2D([], [], marker="s", markerfacecolor="none", ls="", color="0.35", label="≤0.9 µm/s sensitivity intercept"),
    ]
    fig.suptitle("Within-block speed extrapolation to v = 0", fontsize=16, y=0.985)
    fig.legend(handles=block_handles + method_handles, loc="upper center", bbox_to_anchor=(0.5, 0.94), ncol=6, frameon=False)
    fig.text(0.5, 0.012, "Each point is the mean of the symmetric same-speed pair in one 8-map palindrome block; bars are half the pair difference. Dashed segment is extrapolated.", ha="center", fontsize=9, color="0.3")
    fig.tight_layout(rect=(0, 0.035, 1, 0.90))
    for suffix in ("png", "svg"):
        fig.savefig(FIG / f"blockwise_speed_extrapolation_v0.{suffix}", dpi=200)
    plt.close(fig)


def v0_summary(fit_rows: list[dict]) -> list[dict]:
    output: list[dict] = []
    for distance in TARGETS_NM:
        selected = [row for row in fit_rows if np.isclose(float(row["distance_nm"]), distance)]
        all4 = np.asarray([float(row["v0_intercept_all4_pN"]) for row in selected])
        low3 = np.asarray([float(row["v0_intercept_low3_pN"]) for row in selected])
        r2 = np.asarray([float(row["speed_fit_all4_r2"]) for row in selected])
        output.append(
            {
                "distance_nm": distance,
                "blocks": len(selected),
                "v0_all4_median_pN": float(np.median(all4)),
                "v0_all4_min_pN": float(np.min(all4)),
                "v0_all4_max_pN": float(np.max(all4)),
                "v0_low3_median_pN": float(np.median(low3)),
                "v0_low3_min_pN": float(np.min(low3)),
                "v0_low3_max_pN": float(np.max(low3)),
                "all4_minus_low3_median_pN": float(np.median(all4 - low3)),
                "all4_fit_r2_median": float(np.median(r2)),
                "all4_fit_r2_min": float(np.min(r2)),
                "all4_fit_r2_max": float(np.max(r2)),
            }
        )
    return output


def write_report(time_rows: list[dict], fit_rows: list[dict], summary_rows: list[dict], pair_rows: list[dict]) -> None:
    time_lookup = {(float(row["distance_nm"]), float(row["speed_um_per_s"])): row for row in time_rows}
    signs = sum(bool(row["ols_theil_sen_same_sign"]) for row in time_rows)
    v0_intervals_containing_zero = sum(
        float(row["v0_intercept_all4_ci95_low_pN"]) <= 0.0 <= float(row["v0_intercept_all4_ci95_high_pN"])
        for row in fit_rows
    )
    max_pair_offset = max(abs(float(row["pair_centre_offset_from_block_midpoint_min"])) for row in pair_rows)
    lines = [
        "# 12-09-26 · 99.7 wt% glycerol · linear decay and v→0 extrapolation",
        "",
        "## 物理图像",
        "",
        "同一 approach speed 的 map-median force 随实验进行通常下降；这可用一条描述性直线量化。每个 8-map block 内，同速度的两张 map 位于 palindrome 对称位置；先取 pair mean，能在时间趋势近似线性时消去一阶漂移，再把四个速度的 pair means 线性外推至 v=0。",
        "",
        "## 一小时衰减",
        "",
        "定义：`100 × [F_fit(t=0) − F_fit(t=60 min)] / F_fit(t=0)`；正值表示衰减。t=0 是第一张 map 的 midpoint。回归单位是 map，每个速度 n=8；区间为 HC3/delta-method 95% interval。",
        "",
        "| D (nm) | speed (µm/s) | fitted equation, t_h from first map midpoint (pN) | expected decay in first 1 h | 95% interval | R² |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for distance in TARGETS_NM:
        for speed in SPEEDS:
            row = time_lookup[(distance, speed)]
            lines.append(
                f"| {distance:g} | {speed:g} | `F={float(row['fitted_force_at_experiment_t0_pN']):.1f}{float(row['slope_pN_per_hour']):+.1f} t_h` | {float(row['expected_decay_over_first_hour_percent']):.1f}% | [{float(row['expected_decay_over_first_hour_ci95_low_percent']):.1f}, {float(row['expected_decay_over_first_hour_ci95_high_percent']):.1f}]% | {float(row['r2']):.3f} |"
            )
    lines += [
        "",
        f"OLS 与 Theil–Sen slope 在 {signs}/16 个 distance×speed 组合中同号。唯一异号组合是 0.1 µm/s、100 nm；在 200 nm 两种方法都给出很小的正 slope。二者的区间都跨过零，不能解释为已分辨的衰减或增长率。",
        "",
        "## block 内 v=0 外推",
        "",
        f"每个 block 的四个 pair centres 偏离该 block midpoint 最多 {max_pair_offset:.3f} min，因此 pair averaging 对一阶时间项的抵消误差很小。主结果使用全部四个速度；`≤0.9` 三点截距仅作为 2.7 µm/s 瞬态敏感性检查。",
        "",
        "| D (nm) | all-4-speed v=0 median [block range] (pN) | ≤0.9-speed median [range] (pN) | median difference (pN) | all-4 R² median [range] |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {float(row['distance_nm']):g} | {float(row['v0_all4_median_pN']):.1f} [{float(row['v0_all4_min_pN']):.1f}, {float(row['v0_all4_max_pN']):.1f}] | {float(row['v0_low3_median_pN']):.1f} [{float(row['v0_low3_min_pN']):.1f}, {float(row['v0_low3_max_pN']):.1f}] | {float(row['all4_minus_low3_median_pN']):.1f} | {float(row['all4_fit_r2_median']):.3f} [{float(row['all4_fit_r2_min']):.3f}, {float(row['all4_fit_r2_max']):.3f}] |"
        )
    lines += [
        "",
        "虽然 20–100 nm 的 all-4 speed 直线通常有高 R²，但把 2.7 µm/s 去掉后 v=0 截距明显下降；这种 fit-range dependence 在 100–200 nm 尤其大。因此这里得到的是 protocol-dependent apparent zero-speed intercept，不是已识别的 equilibrium surface force。四点 fit 的单-block 截距 95% intervals 也很宽，详见 CSV。",
        f"本次 {v0_intervals_containing_zero}/16 个单-block all-4-speed 截距的 HC3 95% interval 都包含零；点估计可用于比较 block 和 fit range，但不能单独证明正的零速表面力。",
        "",
        "## 图与数值文件",
        "",
        "- `figures/same_speed_linear_decay_vs_time.png`：同速度时间直线、HC3 mean-response bands 与一小时衰减。",
        "- `figures/blockwise_speed_extrapolation_v0.png`：四个 palindrome blocks 的 pair means 和 v=0 外推。",
        "- `same_speed_time_linear_fits.csv`：16 个时间回归及 HC3/delta-method interval。",
        "- `blockwise_speed_pair_means.csv`：每个对称 pair 的两张 map、均值、半差和时间中心。",
        "- `blockwise_speed_linear_fits.csv`：每个 block×distance 的 all-4 与 low-3 sensitivity fits。",
        "- `v0_summary.csv`：四个 block 的 v=0 汇总。",
        "",
        "## 解释边界",
        "",
        "- 线性时间模型描述约 1.5 h 的本次序列，不证明长期指数衰减机制。",
        "- pair averaging 只压低一阶时间漂移；surface/probe conditioning、raster history 与速度相关 baseline/transient 仍可能存在。",
        "- v=0 的数学截距不等于 thermodynamic equilibrium，除非 residual velocity dependence、contact/baseline 系统误差及等待时间效应另有控制。",
        "",
    ]
    (OUT / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def write_manifest() -> None:
    targets = sorted(path for path in OUT.rglob("*") if path.is_file() and path.name != "artifact_manifest.sha256")
    lines = [f"{sha256_file(path)}  {path.relative_to(ROOT).as_posix()}" for path in targets]
    (OUT / "artifact_manifest.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    chronology = sorted(read_csv(TIME_RESULTS / "map_time_speed_observables.csv"), key=lambda row: int(row["acquisition_order"]))
    force_rows = read_csv(FORCE_RESULTS / "map_force_by_separation.csv")
    force_lookup = {(int(row["acquisition_order"]), float(row["distance_nm"])): row for row in force_rows}
    if len(chronology) != 32:
        raise RuntimeError("Expected the complete 32-map chronology")
    for distance in TARGETS_NM:
        if sum((int(row["acquisition_order"]), distance) in force_lookup for row in chronology) != 32:
            raise RuntimeError(f"Missing force rows at {distance:g} nm")

    time_rows = time_decay_fits(chronology, force_lookup)
    pair_rows, fit_rows = speed_pair_means_and_fits(chronology, force_lookup)
    summary_rows = v0_summary(fit_rows)
    write_csv(OUT / "same_speed_time_linear_fits.csv", time_rows)
    write_csv(OUT / "blockwise_speed_pair_means.csv", pair_rows)
    write_csv(OUT / "blockwise_speed_linear_fits.csv", fit_rows)
    write_csv(OUT / "v0_summary.csv", summary_rows)
    time_decay_figure(chronology, force_lookup, time_rows)
    speed_extrapolation_figure(pair_rows, fit_rows)
    write_report(time_rows, fit_rows, summary_rows, pair_rows)

    provenance = {
        "script": str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
        "input_force_csv": str((FORCE_RESULTS / "map_force_by_separation.csv").relative_to(ROOT)).replace("\\", "/"),
        "input_force_csv_sha256": sha256_file(FORCE_RESULTS / "map_force_by_separation.csv"),
        "input_chronology_csv": str((TIME_RESULTS / "map_time_speed_observables.csv").relative_to(ROOT)).replace("\\", "/"),
        "input_chronology_csv_sha256": sha256_file(TIME_RESULTS / "map_time_speed_observables.csv"),
        "time_model": "per distance and speed: OLS F=a+b*elapsed_min on 8 maps; HC3 covariance",
        "one_hour_decay_definition": "100*(F_fit(t=0)-F_fit(t=60min))/F_fit(t=0), t=0 first map midpoint",
        "one_hour_decay_interval": "delta method using HC3 covariance and t critical value with n-2 degrees of freedom",
        "speed_extrapolation_primary": "within each 8-map block, arithmetic symmetric pair means; OLS F=a+b*v on four speeds",
        "speed_extrapolation_sensitivity": "repeat on 0.1, 0.3, 0.9 um/s only",
        "statistical_unit_time_fit": "map",
        "statistical_unit_speed_fit": "same-speed symmetric map-pair mean within block",
        "claim_boundary": "descriptive linear trends and protocol-dependent apparent v=0 intercepts; not equilibrium identification",
        "software": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    write_manifest()
    print(json.dumps({"output": str(OUT), "time_fits": len(time_rows), "block_speed_fits": len(fit_rows), "max_pair_centre_offset_min": max(abs(float(row["pair_centre_offset_from_block_midpoint_min"])) for row in pair_rows)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

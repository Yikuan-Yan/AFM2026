#!/usr/bin/env python3
"""Document blind liquid-refresh candidates from the completed water analysis.

The candidate ranking is an interpretation of measured changes, not an
operation-log reconstruction with validated causal labels or probabilities.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import analyze_08_09_26_fixed_pixel as a


def sampling_diagnostics(maps: list[dict], fixed: list[dict], provenance: dict) -> None:
    """Preserve direct raw evidence for empty bins and the map-25 baseline."""
    sensitivity = provenance["water_InvOLS_nm_per_V"] * 1e-9
    payloads = {}
    missing_rows = [r for r in fixed if r["raw_samples_in_bin"] == 0]
    for number in sorted({24, 25, 26} | {r["acquisition_order"] for r in missing_rows}):
        source = a.base.load_source((a.ROOT / maps[number - 1]["source"]).resolve(), 0)
        curve = a.reference.curve_at_physical_pixel(source)
        far = a.base.fit_far_field_drift(curve.measured_height_m, curve.deflection_V)
        baseline = a.base.baseline_voltage(curve.measured_height_m, curve.deflection_V, far)
        corrected = curve.deflection_V - baseline
        contact = a.pilot.terminal_contact_fit(curve.measured_height_m, corrected, 50)
        displacement = sensitivity * corrected
        distance = (curve.measured_height_m + displacement + contact["intercept_V"] / contact["slope_V_per_m"]) * 1e9
        constant_reference = np.median(curve.deflection_V[:far.n_points])
        constant = a.K * sensitivity * (curve.deflection_V - constant_reference) * 1e12
        fit_force = a.K * sensitivity * (baseline - constant_reference) * 1e12
        payloads[number] = (curve, far, contact, distance, displacement, constant, fit_force)
    gaps = []
    for row in missing_rows:
        number, target = row["acquisition_order"], row["distance_nm"]
        curve, far, contact, distance, displacement, constant, fit_force = payloads[number]
        pre = np.arange(len(distance) - 1) < contact["start"] - 1
        crosses = np.flatnonzero(pre & (distance[:-1] >= target + 2.5) & (distance[1:] < target - 2.5))
        if len(crosses) == 0:
            raise RuntimeError("Empty bin has no verified downward sample jump")
        for i in crosses:
            gaps.append({"acquisition_order": number, "distance_nm": target, "bin_low_nm": target - 2.5, "bin_high_nm": target + 2.5,
                         "sample_index_before": int(i), "sample_index_after": int(i + 1),
                         "distance_before_nm": distance[i], "distance_after_nm": distance[i + 1],
                         "force_before_pN": constant[i] - fit_force[i], "force_after_pN": constant[i + 1] - fit_force[i + 1],
                         "scanner_height_change_nm": (curve.measured_height_m[i + 1] - curve.measured_height_m[i]) * 1e9,
                         "cantilever_displacement_change_nm": (displacement[i + 1] - displacement[i]) * 1e9,
                         "sample_step_ms": curve.duration_s / len(distance) * 1000, "raw_samples_in_bin": 0})
    a.pilot.write_csv(a.OUT / "fixed_pixel_missing_bin_jumps.csv", gaps)
    baseline_rows = []
    for number in (24, 25, 26):
        curve, far, contact, distance, displacement, constant, fit_force = payloads[number]
        halves = [slice(0, far.n_points // 2), slice(far.n_points // 2, far.n_points)]
        slopes = [a.base.robust_line(curve.measured_height_m[part], curve.deflection_V[part])[0] * a.K * sensitivity * 1e5 for part in halves]
        baseline_rows.append({"acquisition_order": number, "point_index": 28, "far_samples": far.n_points,
                              "far_slope_pN_per_100nm_height": a.K * sensitivity * far.slope_V_per_m * 1e5,
                              "far_fit_R2": far.r2, "far_residual_MAD_pN": a.K * sensitivity * far.residual_mad_V * 1e12,
                              "far_first_half_slope_pN_per_100nm": slopes[0], "far_second_half_slope_pN_per_100nm": slopes[1],
                              "contact_R2": contact["r2"], "contact_InvOLS_nm_per_V": contact["sensitivity_m_per_V"] * 1e9})
    a.pilot.write_csv(a.OUT / "map25_fixed_pixel_baseline_QC.csv", baseline_rows)
    curve, far, contact, distance, displacement, constant, fit_force = payloads[25]
    travel = (curve.measured_height_m[0] - curve.measured_height_m) * 1e9
    show = (np.arange(len(distance)) < contact["start"]) & (distance >= 75)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.3))
    ax = axes[0]
    ax.axvspan(travel[0], travel[far.n_points - 1], color="#457b9d", alpha=0.1, label="Far-field fitting window")
    ax.plot(travel[show], constant[show], color="0.45", lw=0.8, label="Raw signal, constant referenced")
    ax.plot(travel[show], fit_force[show], color="#b45309", ls="--", lw=1.3, label="Far line extrapolated")
    ax.set_xlabel("Scanner travel from approach start (nm)")
    ax.set_ylabel("Force equivalent, constant referenced (pN)")
    ax.set_title("Far-field slope extrapolates to a large offset")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.2)
    target = (np.arange(len(distance)) < contact["start"]) & (distance >= 97.5) & (distance < 102.5)
    x = float(np.median(travel[target]))
    lo, hi = float(np.median(fit_force[target])), float(np.median(constant[target]))
    ax.annotate("", xy=(x, hi), xytext=(x, lo), arrowprops={"arrowstyle": "<->", "color": "#9b2226"})
    ax.text(x - 12, (hi + lo) / 2, "~754 pN\nat D = 100 nm", ha="right", va="center", fontsize=9, color="#9b2226")
    ax = axes[1]
    pre = np.arange(len(distance)) < contact["start"]
    line_bins = a.pilot.bin_median(distance, constant - fit_force, pre)
    constant_bins = a.pilot.bin_median(distance, constant, pre)
    displayed_bins = (a.pilot.BIN_CENTERS_NM >= 75) & (a.pilot.BIN_CENTERS_NM <= 220)
    ax.plot(a.pilot.BIN_CENTERS_NM[displayed_bins], line_bins[displayed_bins], color="#9b2226", label="Far-linear: 815.0 pN at 100 nm")
    ax.plot(a.pilot.BIN_CENTERS_NM[displayed_bins], constant_bins[displayed_bins], color="#457b9d", label="Far-constant: 61.2 pN at 100 nm")
    ax.scatter(distance[target], (constant - fit_force)[target], color="#9b2226", s=16)
    ax.scatter(distance[target], constant[target], color="#457b9d", s=16)
    ax.set_xlim(75, 220)
    ax.set_xlabel("Reconstructed separation D (nm)")
    ax.set_ylabel("Force (pN)")
    ax.set_title("The same six samples shift together at 100 nm")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.2)
    fig.suptitle("Map 25 · pixel (3,3) · 1 µm/s · 16:34:23 (UTC+02:00)", fontsize=14)
    fig.text(0.5, 0.015, "The baseline comparison diagnoses extrapolation sensitivity; it does not identify the origin of the raw far-field drift.", ha="center", fontsize=8.5)
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.savefig(a.FIG / "map25_fixed_pixel_baseline_diagnostic.png", dpi=180)
    plt.close(fig)


def main() -> None:
    maps = a.read_derived_csv(a.OUT / "map_inventory_QC.csv")
    pixels = a.read_derived_csv(a.OUT / "pixel_QC.csv")
    fixed = a.read_derived_csv(a.OUT / "fixed_pixel_row3_col3_force_slices.csv")
    adjacent = a.read_derived_csv(a.OUT / "adjacent_map_changes.csv")
    provenance = json.loads((a.OUT / "provenance.json").read_text(encoding="utf-8"))
    sampling_diagnostics(maps, fixed, provenance)
    with np.load(a.OUT / "all_pixel_force_curves.npz", allow_pickle=False) as saved:
        curves = saved["line_pN"].copy()
    by_map = {i: {p["point_index"]: p for p in pixels if p["acquisition_order"] == i} for i in range(1, 37)}
    detected_snap = np.asarray([
        np.median([p["approach_snap_distance_nm"] for p in by_map[i].values()
                   if p["approach_snap_detected"] and np.isfinite(p["approach_snap_distance_nm"])])
        for i in range(1, 37)
    ])
    candidate_rows = []
    for before in range(1, 36):
        after = before + 1
        left, right = maps[before - 1], maps[after - 1]
        delta_s = np.asarray([by_map[after][p]["contact_invOLS_nm_per_V"] - by_map[before][p]["contact_invOLS_nm_per_V"] for p in range(64)])
        delta_f = next(r for r in adjacent if r["before_map"] == before and r["distance_nm"] == 20 and r["baseline"] == "line")
        delta_constant = next(r for r in adjacent if r["before_map"] == before and r["distance_nm"] == 20 and r["baseline"] == "constant")
        candidate_rows.append({
            "before_map": before, "after_map": after,
            "idle_interval_start_instrument_time": left["map_end_time"],
            "idle_interval_end_instrument_time": right["map_start_time"],
            "idle_s": right["idle_before_map_s"],
            "before_speed_um_per_s": left["nominal_speed_um_per_s"],
            "after_speed_um_per_s": right["nominal_speed_um_per_s"],
            "map_median_InvOLS_change_percent": 100 * (right["map_contact_invOLS_median_nm_per_V"] / left["map_contact_invOLS_median_nm_per_V"] - 1),
            "paired_InvOLS_median_change_nm_per_V": float(np.median(delta_s)),
            "InvOLS_positive_pixel_fraction": float(np.mean(delta_s > 0)),
            "detected_snap_distance_before_nm": detected_snap[before - 1],
            "detected_snap_distance_after_nm": detected_snap[after - 1],
            "F20_line_paired_median_change_pN": delta_f["median_paired_change_pN"],
            "F20_constant_paired_median_change_pN": delta_constant["median_paired_change_pN"],
            "F20_positive_pixel_fraction": delta_f["positive_pixel_fraction"],
            "F20_paired_pixel_count": delta_f["paired_pixels"],
            "blind_inference": "primary_refresh_candidate" if before == 18 else "secondary_liquid_or_optical_change_candidate" if before == 12 else "not_selected_as_refresh_candidate",
        })
    a.pilot.write_csv(a.OUT / "refresh_candidate_evidence.csv", candidate_rows)
    # These are descriptive paired comparisons at matching commanded speeds.
    matching = []
    for before, after in ((16, 20), (17, 19), (18, 21), (11, 13), (10, 14), (12, 15)):
        if maps[before - 1]["nominal_speed_um_per_s"] != maps[after - 1]["nominal_speed_um_per_s"]:
            raise RuntimeError("Comparison is not speed matched")
        matching.append({"before_map": before, "after_map": after, "speed_um_per_s": maps[before - 1]["nominal_speed_um_per_s"],
                         "before_detected_snap_distance_median_nm": detected_snap[before - 1], "after_detected_snap_distance_median_nm": detected_snap[after - 1],
                         "before_map_contact_InvOLS_nm_per_V": maps[before - 1]["map_contact_invOLS_median_nm_per_V"],
                         "after_map_contact_InvOLS_nm_per_V": maps[after - 1]["map_contact_invOLS_median_nm_per_V"]})
    a.pilot.write_csv(a.OUT / "refresh_same_speed_comparisons.csv", matching)
    fig, axes = plt.subplots(2, 2, figsize=(15, 9), sharex=True)
    plot_values = [np.nanmedian(curves[:, :, a.TARGET_INDICES[0]], axis=1) / 1000,
                   [r["map_contact_invOLS_median_nm_per_V"] for r in maps], detected_snap,
                   [r["idle_before_map_s"] for r in maps]]
    labels = ["20 nm apparent force: map median (nN)", "Hard-contact InvOLS: map median (nm/V)",
              "Detected approach snap-in distance (nm)", "Idle time before map (s)"]
    for panel, (ax, values, label) in enumerate(zip(axes.flat, plot_values, labels, strict=True)):
        time_basis = "map_start_time" if panel == 3 else "map_midpoint_time"
        times = a.measurement_times(maps, time_basis)
        ax.plot(times, values, color="0.45", lw=1)
        ax.scatter(times, values, c=[a.COLORS[r["nominal_speed_um_per_s"]] for r in maps], s=32, zorder=3)
        for before, color, style, legend in ((18, "#9b2226", "-", "Primary interval: #18 → #19"), (12, "#a97613", ":", "Secondary interval: #12 → #13")):
            begin, end = a.interval_between_maps(maps, before)
            ax.axvspan(begin, end, color=color, alpha=0.07)
            ax.axvline((begin + end) / 2, color=color, lw=1.8, ls=style, label=legend)
        ax.set_ylabel(label)
        a.axes_chronology(ax, maps, time_basis, show_blocks=False)
    axes[0, 1].legend(loc="best", frameon=False, fontsize=9)
    fig.suptitle("Blind liquid-refresh inference: measured changes and candidate intervals", fontsize=15)
    fig.text(0.5, 0.015, "20 nm values before map 19 often include snap-in trajectories. Candidate labels are hypotheses; identical scanner coordinates do not verify sample stability.", ha="center", fontsize=8.5)
    fig.tight_layout(rect=(0, 0.035, 1, 0.975))
    fig.savefig(a.FIG / "liquid_refresh_candidate_diagnostics.png", dpi=200)
    plt.close(fig)

    # Inspect the actual trajectory underlying the sparse near-contact bins.
    selected = (12, 13, 14, 18, 19, 25)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True, sharey=True)
    sensitivity = provenance["water_InvOLS_nm_per_V"] * 1e-9
    for ax, number in zip(axes.flat, selected, strict=True):
        source = a.base.load_source((a.ROOT / maps[number - 1]["source"]).resolve(), 0)
        curve = a.reference.curve_at_physical_pixel(source)
        far = a.base.fit_far_field_drift(curve.measured_height_m, curve.deflection_V)
        corrected = curve.deflection_V - a.base.baseline_voltage(curve.measured_height_m, curve.deflection_V, far)
        contact = a.pilot.terminal_contact_fit(curve.measured_height_m, corrected, 50)
        height0 = -contact["intercept_V"] / contact["slope_V_per_m"]
        displacement = sensitivity * corrected
        distance = (curve.measured_height_m + displacement - height0) * 1e9
        force = a.K * displacement * 1e9
        mask = np.arange(len(distance)) < contact["start"]
        ax.plot(distance[mask], force[mask], color="0.5", lw=0.75, marker=".", ms=2)
        rows = [r for r in fixed if r["acquisition_order"] == number]
        ax.scatter([r["distance_nm"] for r in rows], [r["force_linear_drift_corrected_pN"] / 1000 for r in rows],
                   c=a.COLORS[maps[number - 1]["nominal_speed_um_per_s"]], marker="s", s=33, zorder=3)
        ax.axhline(0, color="0.7", lw=0.6)
        for d in (20, 50):
            ax.axvline(d, color="0.8", lw=0.5, ls=":")
        ax.set_title(f"Map {number} · {maps[number - 1]['nominal_speed_um_per_s']:g} µm/s")
        ax.set_xlim(-5, 220)
        ax.set_ylim(-8, 3.5)
        ax.set_xlabel("Reconstructed separation D (nm)")
        ax.set_ylabel("Apparent force (nN)")
        ax.grid(alpha=0.15)
    fig.suptitle("Fixed pixel (3,3): raw sampled trajectories and 5 nm force slices", fontsize=15)
    fig.text(0.5, 0.015, "Gray: successive raw samples before the terminal fitting window. Squares: requested bins. Lines across jumps do not imply intermediate measurements.", ha="center", fontsize=8.5)
    fig.tight_layout(rect=(0, 0.04, 1, 0.965))
    fig.savefig(a.FIG / "fixed_pixel_raw_FD_diagnostics.png", dpi=180)
    plt.close(fig)

    report = f"""# 08-09-26 水中 AFM：固定像素图与换液位置盲测

已下载共享目录全部 **36 个 JPK force map，共 100,862,731 bytes**；逐文件大小、SHA-256 和 ZIP CRC 已核对。共处理 **2,304 条 approach 及其 retract**，两支 parser 均无跳过曲线。原始数据位于 `raw/keeper_6f30d59535114f89b568/08-09-26/`。

## 图与处理口径

- [主图：固定像素的四个距离切片](figures/fixed_pixel_row3_col3_force_slices.png)，另有 [SVG](figures/fixed_pixel_row3_col3_force_slices.svg) 和 [逐点 CSV](fixed_pixel_row3_col3_force_slices.csv)。横轴已改为目标像素实际 approach 段开始时间，显示仪器时区 UTC+02:00；该段持续0.25–1 s。真实停顿会直接表现为横向间隔，颜色仍表示速度。整map中位数与map级诊断采用整张map实际采集区间的中点，停顿诊断的点则对应下一张map的开始时间。
- 主图的 **median 是单条固定像素曲线内，每个5 nm距离箱的原始采样中位数**，没有跨64个像素平均。另存的 `map_median_force_slices.png` 才是先逐像素取距离箱median、再跨该map中有数据的像素取median。
- 沿用参考图的零基 `(row,column)=(3,3)`，即第 4 行第 4 列；8×8 蛇形扫描对应 archive index 28。全部 map 的扫描区域为 20×20 µm，记录的目标像素坐标相同：x={maps[0]['fixed_pixel_x_um']:.6f} µm、y={maps[0]['fixed_pixel_y_um']:.6f} µm。此检查不能测出样品自身漂移。
- 采用用户指定的 **cantilever 2，k=0.2736±0.0063 N/m**。文件头的 0.2384 N/m 是错误的 cantilever 3 参数，未用于力换算。单独修正 k 使相同重建间距处的力乘以 **{a.K / 0.2384:.9f}**；给定 k 不确定度对应共同的 **{100 * a.K_UNCERTAINTY / a.K:.3f}%** 力标度不确定度，其统计置信含义未指定。
- InvOLS 用本批水中原始硬接触段重算：50 nm 末端窗口、R²≥0.995、35–90 nm/V 初筛、每张 map 至少 56 条合格接触，再作 map 内 MAD 筛选；全局保留 **{provenance['water_InvOLS_retained_contacts']}/{provenance['total_contacts']}** 条，得到 **{provenance['water_InvOLS_nm_per_V']:.6f} nm/V**。40/60 nm 窗口对同一保留样本群的全局结果为 {provenance['contact_sensitivity_checks_nm_per_V']['40']:.6f}/{provenance['contact_sensitivity_checks_nm_per_V']['60']:.6f} nm/V。
- 逐曲线取起始 20% 原始采样拟合远场直线并扣除；以校正后的硬接触直线确定零点。`D=h+InvOLS×V_corrected−h_contact`，`F=k×InvOLS×V_corrected`。主图在 20、50、100、200 nm 的 `[D−2.5,D+2.5)` 区间取原始采样中位数，无插值或外推；另保留远场常数基线结果。
- 六组速度序列是 `[2,1,4,4,1,2]`、`[1,4,2,2,4,1]`、`[4,2,1,1,2,4]`，随后重复一次；单位 µm/s。主图竖虚线仍只表示速度回文分组，位于相邻分组之间实际停顿区间的中点；图底部的×表示无采样的距离箱。折线断开来自NaN空箱，竖虚线不是数据断裂或换液证明。

## 我的盲测判断

**首选：第 18→19 张之间换液。第二候选：第 12→13 张之间有一次补液、液层或光路状态变化；单凭这些数据，我对第二处是否确为换液的把握较低。** 如果只允许报一个位置，我选 18→19。

| 候选区间 | 无采集间隔，仪器时间 UTC+02:00 | 直接证据 | 判断边界 |
|---|---|---|---|
| **18→19** | **16:18:37.506–16:23:55.043，317.537 s** | map 接触 InvOLS 84.757→81.482 nm/V，降低3.86%；64/64 个配对像素均降低。20 nm 配对力变化中位数 +4.014 nN，63/64 个像素增大；常数基线下仍为 +3.960 nN。检测通过的 snap-in 距离中位数 **{detected_snap[17]:.2f}→{detected_snap[18]:.2f} nm**，后续保持短距离状态。 | 多项变化集中且持续，是最强的液体环境改变候选；数据本身不记录具体操作。近接触力跳变包含 snap-in 分支改变，不能当作平衡力增加4 nN。 |
| **12→13** | **15:45:44.700–15:49:51.191，246.491 s** | 接触 InvOLS 81.688→85.543 nm/V，增加4.72%；64/64 个配对像素均增加，13–18 张持续处在较高水平。同速度的11→13、10→14、12→15对照也保留此灵敏度变化。 | 接触/光学响应台阶可靠；snap-in 状态改变远弱于18→19，速度也同时改变，因此换液、补液或光路变化不能唯一分辨。 |

18→19 的结论还有同速度对照：1 µm/s 的16→20、2 µm/s 的17→19、4 µm/s 的18→21，检测通过的 snap-in 距离分别约 **{detected_snap[15]:.1f}→{detected_snap[19]:.1f}**、**{detected_snap[16]:.1f}→{detected_snap[18]:.1f}**、**{detected_snap[17]:.1f}→{detected_snap[20]:.1f} nm**。因此，这一变化不能仅用相邻图速度从4改成2 µm/s解释。相邻图的末端载荷仍约30.6 nN，retract pull-off没有同样的即时台阶；这些是限制证据。

第25张固定像素出现的高尖峰不列为换液候选：该点实际采集于16:34:23.270，速度1 µm/s；它的远场拟合斜率为+110.413 pN/100 nm，而前后第24/26张同像素分别为−6.224/−13.142 pN/100 nm。将这一远场正斜率向近表面外推并扣除，会把近表面整段力抬高。在100 nm处，远场直线基线给出815.0 pN，常数基线为61.2 pN，差753.8 pN；同图全map直线基线中位数仅17.8 pN。该箱6个原始采样共同受到基线平移，取median无法消除这种系统偏移。50 nm处两种基线分别为1033.7/224.8 pN，200 nm处为738.6/94.6 pN。数据支持局部远场斜率及其外推造成的敏感性；远场原始信号为何倾斜仍未确定，不能直接断言是具体仪器故障。

这一远场趋势也不是一个孤立坏点：原始200点窗口的直线R²约0.856，残差MAD约26.81 pN；前后半窗口斜率均为正，约+98.4/+159.3 pN/100 nm。[第25张基线诊断图](figures/map25_fixed_pixel_baseline_diagnostic.png)直接展示了外推偏移，[QC数值](map25_fixed_pixel_baseline_QC.csv)保留24/25/26张的对照。

第8–12张前的停顿约127–176 s，第13–18张前约242–265 s；这是明显的等待安排变化，但停顿本身不足以标记一次换液。第19张之后的逐图停顿恢复至约9–40 s。

## 必须保留的数值限制

- 144个固定像素目标箱中，139个有采样。20 nm处第1、3、5张为空；50 nm处第15、17张为空。jump-in附近相邻采样可以跨过整个5 nm距离箱；这些曲线和map都存在，只是目标箱内无点。主图保留空缺，CSV同时保留每个箱的原始采样数；没有用插值把缺口接起来。
- 例如第15张在1 ms内，重建距离由56.071跳到43.576 nm，跳过50 nm箱的47.5–52.5 nm；第1张在0.5 ms内由24.725跳到15.340 nm，跳过20 nm箱的17.5–22.5 nm。扫描器本身只移动约1 nm，较大的距离变化来自悬臂偏转突变。全部五处的相邻原始采样位置与时间步长见 [空箱跳跃明细](fixed_pixel_missing_bin_jumps.csv)。
- 第4、10、11、12张的20 nm箱，以及第18张的50 nm箱仅有1个原始点。第18张50 nm的−3204 pN尤其不稳健：改用map局部InvOLS或60 nm接触窗口后，该箱为空。其尖峰不能作为单点精确力结论。
- 前18张常在约20–60 nm发生jump-in，20 nm切片大量落在失稳/跳入后的轨迹；第19张之后的20 nm位于较平稳的排斥支。主图保留与参考图相同的观测定义，但两段20 nm值不等同于可直接比较的稳定平衡力。
- 固定点28的接触拟合R²均超过0.9994，局部InvOLS范围79.448–87.017 nm/V。光学响应有时间变化；全局InvOLS用于保持参考图口径，同时提供每张map的局部InvOLS及40/60 nm窗口敏感性列。
- 100–200 nm处基线选择会显著改变数值。没有对64个相关像素作独立实验的显著性检验，也没有把上述候选变成换液的已验证标签、概率或零速平衡力结论。
- 早期部分retract曲线未满足自由端基线/量程检查，pull-off有删失；完整计数在 `map_inventory_QC.csv`。pull-off图是有效记录的条件中位数，不能据此断言未观测的黏附值。

## 核对文件与复现

- [换液候选诊断图](figures/liquid_refresh_candidate_diagnostics.png)、[原始采样轨迹](figures/fixed_pixel_raw_FD_diagnostics.png)、[全map中位数切片](figures/map_median_force_slices.png)、[基线敏感性图](figures/fixed_pixel_row3_col3_baseline_sensitivity.png)。
- [候选证据](refresh_candidate_evidence.csv)保留全部35个相邻区间；[同速度对照](refresh_same_speed_comparisons.csv)、[map级QC](map_inventory_QC.csv)、[逐像素QC](pixel_QC.csv)、[相邻map差值](adjacent_map_changes.csv)提供可复核数值。
- [provenance.json](provenance.json)记录用户标定、下载来源、原始输入hash、算法、软件版本与139/144覆盖情况；原始重建与参考绘图函数在全部36个固定像素上的最大绝对差为0 pN。这是实现一致性检查，不是经验标定或科学解释的独立验证。

从仓库根目录执行：

```powershell
python -X utf8 analysis/download_08_09_26_keeper.py
python -X utf8 analysis/analyze_08_09_26_fixed_pixel.py
python -X utf8 analysis/assess_08_09_26_refresh.py
```

若完整像素重建已经完成，重试绘图/导出可使用 `analyze_08_09_26_fixed_pixel.py --reuse-reconstruction`；该模式重新核验原始输入hash和标定，并从原始数据重算全部36个固定像素以核对已保存矩阵。
"""
    (a.OUT / "REPORT.md").write_text(report, encoding="utf-8")
    (a.OUT / "refresh_inference_provenance.json").write_text(json.dumps({
        "claim_status": "blind_hypothesis_before_user_discloses_operations",
        "primary_candidate_between_maps": [18, 19], "secondary_candidate_between_maps": [12, 13],
        "secondary_scope": "liquid level or optical/contact-response change; actual liquid replacement not uniquely identified",
        "snap_summary_filter": "approach_snap_detected is True and distance finite; avoids unconfirmed detector candidates",
        "script_sha256": a.pilot.sha256_file(Path(__file__).resolve()),
        "input_hashes": {name: a.pilot.sha256_file(a.OUT / name) for name in ("map_inventory_QC.csv", "pixel_QC.csv", "fixed_pixel_row3_col3_force_slices.csv", "adjacent_map_changes.csv", "all_pixel_force_curves.npz", "provenance.json")},
    }, indent=2) + "\n", encoding="utf-8")
    artifacts = [p for p in a.OUT.rglob("*") if p.is_file() and p.name != "artifact_manifest.sha256"]
    artifacts += [Path(__file__).resolve(), Path(a.__file__).resolve(), a.ROOT / "analysis" / "download_08_09_26_keeper.py"]
    a.pilot.create_manifest(artifacts, a.OUT / "artifact_manifest.sha256")
    print(f"Wrote {a.OUT / 'REPORT.md'}", flush=True)


if __name__ == "__main__":
    main()

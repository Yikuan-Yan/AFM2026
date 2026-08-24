#!/usr/bin/env python3
"""Classical distribution-separation tests for the three 0 wt% AFM maps.

Tests are evaluated at every 5 nm bin from 25 to 250 nm.  The maps cover the
same physical raster, so paired tests are appropriate.  Two inferential levels
are reported:

1. pixel-paired tests (n=255), which are conventional but pseudoreplicate
   spatially correlated pixels;
2. primary 4x4-pixel block-paired tests (n=16), which first reduce every map to
   matched spatial-block medians.

Welch, Mann-Whitney, and two-sample KS tests are also included because they ask
about the marginal distributions themselves, but their independent-sample
assumption is not literally satisfied by matched, spatially correlated maps.
Every p-value family is Holm-adjusted across three pairs x 46 distances; BH-FDR
q-values are provided as a secondary multiplicity summary.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import platform

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats

import analyze_velocity_joint_fit as joint
import analyze_zero_wt_distribution_separation as separation


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "analysis" / "velocity_joint_fit_results"
FIGURES_DIR = RESULTS_DIR / "figures"

DISTANCE_MIN_NM = 25.0
DISTANCE_MAX_NM = 250.0
SELECTED_DISTANCES_NM = (25.0, 50.0, 100.0, 150.0, 200.0, 250.0)
SPEEDS_UM_PER_S = (1, 2, 4)
COMPARISONS = ((1, 2), (2, 4), (1, 4))
BLOCK_SIZE_PIXELS = 4
ALPHA = 0.05

DETAIL_CSV = RESULTS_DIR / "zero_wt_classical_distribution_tests.csv"
SELECTED_CSV = RESULTS_DIR / "zero_wt_classical_distribution_tests_selected.csv"
OMNIBUS_CSV = RESULTS_DIR / "zero_wt_classical_omnibus_tests.csv"
SUMMARY_CSV = RESULTS_DIR / "zero_wt_classical_test_summary.csv"
PIXEL_FORCE_CSV = RESULTS_DIR / "zero_wt_map_pixel_forces_25_250nm.csv"
REPORT_MD = RESULTS_DIR / "ZERO_WT_CLASSICAL_DISTRIBUTION_TESTS.md"
PROVENANCE_JSON = RESULTS_DIR / "zero_wt_classical_distribution_tests_provenance.json"
FIGURE_STEM = FIGURES_DIR / "zero_wt_classical_distribution_test_pvalues"


PAIRWISE_P_FIELDS = (
    "pixel_paired_t_p",
    "pixel_wilcoxon_p",
    "block_paired_t_p",
    "block_wilcoxon_p",
    "welch_t_p",
    "mann_whitney_p",
    "ks_2sample_p",
)
OMNIBUS_P_FIELDS = (
    "pixel_repeated_measures_anova_p",
    "pixel_friedman_p",
    "block_repeated_measures_anova_p",
    "block_friedman_p",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def holm_adjust(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    if p.ndim != 1 or np.any(~np.isfinite(p)) or np.any((p < 0.0) | (p > 1.0)):
        raise ValueError("Holm correction requires finite one-dimensional p-values")
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 0.0
    count = p.size
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (count - rank) * p[index]))
        adjusted[index] = running
    return adjusted


def benjamini_hochberg(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    if p.ndim != 1 or np.any(~np.isfinite(p)) or np.any((p < 0.0) | (p > 1.0)):
        raise ValueError("BH correction requires finite one-dimensional p-values")
    order = np.argsort(p)
    sorted_p = p[order]
    count = p.size
    adjusted_sorted = sorted_p * count / np.arange(1, count + 1)
    adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
    adjusted_sorted = np.clip(adjusted_sorted, 0.0, 1.0)
    adjusted = np.empty_like(p)
    adjusted[order] = adjusted_sorted
    return adjusted


def repeated_measures_anova(values: np.ndarray) -> tuple[float, float, int, int]:
    """One-factor repeated-measures ANOVA for complete subject x condition data."""

    data = np.asarray(values, dtype=float)
    if data.ndim != 2 or data.shape[0] < 3 or data.shape[1] < 2:
        raise ValueError("repeated-measures ANOVA needs subjects x conditions")
    if np.any(~np.isfinite(data)):
        raise ValueError("repeated-measures ANOVA does not silently omit values")
    subject_count, condition_count = data.shape
    grand = float(np.mean(data))
    subject_means = np.mean(data, axis=1)
    condition_means = np.mean(data, axis=0)
    total_ss = float(np.sum((data - grand) ** 2))
    subject_ss = float(condition_count * np.sum((subject_means - grand) ** 2))
    condition_ss = float(subject_count * np.sum((condition_means - grand) ** 2))
    error_ss = total_ss - subject_ss - condition_ss
    scale = max(total_ss, np.finfo(float).tiny)
    if error_ss < -1e-12 * scale:
        raise FloatingPointError(f"negative repeated-measures error SS: {error_ss}")
    error_ss = max(error_ss, 0.0)
    condition_df = condition_count - 1
    error_df = (subject_count - 1) * condition_df
    condition_ms = condition_ss / condition_df
    error_ms = error_ss / error_df
    if error_ms <= np.finfo(float).tiny:
        statistic = math.inf if condition_ms > 0.0 else 0.0
    else:
        statistic = condition_ms / error_ms
    p_value = float(stats.f.sf(statistic, condition_df, error_df))
    return float(statistic), p_value, condition_df, error_df


def t_confidence_interval(difference: np.ndarray) -> tuple[float, float, float, float]:
    values = np.asarray(difference, dtype=float)
    if values.ndim != 1 or values.size < 2 or np.any(~np.isfinite(values)):
        raise ValueError("paired t interval requires complete finite differences")
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=1))
    se = sd / math.sqrt(values.size)
    critical = float(stats.t.ppf(0.975, values.size - 1))
    low = mean - critical * se
    high = mean + critical * se
    effect_dz = mean / sd if sd > np.finfo(float).tiny else math.copysign(math.inf, mean)
    return mean, float(low), float(high), float(effect_dz)


def spatial_block_medians(
    forces: dict[int, np.ndarray],
    coordinates: list[tuple[int, int]],
) -> tuple[dict[int, np.ndarray], int]:
    _, members = separation.block_members(coordinates, BLOCK_SIZE_PIXELS)
    output = {
        speed: np.asarray(
            [np.median(forces[speed][indices], axis=0) for indices in members],
            dtype=float,
        )
        for speed in SPEEDS_UM_PER_S
    }
    expected_shape = (len(members), next(iter(forces.values())).shape[1])
    if any(values.shape != expected_shape for values in output.values()):
        raise AssertionError("block force shape mismatch")
    return output, len(members)


def add_multiplicity(rows: list[dict], fields: tuple[str, ...]) -> None:
    for field in fields:
        p = np.asarray([float(row[field]) for row in rows], dtype=float)
        holm = holm_adjust(p)
        bh = benjamini_hochberg(p)
        for index, row in enumerate(rows):
            row[f"{field}_holm"] = float(holm[index])
            row[f"{field}_bh_fdr"] = float(bh[index])


def significance_segments(
    records: list[dict],
    adjusted_field: str,
    sign_field: str | None = None,
) -> str:
    ordered = sorted(records, key=lambda row: float(row["distance_nm"]))
    status: list[int] = []
    for row in ordered:
        significant = float(row[adjusted_field]) < ALPHA
        if not significant:
            status.append(0)
        elif sign_field is None:
            status.append(1)
        else:
            value = float(row[sign_field])
            status.append(1 if value > 0.0 else -1 if value < 0.0 else 0)
    segments: list[tuple[float, float, int]] = []
    start: int | None = None
    current = 0
    for index, value in enumerate(status):
        if value != current:
            if current != 0 and start is not None:
                segments.append(
                    (
                        float(ordered[start]["distance_nm"]),
                        float(ordered[index - 1]["distance_nm"]),
                        current,
                    )
                )
            start = index if value != 0 else None
            current = value
    if current != 0 and start is not None:
        segments.append(
            (
                float(ordered[start]["distance_nm"]),
                float(ordered[-1]["distance_nm"]),
                current,
            )
        )
    if not segments:
        return "none"
    output = []
    for low, high, direction in segments:
        suffix = ""
        if sign_field is not None:
            suffix = "+" if direction > 0 else "-"
        output.append(f"{low:.0f}-{high:.0f} nm{suffix}")
    return "; ".join(output)


def safe_minus_log10(values: np.ndarray) -> np.ndarray:
    minimum = np.finfo(float).tiny
    return -np.log10(np.maximum(np.asarray(values, dtype=float), minimum))


def plot_pvalues(rows: list[dict]) -> None:
    colors = {
        "pixel_paired_t_p_holm": "#9ecae1",
        "block_paired_t_p_holm": "#08519c",
        "block_wilcoxon_p_holm": "#31a354",
        "ks_2sample_p_holm": "#e6550d",
    }
    labels = {
        "pixel_paired_t_p_holm": "pixel paired t (n=255; naive)",
        "block_paired_t_p_holm": "4x4-block paired t (n=16)",
        "block_wilcoxon_p_holm": "4x4-block Wilcoxon",
        "ks_2sample_p_holm": "two-sample KS (naive)",
    }
    figure, axes = plt.subplots(3, 1, figsize=(10.5, 9.0), sharex=True)
    for axis, (slow, fast) in zip(axes, COMPARISONS):
        selected = [
            row
            for row in rows
            if int(row["slow_speed_um_per_s"]) == slow
            and int(row["fast_speed_um_per_s"]) == fast
        ]
        selected.sort(key=lambda row: float(row["distance_nm"]))
        distance = np.asarray([float(row["distance_nm"]) for row in selected])
        for field, color in colors.items():
            adjusted_p = np.asarray([float(row[field]) for row in selected])
            axis.plot(
                distance,
                safe_minus_log10(adjusted_p),
                color=color,
                linewidth=2.0 if field == "block_paired_t_p_holm" else 1.35,
                linestyle="--" if "wilcoxon" in field else "-",
                label=labels[field],
            )
        axis.axhline(
            -math.log10(ALPHA),
            color="black",
            linewidth=0.9,
            linestyle=":",
            label="Holm p=0.05" if axis is axes[0] else None,
        )
        axis.set_ylabel(f"{fast} vs {slow}\n-log10(Holm p)")
        axis.grid(alpha=0.22)
    axes[0].legend(fontsize=8, ncol=2, loc="upper right")
    axes[-1].set_xlabel("Sphere-plane separation, D (nm)")
    figure.suptitle(
        "0 wt% classical distribution-separation tests",
        fontsize=14,
        y=0.98,
    )
    figure.text(
        0.5,
        0.012,
        (
            "Holm correction is over 3 pairwise comparisons x 46 distance bins "
            "within each test family. Pixel and KS tests ignore spatial dependence."
        ),
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="0.30",
    )
    figure.subplots_adjust(left=0.11, right=0.985, bottom=0.075, top=0.93, hspace=0.13)
    figure.savefig(FIGURE_STEM.with_suffix(".png"), dpi=300)
    figure.savefig(FIGURE_STEM.with_suffix(".svg"))
    plt.close(figure)


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/5] reconstruct 0 wt% line-corrected map force distributions", flush=True)
    sources, maps = separation.load_zero_wt_maps()
    coordinates = sorted(
        set.intersection(*(set(maps[speed].points) for speed in SPEEDS_UM_PER_S))
    )
    distance_mask = (
        (joint.BIN_CENTERS_NM >= DISTANCE_MIN_NM)
        & (joint.BIN_CENTERS_NM <= DISTANCE_MAX_NM)
    )
    distance_nm = joint.BIN_CENTERS_NM[distance_mask]
    forces = {
        speed: np.asarray(
            [
                maps[speed].points[coordinate].force_linear_corrected_pN[
                    distance_mask
                ]
                for coordinate in coordinates
            ],
            dtype=float,
        )
        for speed in SPEEDS_UM_PER_S
    }
    complete = np.ones(len(coordinates), dtype=bool)
    for speed in SPEEDS_UM_PER_S:
        complete &= np.all(np.isfinite(forces[speed]), axis=1)
    coordinates = [
        coordinate for coordinate, keep in zip(coordinates, complete) if keep
    ]
    forces = {speed: values[complete] for speed, values in forces.items()}
    pixel_count = len(coordinates)
    if pixel_count != 255:
        raise RuntimeError(f"expected 255 common complete pixels, got {pixel_count}")
    expected_shape = (pixel_count, distance_nm.size)
    if any(values.shape != expected_shape for values in forces.values()):
        raise AssertionError("pixel force shape mismatch")
    block_forces, block_count = spatial_block_medians(forces, coordinates)
    if block_count != 16:
        raise RuntimeError(f"expected 16 spatial blocks, got {block_count}")

    pixel_force_rows: list[dict] = []
    for pixel_index, (row_index, column_index) in enumerate(coordinates):
        for speed in SPEEDS_UM_PER_S:
            for distance_index, distance in enumerate(distance_nm):
                pixel_force_rows.append(
                    {
                        "concentration_wt_percent": 0,
                        "speed_um_per_s": speed,
                        "map_timestamp": maps[speed].timestamp,
                        "map_row": row_index,
                        "map_column": column_index,
                        "distance_nm": float(distance),
                        "force_linear_drift_corrected_pN": float(
                            forces[speed][pixel_index, distance_index]
                        ),
                    }
                )

    print("[2/5] pairwise classical tests", flush=True)
    rows: list[dict] = []
    for slow, fast in COMPARISONS:
        for distance_index, distance in enumerate(distance_nm):
            slow_pixel = forces[slow][:, distance_index]
            fast_pixel = forces[fast][:, distance_index]
            pixel_difference = fast_pixel - slow_pixel
            slow_block = block_forces[slow][:, distance_index]
            fast_block = block_forces[fast][:, distance_index]
            block_difference = fast_block - slow_block

            pixel_t = stats.ttest_rel(fast_pixel, slow_pixel)
            block_t = stats.ttest_rel(fast_block, slow_block)
            welch = stats.ttest_ind(fast_pixel, slow_pixel, equal_var=False)
            pixel_wilcoxon = stats.wilcoxon(
                pixel_difference,
                zero_method="wilcox",
                correction=False,
                alternative="two-sided",
                method="approx",
            )
            block_wilcoxon = stats.wilcoxon(
                block_difference,
                zero_method="wilcox",
                correction=False,
                alternative="two-sided",
                method="auto",
            )
            mann_whitney = stats.mannwhitneyu(
                fast_pixel,
                slow_pixel,
                alternative="two-sided",
                method="asymptotic",
            )
            ks = stats.ks_2samp(
                fast_pixel,
                slow_pixel,
                alternative="two-sided",
                method="auto",
            )
            pixel_mean, pixel_low, pixel_high, pixel_dz = t_confidence_interval(
                pixel_difference
            )
            block_mean, block_low, block_high, block_dz = t_confidence_interval(
                block_difference
            )
            shapiro = stats.shapiro(block_difference)

            values_to_check = np.asarray(
                [
                    pixel_t.statistic,
                    pixel_t.pvalue,
                    block_t.statistic,
                    block_t.pvalue,
                    welch.statistic,
                    welch.pvalue,
                    pixel_wilcoxon.statistic,
                    pixel_wilcoxon.pvalue,
                    block_wilcoxon.statistic,
                    block_wilcoxon.pvalue,
                    mann_whitney.statistic,
                    mann_whitney.pvalue,
                    ks.statistic,
                    ks.pvalue,
                    shapiro.statistic,
                    shapiro.pvalue,
                ],
                dtype=float,
            )
            if np.any(~np.isfinite(values_to_check)):
                raise FloatingPointError(
                    f"non-finite test result for {fast}-{slow} at {distance} nm"
                )
            rows.append(
                {
                    "concentration_wt_percent": 0,
                    "slow_speed_um_per_s": slow,
                    "fast_speed_um_per_s": fast,
                    "comparison": f"F_{fast}_minus_F_{slow}",
                    "distance_nm": float(distance),
                    "paired_pixel_count": pixel_count,
                    "spatial_block_size_pixels": BLOCK_SIZE_PIXELS,
                    "paired_block_count": block_count,
                    "pixel_mean_difference_pN": pixel_mean,
                    "pixel_mean_difference_t_ci95_low_pN": pixel_low,
                    "pixel_mean_difference_t_ci95_high_pN": pixel_high,
                    "pixel_paired_cohen_dz": pixel_dz,
                    "pixel_paired_t_statistic": float(pixel_t.statistic),
                    "pixel_paired_t_df": pixel_count - 1,
                    "pixel_paired_t_p": float(pixel_t.pvalue),
                    "pixel_wilcoxon_statistic": float(pixel_wilcoxon.statistic),
                    "pixel_wilcoxon_p": float(pixel_wilcoxon.pvalue),
                    "block_mean_difference_pN": block_mean,
                    "block_mean_difference_t_ci95_low_pN": block_low,
                    "block_mean_difference_t_ci95_high_pN": block_high,
                    "block_paired_cohen_dz": block_dz,
                    "block_paired_t_statistic": float(block_t.statistic),
                    "block_paired_t_df": block_count - 1,
                    "block_paired_t_p": float(block_t.pvalue),
                    "block_wilcoxon_statistic": float(block_wilcoxon.statistic),
                    "block_wilcoxon_p": float(block_wilcoxon.pvalue),
                    "block_difference_shapiro_W": float(shapiro.statistic),
                    "block_difference_shapiro_p": float(shapiro.pvalue),
                    "welch_t_statistic": float(welch.statistic),
                    "welch_t_p": float(welch.pvalue),
                    "mann_whitney_U": float(mann_whitney.statistic),
                    "mann_whitney_p": float(mann_whitney.pvalue),
                    "ks_2sample_statistic": float(ks.statistic),
                    "ks_2sample_p": float(ks.pvalue),
                    "inference_scope": (
                        "block tests: conditional mapped-area separation; pixel/"
                        "marginal tests ignore spatial dependence; one map per speed"
                    ),
                }
            )
    add_multiplicity(rows, PAIRWISE_P_FIELDS)

    print("[3/5] three-distribution omnibus tests", flush=True)
    omnibus_rows: list[dict] = []
    for distance_index, distance in enumerate(distance_nm):
        pixel_matrix = np.column_stack(
            [forces[speed][:, distance_index] for speed in SPEEDS_UM_PER_S]
        )
        block_matrix = np.column_stack(
            [block_forces[speed][:, distance_index] for speed in SPEEDS_UM_PER_S]
        )
        pixel_anova = repeated_measures_anova(pixel_matrix)
        block_anova = repeated_measures_anova(block_matrix)
        pixel_friedman = stats.friedmanchisquare(
            *(pixel_matrix[:, index] for index in range(pixel_matrix.shape[1]))
        )
        block_friedman = stats.friedmanchisquare(
            *(block_matrix[:, index] for index in range(block_matrix.shape[1]))
        )
        omnibus_rows.append(
            {
                "concentration_wt_percent": 0,
                "distance_nm": float(distance),
                "speed_groups_um_per_s": "1;2;4",
                "paired_pixel_count": pixel_count,
                "paired_block_count": block_count,
                "pixel_repeated_measures_anova_F": pixel_anova[0],
                "pixel_repeated_measures_anova_df_condition": pixel_anova[2],
                "pixel_repeated_measures_anova_df_error": pixel_anova[3],
                "pixel_repeated_measures_anova_p": pixel_anova[1],
                "pixel_friedman_chi2": float(pixel_friedman.statistic),
                "pixel_friedman_p": float(pixel_friedman.pvalue),
                "block_repeated_measures_anova_F": block_anova[0],
                "block_repeated_measures_anova_df_condition": block_anova[2],
                "block_repeated_measures_anova_df_error": block_anova[3],
                "block_repeated_measures_anova_p": block_anova[1],
                "block_friedman_chi2": float(block_friedman.statistic),
                "block_friedman_p": float(block_friedman.pvalue),
                "inference_scope": (
                    "block omnibus tests are primary conditional mapped-area tests; "
                    "not replicated map-level velocity inference"
                ),
            }
        )
    add_multiplicity(omnibus_rows, OMNIBUS_P_FIELDS)

    selected_rows = [
        row for row in rows if float(row["distance_nm"]) in SELECTED_DISTANCES_NM
    ]
    summary_rows: list[dict] = []
    for slow, fast in COMPARISONS:
        selected_pair = [
            row
            for row in rows
            if int(row["slow_speed_um_per_s"]) == slow
            and int(row["fast_speed_um_per_s"]) == fast
        ]
        summary_rows.append(
            {
                "concentration_wt_percent": 0,
                "slow_speed_um_per_s": slow,
                "fast_speed_um_per_s": fast,
                "comparison": f"F_{fast}_minus_F_{slow}",
                "pixel_paired_t_holm_significant_segments": significance_segments(
                    selected_pair,
                    "pixel_paired_t_p_holm",
                    "pixel_mean_difference_pN",
                ),
                "block_paired_t_holm_significant_segments": significance_segments(
                    selected_pair,
                    "block_paired_t_p_holm",
                    "block_mean_difference_pN",
                ),
                "block_wilcoxon_holm_significant_segments": significance_segments(
                    selected_pair,
                    "block_wilcoxon_p_holm",
                    "block_mean_difference_pN",
                ),
                "ks_holm_significant_segments": significance_segments(
                    selected_pair,
                    "ks_2sample_p_holm",
                    None,
                ),
                "minimum_block_paired_t_holm_p": min(
                    float(row["block_paired_t_p_holm"]) for row in selected_pair
                ),
                "minimum_ks_holm_p": min(
                    float(row["ks_2sample_p_holm"]) for row in selected_pair
                ),
                "inference_scope": (
                    "classical separation across distance; Holm correction within "
                    "each test family"
                ),
            }
        )

    print("[4/5] outputs and figure", flush=True)
    write_csv(DETAIL_CSV, rows)
    write_csv(SELECTED_CSV, selected_rows)
    write_csv(OMNIBUS_CSV, omnibus_rows)
    write_csv(SUMMARY_CSV, summary_rows)
    write_csv(PIXEL_FORCE_CSV, pixel_force_rows)
    plot_pvalues(rows)

    report_lines = [
        "# 0 wt%三组F-D分布的经典假设检验",
        "",
        "## 检验对象",
        "",
        (
            "每个距离D上直接检验1、2、4 um/s三组line-corrected force分布。"
            f"同物理位置配对pixel数为{pixel_count}；主检验先取4x4空间block中位数，"
            f"得到{block_count}个配对block。"
        ),
        "",
        (
            "pairwise paired t-test检验均值差是否为0；Wilcoxon检验配对差的中心"
            "是否为0；Welch/Mann-Whitney/KS把两组边缘分布作为独立样本比较。"
            "所有pairwise检验分别在3 pairs x 46 distances内作Holm校正。"
        ),
        "",
        "## 显著距离区间",
        "",
        "| comparison | pixel paired t | 4x4-block paired t | block Wilcoxon | two-sample KS |",
        "|:---|:---|:---|:---|:---|",
    ]
    for row in summary_rows:
        report_lines.append(
            f"| {int(row['fast_speed_um_per_s'])}-{int(row['slow_speed_um_per_s'])} um/s "
            f"| {row['pixel_paired_t_holm_significant_segments']} "
            f"| {row['block_paired_t_holm_significant_segments']} "
            f"| {row['block_wilcoxon_holm_significant_segments']} "
            f"| {row['ks_holm_significant_segments']} |"
        )
    block_anova_segments = significance_segments(
        omnibus_rows,
        "block_repeated_measures_anova_p_holm",
        None,
    )
    block_friedman_segments = significance_segments(
        omnibus_rows,
        "block_friedman_p_holm",
        None,
    )
    report_lines.extend(
        [
            "",
            "符号+/-表示fast-minus-slow的均值方向；KS只判断CDF是否相同，不定义方向。",
            "",
            "## 三组总体检验",
            "",
            (
                "4x4-block repeated-measures ANOVA在Holm校正后显著的距离区间为 "
                f"`{block_anova_segments}`；block Friedman检验为 "
                f"`{block_friedman_segments}`。总体原假设是1、2、4 um/s三组"
                "位置参数相同。"
            ),
            "",
            "## 选定距离的4x4-block paired t-test",
            "",
            "| comparison | D (nm) | mean difference (pN) | 95% t-CI (pN) | Cohen dz | raw p | Holm p | block Wilcoxon Holm p | KS Holm p |",
            "|:---|---:|---:|:---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in selected_rows:
        report_lines.append(
            f"| {int(row['fast_speed_um_per_s'])}-{int(row['slow_speed_um_per_s'])} "
            f"| {float(row['distance_nm']):.0f} "
            f"| {float(row['block_mean_difference_pN']):.4g} "
            f"| [{float(row['block_mean_difference_t_ci95_low_pN']):.4g}, "
            f"{float(row['block_mean_difference_t_ci95_high_pN']):.4g}] "
            f"| {float(row['block_paired_cohen_dz']):.3g} "
            f"| {float(row['block_paired_t_p']):.3g} "
            f"| {float(row['block_paired_t_p_holm']):.3g} "
            f"| {float(row['block_wilcoxon_p_holm']):.3g} "
            f"| {float(row['ks_2sample_p_holm']):.3g} |"
        )
    report_lines.extend(
        [
            "",
            "## 推断限制",
            "",
            (
                "pixel-level paired t/Wilcoxon与Welch/Mann-Whitney/KS均把空间pixel"
                "视为独立或近似独立，因此p值通常过小。4x4-block结果是本表推荐的"
                "传统检验，但它仍只说明这三张已测map在该区域内分离。每个速度只有"
                "一张map，不能据此得到重复实验层面的速度因果t-test。"
            ),
            "",
        ]
    )
    REPORT_MD.write_text("\n".join(report_lines), encoding="utf-8")

    provenance = {
        "analysis_script": str(Path(__file__).relative_to(ROOT)),
        "analysis_script_sha256": sha256_file(Path(__file__)),
        "dependency_script_sha256": {
            "analyze_zero_wt_distribution_separation.py": sha256_file(
                ROOT / "analysis" / "analyze_zero_wt_distribution_separation.py"
            ),
            "analyze_velocity_joint_fit.py": sha256_file(
                ROOT / "analysis" / "analyze_velocity_joint_fit.py"
            ),
        },
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "force_unit": "pN",
        "distance_unit": "nm",
        "distance_window_nm": [DISTANCE_MIN_NM, DISTANCE_MAX_NM],
        "distance_step_nm": float(np.diff(distance_nm)[0]),
        "paired_pixel_count": pixel_count,
        "spatial_block_size_pixels": BLOCK_SIZE_PIXELS,
        "paired_block_count": block_count,
        "map_acquisition_order_um_per_s": [2, 1, 4],
        "multiple_comparison_scope": "3 pairs x 46 distances separately per test family",
        "primary_test": "paired t-test on matched 4x4 spatial-block medians",
        "force_branch": "per-curve far-linear-drift-corrected approach",
        "sensitivity_nm_per_V": maps[1].sensitivity_nm_per_V,
        "spring_constant_N_per_m": separation.base.SPRING_CONSTANT_N_PER_M,
        "raw_sources": [
            {
                "path": str(source.path.relative_to(ROOT)),
                "sha256": source.sha256,
                "source_type": source.source_type,
                "timestamp": source.timestamp.isoformat(),
            }
            for source in sources
        ],
        "inference_limit": (
            "one map per speed; classical tests quantify conditional mapped-area "
            "separation, not independent map-level velocity replication"
        ),
    }
    PROVENANCE_JSON.write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    print("[5/5] summary", flush=True)
    print(
        "three-group block rm-ANOVA: "
        f"{block_anova_segments}; block Friedman: {block_friedman_segments}"
    )
    for row in summary_rows:
        print(
            f"{int(row['fast_speed_um_per_s'])}-{int(row['slow_speed_um_per_s'])} um/s: "
            f"block paired-t={row['block_paired_t_holm_significant_segments']}; "
            f"block Wilcoxon={row['block_wilcoxon_holm_significant_segments']}; "
            f"KS={row['ks_holm_significant_segments']}"
        )
    print(f"wrote {DETAIL_CSV}")
    print(f"wrote {SELECTED_CSV}")
    print(f"wrote {OMNIBUS_CSV}")
    print(f"wrote {SUMMARY_CSV}")
    print(f"wrote {PIXEL_FORCE_CSV}")
    print(f"wrote {REPORT_MD}")
    print(f"wrote {FIGURE_STEM.with_suffix('.png')}")
    print(f"wrote {FIGURE_STEM.with_suffix('.svg')}")
    print(f"wrote {PROVENANCE_JSON}")


if __name__ == "__main__":
    main()

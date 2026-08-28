#!/usr/bin/env python3
"""Paired distribution tests after palindrome averaging of the 27-08-26 maps.

The analysis keeps three statistical levels separate:

1. 64 same-position pixel differences within one palindrome block;
2. 16 matched 2x2-pixel tile means, reducing local spatial pseudoreplication;
3. three block-level spatial medians for Test A (blocks 1-3) or Test B
   (blocks 3-5), which are the actual experimental repetitions.

The primary transformation is F_sym=(F_early+F_late)/2.  A sensitivity branch
also estimates a common within-block linear drift from the early/late
differences and translates each speed pair to a common clock-time center.  It
does not claim to remove nonlinear relaxation.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import platform

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import scipy
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = ROOT / "analysis" / "palindrome_27_08_26_full_results"
RESULTS = ROOT / "analysis" / "palindrome_27_08_26_distribution_results"
FIGURES = RESULTS / "figures"

DISTANCE_MIN_NM = 20.0
DISTANCE_MAX_NM = 200.0
TARGET_DISTANCES_NM = (20.0, 50.0, 100.0, 200.0)
SPEEDS = (0.05, 0.1, 0.2)
COMPARISONS = ((0.05, 0.1), (0.1, 0.2), (0.05, 0.2))
BASELINES = ("linear_drift_corrected", "far_constant_referenced")
PREPROCESSING = ("palindrome_mean", "linear_center_aligned")
TILE_SIZE_PIXELS = 2
ALPHA = 0.05

COLORS = {0.05: "#2a9d8f", 0.1: "#e9b949", 0.2: "#e76f51"}
PREP_COLORS = {
    "palindrome_mean": "#264653",
    "linear_center_aligned": "#9b5de5",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def holm_adjust(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    if p.ndim != 1 or p.size == 0 or np.any(~np.isfinite(p)):
        raise ValueError("Holm adjustment requires a nonempty finite vector")
    if np.any((p < 0.0) | (p > 1.0)):
        raise ValueError("p-values outside [0,1]")
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 0.0
    count = p.size
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (count - rank) * p[index]))
        adjusted[index] = running
    return adjusted


def bh_adjust(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    if p.ndim != 1 or p.size == 0 or np.any(~np.isfinite(p)):
        raise ValueError("BH adjustment requires a nonempty finite vector")
    order = np.argsort(p)
    sorted_p = p[order]
    count = p.size
    adjusted_sorted = sorted_p * count / np.arange(1, count + 1)
    adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
    adjusted = np.empty_like(p)
    adjusted[order] = np.clip(adjusted_sorted, 0.0, 1.0)
    return adjusted


def add_multiplicity(
    frame: pd.DataFrame,
    group_columns: list[str],
    p_column: str,
    suffix: str,
) -> None:
    holm_column = f"{p_column}_holm_{suffix}"
    bh_column = f"{p_column}_bh_{suffix}"
    frame[holm_column] = np.nan
    frame[bh_column] = np.nan
    for _, indices in frame.groupby(group_columns, sort=False).groups.items():
        index = np.asarray(list(indices), dtype=int)
        finite = np.isfinite(frame.loc[index, p_column].to_numpy(dtype=float))
        selected = index[finite]
        if not selected.size:
            continue
        values = frame.loc[selected, p_column].to_numpy(dtype=float)
        frame.loc[selected, holm_column] = holm_adjust(values)
        frame.loc[selected, bh_column] = bh_adjust(values)


def paired_t_summary(difference: np.ndarray) -> dict[str, float | int]:
    values = np.asarray(difference, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2:
        return {
            "sample_count": int(values.size),
            "mean_difference_pN": float("nan"),
            "median_difference_pN": float("nan"),
            "difference_sd_pN": float("nan"),
            "mean_95CI_low_pN": float("nan"),
            "mean_95CI_high_pN": float("nan"),
            "paired_t_statistic": float("nan"),
            "paired_t_df": float("nan"),
            "paired_t_p_two_sided": float("nan"),
            "cohen_dz": float("nan"),
            "fraction_high_greater_than_low": float("nan"),
        }
    count = int(values.size)
    mean = float(np.mean(values))
    median = float(np.median(values))
    sd = float(np.std(values, ddof=1))
    se = sd / math.sqrt(count)
    critical = float(stats.t.ppf(0.975, count - 1))
    scale = max(float(np.max(np.abs(values))), 1.0)
    if sd <= 32.0 * np.finfo(float).eps * scale:
        if abs(mean) <= 32.0 * np.finfo(float).eps * scale:
            statistic = 0.0
            p_value = 1.0
            effect = 0.0
        else:
            statistic = math.copysign(math.inf, mean)
            p_value = 0.0
            effect = math.copysign(math.inf, mean)
    else:
        statistic = mean / se
        p_value = float(2.0 * stats.t.sf(abs(statistic), count - 1))
        effect = mean / sd
    return {
        "sample_count": count,
        "mean_difference_pN": mean,
        "median_difference_pN": median,
        "difference_sd_pN": sd,
        "mean_95CI_low_pN": mean - critical * se,
        "mean_95CI_high_pN": mean + critical * se,
        "paired_t_statistic": statistic,
        "paired_t_df": count - 1,
        "paired_t_p_two_sided": p_value,
        "cohen_dz": effect,
        "fraction_high_greater_than_low": float(
            np.mean(values > 0.0) + 0.5 * np.mean(values == 0.0)
        ),
    }


def exact_sign_test(difference: np.ndarray) -> tuple[int, float]:
    values = np.asarray(difference, dtype=float)
    values = values[np.isfinite(values) & (values != 0.0)]
    if not values.size:
        return 0, float("nan")
    positive = int(np.count_nonzero(values > 0.0))
    tail = min(positive, values.size - positive)
    p_value = min(1.0, 2.0 * float(stats.binom.cdf(tail, values.size, 0.5)))
    return int(values.size), p_value


def tile_difference_means(
    difference: np.ndarray,
    coordinates: np.ndarray,
) -> tuple[np.ndarray, int, int]:
    values = np.asarray(difference, dtype=float)
    if values.shape != (64,) or coordinates.shape != (64, 2):
        raise ValueError("tile input shape mismatch")
    groups: dict[tuple[int, int], list[int]] = {}
    for index, (row, column) in enumerate(coordinates):
        key = (int(row) // TILE_SIZE_PIXELS, int(column) // TILE_SIZE_PIXELS)
        groups.setdefault(key, []).append(index)
    if len(groups) != 16 or any(len(members) != 4 for members in groups.values()):
        raise AssertionError("expected sixteen 2x2 tiles")
    output: list[float] = []
    counts: list[int] = []
    for key in sorted(groups):
        selected = values[groups[key]]
        finite = selected[np.isfinite(selected)]
        if finite.size >= 3:
            output.append(float(np.mean(finite)))
            counts.append(int(finite.size))
    if len(output) < 12:
        raise RuntimeError("fewer than twelve usable spatial tiles")
    return np.asarray(output), min(counts), max(counts)


def p_text(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value < 1e-4:
        return f"{value:.1e}"
    return f"{value:.4f}"


def significant_segments(records: pd.DataFrame, p_column: str) -> str:
    ordered = records.sort_values("distance_nm")
    active = (
        ordered[p_column].to_numpy(dtype=float) < ALPHA
    )
    sign = np.sign(ordered.mean_difference_pN.to_numpy(dtype=float)).astype(int)
    distance = ordered.distance_nm.to_numpy(dtype=float)
    segments: list[str] = []
    start: int | None = None
    current = 0
    for index, (is_active, direction) in enumerate(zip(active, sign, strict=True)):
        state = int(direction) if is_active else 0
        if state != current:
            if current != 0 and start is not None:
                suffix = "+" if current > 0 else "−"
                segments.append(
                    f"{distance[start]:.0f}–{distance[index-1]:.0f} nm {suffix}"
                )
            start = index if state != 0 else None
            current = state
    if current != 0 and start is not None:
        suffix = "+" if current > 0 else "−"
        segments.append(f"{distance[start]:.0f}–{distance[-1]:.0f} nm {suffix}")
    return "; ".join(segments) if segments else "none"


def load_inputs() -> tuple[
    np.ndarray,
    np.ndarray,
    dict[str, np.ndarray],
    pd.DataFrame,
    pd.DataFrame,
    np.ndarray,
    dict,
]:
    map_frame = pd.read_csv(UPSTREAM / "map_inventory_QC.csv")
    pair_frame = pd.read_csv(UPSTREAM / "pair_inventory.csv")
    pixel_frame = pd.read_csv(UPSTREAM / "pixel_QC.csv")
    provenance = json.loads((UPSTREAM / "provenance.json").read_text())
    with np.load(UPSTREAM / "pixel_force_curves.npz") as archive:
        distance = archive["distance_nm"].astype(float)
        sources = archive["source"].astype(str)
        arrays = {
            "linear_drift_corrected": archive[
                "force_linear_drift_corrected_pN"
            ].astype(float),
            "far_constant_referenced": archive[
                "force_far_constant_referenced_pN"
            ].astype(float),
        }
    if len(map_frame) != 29 or sources.shape != (29,):
        raise AssertionError("upstream 29-map contract failed")
    if any(values.shape != (29, 64, 100) for values in arrays.values()):
        raise AssertionError("upstream pixel-force shape failed")
    selected = (distance >= DISTANCE_MIN_NM) & (distance <= DISTANCE_MAX_NM)
    selected_distance = distance[selected]
    if selected_distance.shape != (37,) or not np.allclose(np.diff(selected_distance), 5.0):
        raise AssertionError("expected 37 five-nanometre bins from 20 to 200 nm")
    arrays = {key: values[:, :, selected] for key, values in arrays.items()}

    coordinate_reference: np.ndarray | None = None
    for source in sources:
        rows = pixel_frame[pixel_frame.source == source].sort_values("point_index")
        if len(rows) != 64 or not np.array_equal(
            rows.point_index.to_numpy(dtype=int), np.arange(64)
        ):
            raise AssertionError(f"invalid point-index map: {source}")
        coordinates = rows[["row", "column"]].to_numpy(dtype=int)
        if coordinate_reference is None:
            coordinate_reference = coordinates
        elif not np.array_equal(coordinates, coordinate_reference):
            raise AssertionError("point-index to physical-pixel mapping changed")
    assert coordinate_reference is not None
    if len(np.unique(coordinate_reference, axis=0)) != 64:
        raise AssertionError("physical pixel coordinates are not unique")
    return (
        selected_distance,
        sources,
        arrays,
        map_frame,
        pair_frame,
        coordinate_reference,
        provenance,
    )


def build_palindrome_forces(
    distance: np.ndarray,
    sources: np.ndarray,
    arrays: dict[str, np.ndarray],
    map_frame: pd.DataFrame,
    pair_frame: pd.DataFrame,
) -> tuple[
    dict[tuple[int, float, str, str], np.ndarray],
    dict[tuple[int, float], dict],
    pd.DataFrame,
    pd.DataFrame,
]:
    source_index = {source: index for index, source in enumerate(sources)}
    map_lookup = map_frame.set_index("source")
    complete = pair_frame[pair_frame.pair_status == "complete_pair"].copy()
    if len(complete) != 14:
        raise AssertionError("expected fourteen complete palindrome pairs")

    pair_metadata: dict[tuple[int, float], dict] = {}
    raw_sym: dict[tuple[int, float, str], np.ndarray] = {}
    early_late: dict[tuple[int, float, str], tuple[np.ndarray, np.ndarray]] = {}
    timing_rows: list[dict] = []
    for row in complete.itertuples(index=False):
        block = int(row.block)
        speed = float(row.nominal_speed_um_per_s)
        early_source = str(row.early_source)
        late_source = str(row.late_source)
        early_map = map_lookup.loc[early_source]
        late_map = map_lookup.loc[late_source]
        early_time = float(early_map.map_protocol_midpoint_epoch_s)
        late_time = float(late_map.map_protocol_midpoint_epoch_s)
        if not late_time > early_time:
            raise AssertionError("palindrome late map does not follow early map")
        center = 0.5 * (early_time + late_time)
        half_span_min = 0.5 * (late_time - early_time) / 60.0
        actual_speed = float(
            np.median(
                [
                    early_map.gap_speed_20_200nm_median_um_per_s,
                    late_map.gap_speed_20_200nm_median_um_per_s,
                ]
            )
        )
        pair_metadata[(block, speed)] = {
            "block": block,
            "nominal_speed_um_per_s": speed,
            "actual_gap_speed_um_per_s": actual_speed,
            "early_source": early_source,
            "late_source": late_source,
            "early_midpoint_epoch_s": early_time,
            "late_midpoint_epoch_s": late_time,
            "pair_center_epoch_s": center,
            "pair_half_span_min": half_span_min,
        }
        for baseline in BASELINES:
            early = arrays[baseline][source_index[early_source]]
            late = arrays[baseline][source_index[late_source]]
            raw_sym[(block, speed, baseline)] = 0.5 * (early + late)
            early_late[(block, speed, baseline)] = (early, late)

    for block in range(1, 6):
        members = [value for key, value in pair_metadata.items() if key[0] == block]
        center_mean = float(np.mean([item["pair_center_epoch_s"] for item in members]))
        for item in sorted(members, key=lambda value: value["nominal_speed_um_per_s"]):
            timing_rows.append(
                {
                    **item,
                    "pair_center_offset_from_block_mean_min": (
                        item["pair_center_epoch_s"] - center_mean
                    )
                    / 60.0,
                    "block_pair_center_spread_min": (
                        max(value["pair_center_epoch_s"] for value in members)
                        - min(value["pair_center_epoch_s"] for value in members)
                    )
                    / 60.0,
                    "time_semantics": "JPK map start plus half stored approach+retract protocol duration; XY and baseline overhead excluded",
                }
            )

    force: dict[tuple[int, float, str, str], np.ndarray] = {}
    alignment_rows: list[dict] = []
    for baseline in BASELINES:
        for block in range(1, 6):
            speeds = sorted(speed for b, speed in pair_metadata if b == block)
            centers = np.asarray(
                [pair_metadata[(block, speed)]["pair_center_epoch_s"] for speed in speeds]
            )
            common_center = float(np.mean(centers))
            numerator = np.zeros((64, distance.size), dtype=float)
            denominator = np.zeros_like(numerator)
            contributing_pairs = np.zeros_like(numerator, dtype=int)
            for speed in speeds:
                early, late = early_late[(block, speed, baseline)]
                dt_min = (
                    pair_metadata[(block, speed)]["late_midpoint_epoch_s"]
                    - pair_metadata[(block, speed)]["early_midpoint_epoch_s"]
                ) / 60.0
                difference = late - early
                valid = np.isfinite(difference)
                numerator[valid] += dt_min * difference[valid]
                denominator[valid] += dt_min**2
                contributing_pairs[valid] += 1
            drift = np.full_like(numerator, np.nan)
            supported = (denominator > 0.0) & (contributing_pairs >= 2)
            drift[supported] = numerator[supported] / denominator[supported]
            for speed in speeds:
                sym = raw_sym[(block, speed, baseline)]
                offset_min = (
                    common_center
                    - pair_metadata[(block, speed)]["pair_center_epoch_s"]
                ) / 60.0
                aligned = sym + drift * offset_min
                force[(block, speed, baseline, "palindrome_mean")] = sym
                force[(block, speed, baseline, "linear_center_aligned")] = aligned
            for index, value in enumerate(distance):
                valid_drift = drift[:, index][np.isfinite(drift[:, index])]
                adjustments = []
                for speed in speeds:
                    sym = raw_sym[(block, speed, baseline)][:, index]
                    aligned = force[
                        (block, speed, baseline, "linear_center_aligned")
                    ][:, index]
                    delta = aligned - sym
                    adjustments.extend(delta[np.isfinite(delta)].tolist())
                alignment_rows.append(
                    {
                        "block": block,
                        "baseline_method": baseline,
                        "distance_nm": float(value),
                        "speed_pairs_used_for_drift": len(speeds),
                        "drift_supported_pixels": int(valid_drift.size),
                        "common_drift_median_pN_per_min": float(np.median(valid_drift)),
                        "common_drift_q25_pN_per_min": float(
                            np.quantile(valid_drift, 0.25)
                        ),
                        "common_drift_q75_pN_per_min": float(
                            np.quantile(valid_drift, 0.75)
                        ),
                        "absolute_center_adjustment_median_pN": float(
                            np.median(np.abs(adjustments))
                        ),
                        "absolute_center_adjustment_max_pN": float(
                            np.max(np.abs(adjustments))
                        ),
                        "claim_status": "linear common-drift sensitivity only; nonlinear relaxation remains",
                    }
                )
    return force, pair_metadata, pd.DataFrame(timing_rows), pd.DataFrame(alignment_rows)


def distribution_tests(
    distance: np.ndarray,
    force: dict[tuple[int, float, str, str], np.ndarray],
    pair_metadata: dict[tuple[int, float], dict],
    coordinates: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pixel_rows: list[dict] = []
    tile_rows: list[dict] = []
    for preprocessing in PREPROCESSING:
        for baseline in BASELINES:
            for block in range(1, 6):
                available = {
                    speed
                    for b, speed, method, prep in force
                    if b == block and method == baseline and prep == preprocessing
                }
                for low_speed, high_speed in COMPARISONS:
                    if low_speed not in available or high_speed not in available:
                        continue
                    low = force[(block, low_speed, baseline, preprocessing)]
                    high = force[(block, high_speed, baseline, preprocessing)]
                    for index, value in enumerate(distance):
                        low_slice = low[:, index]
                        high_slice = high[:, index]
                        complete = np.isfinite(low_slice) & np.isfinite(high_slice)
                        difference = high_slice[complete] - low_slice[complete]
                        summary = paired_t_summary(difference)
                        common = {
                            "block": block,
                            "baseline_method": baseline,
                            "preprocessing": preprocessing,
                            "distance_nm": float(value),
                            "low_nominal_speed_um_per_s": low_speed,
                            "high_nominal_speed_um_per_s": high_speed,
                            "low_actual_gap_speed_um_per_s": pair_metadata[
                                (block, low_speed)
                            ]["actual_gap_speed_um_per_s"],
                            "high_actual_gap_speed_um_per_s": pair_metadata[
                                (block, high_speed)
                            ]["actual_gap_speed_um_per_s"],
                            "low_force_mean_pN": float(np.mean(low_slice[complete])),
                            "high_force_mean_pN": float(np.mean(high_slice[complete])),
                            "low_force_median_pN": float(np.median(low_slice[complete])),
                            "high_force_median_pN": float(np.median(high_slice[complete])),
                            "pair_center_difference_high_minus_low_min": (
                                pair_metadata[(block, high_speed)][
                                    "pair_center_epoch_s"
                                ]
                                - pair_metadata[(block, low_speed)][
                                    "pair_center_epoch_s"
                                ]
                            )
                            / 60.0,
                            "claim_status": "conditional same-area paired comparison; pixels are not independent experimental repeats",
                        }
                        pixel_rows.append({**common, **summary})

                        full_difference = np.full(64, np.nan)
                        full_difference[complete] = difference
                        tile_values, minimum, maximum = tile_difference_means(
                            full_difference, coordinates
                        )
                        tile_summary = paired_t_summary(tile_values)
                        tile_rows.append(
                            {
                                **common,
                                **tile_summary,
                                "tile_size_pixels": TILE_SIZE_PIXELS,
                                "minimum_complete_pixels_per_retained_tile": minimum,
                                "maximum_complete_pixels_per_retained_tile": maximum,
                                "claim_status": "conditional same-area 2x2 spatial-tile paired comparison; not independent map replication",
                            }
                        )
    pixel = pd.DataFrame(pixel_rows)
    tile = pd.DataFrame(tile_rows)
    contrast_groups = [
        "preprocessing",
        "baseline_method",
        "block",
        "low_nominal_speed_um_per_s",
        "high_nominal_speed_um_per_s",
    ]
    family_groups = ["preprocessing", "baseline_method", "block"]
    for frame in (pixel, tile):
        add_multiplicity(
            frame,
            contrast_groups,
            "paired_t_p_two_sided",
            "across_37_distances",
        )
        add_multiplicity(
            frame,
            family_groups,
            "paired_t_p_two_sided",
            "across_all_contrasts_and_distances",
        )
    return pixel, tile


def experimental_unit_tests(
    distance: np.ndarray,
    force: dict[tuple[int, float, str, str], np.ndarray],
) -> pd.DataFrame:
    tests = {
        "map1_3_refresh_affected": (1, 2, 3),
        "map3_5_no_refresh": (3, 4, 5),
    }
    rows: list[dict] = []
    for preprocessing in PREPROCESSING:
        for baseline in BASELINES:
            for test, blocks in tests.items():
                for low_speed, high_speed in COMPARISONS:
                    for index, value in enumerate(distance):
                        differences: list[float] = []
                        included_blocks: list[int] = []
                        low_values: list[float] = []
                        high_values: list[float] = []
                        for block in blocks:
                            low_key = (block, low_speed, baseline, preprocessing)
                            high_key = (block, high_speed, baseline, preprocessing)
                            if low_key not in force or high_key not in force:
                                continue
                            low_slice = force[low_key][:, index]
                            high_slice = force[high_key][:, index]
                            complete = np.isfinite(low_slice) & np.isfinite(high_slice)
                            if np.count_nonzero(complete) < 56:
                                continue
                            low_summary = float(np.median(low_slice[complete]))
                            high_summary = float(np.median(high_slice[complete]))
                            low_values.append(low_summary)
                            high_values.append(high_summary)
                            differences.append(high_summary - low_summary)
                            included_blocks.append(block)
                        vector = np.asarray(differences, dtype=float)
                        if vector.size >= 3:
                            summary = paired_t_summary(vector)
                            sign_count, sign_p = exact_sign_test(vector)
                            status = "block_level_paired_t; n3_low_power"
                        else:
                            summary = paired_t_summary(np.asarray([], dtype=float))
                            summary["sample_count"] = int(vector.size)
                            sign_count, sign_p = exact_sign_test(vector)
                            status = "insufficient_complete_blocks_for_t_test"
                        rows.append(
                            {
                                "test": test,
                                "blocks_planned": ",".join(str(block) for block in blocks),
                                "blocks_included": ",".join(
                                    str(block) for block in included_blocks
                                ),
                                "baseline_method": baseline,
                                "preprocessing": preprocessing,
                                "distance_nm": float(value),
                                "low_nominal_speed_um_per_s": low_speed,
                                "high_nominal_speed_um_per_s": high_speed,
                                "block_summary": "spatial_median_of_F_sym_pN",
                                "low_block_values_pN": ";".join(
                                    f"{item:.12g}" for item in low_values
                                ),
                                "high_block_values_pN": ";".join(
                                    f"{item:.12g}" for item in high_values
                                ),
                                "block_differences_high_minus_low_pN": ";".join(
                                    f"{item:.12g}" for item in differences
                                ),
                                **summary,
                                "exact_sign_nonzero_block_count": sign_count,
                                "exact_sign_p_two_sided": sign_p,
                                "fit_status": status,
                                "claim_status": "block is the experimental replicate; Test A includes the refresh boundary and Test B has only two blocks for contrasts involving 0.1 um/s",
                            }
                        )
    frame = pd.DataFrame(rows)
    add_multiplicity(
        frame,
        ["test", "preprocessing", "baseline_method"],
        "paired_t_p_two_sided",
        "across_valid_contrasts_and_37_distances",
    )
    return frame


def save_symmetric_npz(
    distance: np.ndarray,
    force: dict[tuple[int, float, str, str], np.ndarray],
) -> None:
    keys = sorted(force)
    stack = np.stack([force[key] for key in keys], axis=0)
    np.savez_compressed(
        RESULTS / "palindrome_symmetric_pixel_forces.npz",
        distance_nm=distance,
        block=np.asarray([key[0] for key in keys], dtype=int),
        nominal_speed_um_per_s=np.asarray([key[1] for key in keys], dtype=float),
        baseline_method=np.asarray([key[2] for key in keys]),
        preprocessing=np.asarray([key[3] for key in keys]),
        force_pN=stack,
    )


def plot_force_distributions(
    distance: np.ndarray,
    force: dict[tuple[int, float, str, str], np.ndarray],
) -> None:
    fig, axes = plt.subplots(3, 4, figsize=(15.8, 11.0), constrained_layout=True)
    for row_index, block in enumerate((3, 4, 5)):
        for column, target in enumerate(TARGET_DISTANCES_NM):
            ax = axes[row_index, column]
            index = int(np.flatnonzero(distance == target)[0])
            speeds = [
                speed
                for speed in SPEEDS
                if (block, speed, "linear_drift_corrected", "palindrome_mean")
                in force
            ]
            values = [
                force[(block, speed, "linear_drift_corrected", "palindrome_mean")][
                    :, index
                ]
                for speed in speeds
            ]
            values = [item[np.isfinite(item)] for item in values]
            positions = np.arange(len(speeds), dtype=float)
            violin = ax.violinplot(
                values,
                positions=positions,
                widths=0.72,
                showmeans=False,
                showmedians=False,
                showextrema=False,
            )
            for body, speed in zip(violin["bodies"], speeds, strict=True):
                body.set_facecolor(COLORS[speed])
                body.set_edgecolor(COLORS[speed])
                body.set_alpha(0.28)
            for position, speed, item in zip(positions, speeds, values, strict=True):
                q25, median, q75 = np.quantile(item, [0.25, 0.5, 0.75])
                ax.vlines(position, q25, q75, color=COLORS[speed], lw=4)
                ax.scatter(position, median, color=COLORS[speed], s=28, zorder=3)
            ax.axhline(0.0, color="0.55", lw=0.8)
            ax.set_xticks(positions, [f"{speed:g}" for speed in speeds])
            ax.set_title(f"map{block}, D={target:.0f} nm")
            ax.grid(alpha=0.22, axis="y")
            if row_index == 2:
                ax.set_xlabel("Nominal approach speed (µm/s)")
            if column == 0:
                ax.set_ylabel("Palindrome-mean force (pN)")
            if block == 5:
                ax.text(
                    0.98,
                    0.96,
                    "0.1 late map missing",
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    fontsize=8,
                    color="0.35",
                )
    fig.suptitle(
        "Primary Test B: same-pixel force distributions after palindrome averaging\n"
        "far-linear baseline; point distributions are spatial observations, not independent map replicates",
        fontsize=16,
    )
    fig.savefig(FIGURES / "primary_palindrome_force_distributions.png", dpi=240)
    plt.close(fig)


def plot_endpoint_differences(
    distance: np.ndarray,
    force: dict[tuple[int, float, str, str], np.ndarray],
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.2), constrained_layout=True)
    for ax, target in zip(axes.flat, TARGET_DISTANCES_NM, strict=True):
        index = int(np.flatnonzero(distance == target)[0])
        for block_position, block in enumerate((3, 4, 5), start=1):
            for prep_index, preprocessing in enumerate(PREPROCESSING):
                low = force[(block, 0.05, "linear_drift_corrected", preprocessing)][
                    :, index
                ]
                high = force[(block, 0.2, "linear_drift_corrected", preprocessing)][
                    :, index
                ]
                valid = np.isfinite(low) & np.isfinite(high)
                difference = high[valid] - low[valid]
                position = block_position + (-0.14 if prep_index == 0 else 0.14)
                violin = ax.violinplot(
                    [difference],
                    positions=[position],
                    widths=0.24,
                    showmeans=False,
                    showmedians=False,
                    showextrema=False,
                )
                body = violin["bodies"][0]
                body.set_facecolor(PREP_COLORS[preprocessing])
                body.set_edgecolor(PREP_COLORS[preprocessing])
                body.set_alpha(0.30)
                q25, median, q75 = np.quantile(difference, [0.25, 0.5, 0.75])
                ax.vlines(
                    position,
                    q25,
                    q75,
                    color=PREP_COLORS[preprocessing],
                    lw=3,
                )
                ax.scatter(
                    position,
                    median,
                    color=PREP_COLORS[preprocessing],
                    s=25,
                    zorder=3,
                )
        ax.axhline(0.0, color="black", lw=1.0)
        ax.set_xticks((1, 2, 3), ("map3", "map4", "map5"))
        ax.set_title(f"D={target:.0f} nm")
        ax.set_ylabel("F_sym(0.2) − F_sym(0.05) (pN)")
        ax.grid(alpha=0.22, axis="y")
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color=PREP_COLORS[prep],
            label=prep.replace("_", " "),
            lw=0,
        )
        for prep in PREPROCESSING
    ]
    fig.legend(
        handles=handles,
        loc="outside lower center",
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    fig.suptitle(
        "Primary Test B: pixel-paired endpoint differences after palindrome averaging",
        fontsize=16,
    )
    fig.savefig(FIGURES / "primary_endpoint_paired_difference_distributions.png", dpi=240)
    plt.close(fig)


def plot_ttest_heatmap(distance: np.ndarray, pixel: pd.DataFrame, tile: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14.2, 7.8), constrained_layout=True)
    levels = (("64 pixels", pixel), ("16 spatial 2×2 tiles", tile))
    for row, preprocessing in enumerate(PREPROCESSING):
        for column, (label, frame) in enumerate(levels):
            ax = axes[row, column]
            values = np.full((3, distance.size), np.nan)
            significant = np.zeros_like(values, dtype=bool)
            for block_index, block in enumerate((3, 4, 5)):
                selected = frame[
                    (frame.preprocessing == preprocessing)
                    & (frame.baseline_method == "linear_drift_corrected")
                    & (frame.block == block)
                    & (frame.low_nominal_speed_um_per_s == 0.05)
                    & (frame.high_nominal_speed_um_per_s == 0.2)
                ].sort_values("distance_nm")
                if len(selected) != distance.size:
                    raise AssertionError("endpoint t-test heatmap row incomplete")
                values[block_index] = selected.paired_t_statistic.to_numpy(dtype=float)
                significant[block_index] = (
                    selected[
                        "paired_t_p_two_sided_holm_across_37_distances"
                    ].to_numpy(dtype=float)
                    < ALPHA
                )
            clipped = np.clip(values, -12.0, 12.0)
            image = ax.imshow(
                clipped,
                aspect="auto",
                origin="upper",
                cmap="coolwarm",
                norm=TwoSlopeNorm(vmin=-12.0, vcenter=0.0, vmax=12.0),
                extent=(distance[0] - 2.5, distance[-1] + 2.5, 2.5, -0.5),
            )
            y, x = np.where(significant)
            ax.scatter(distance[x], y, s=12, color="black", marker="o")
            ax.set_yticks((0, 1, 2), ("map3", "map4", "map5"))
            ax.set_xticks(TARGET_DISTANCES_NM)
            ax.set_xlabel("Separation D (nm)")
            ax.set_title(f"{preprocessing.replace('_', ' ')} — {label}")
            if column == 0:
                ax.set_ylabel("Palindrome block")
    colorbar = fig.colorbar(image, ax=axes, shrink=0.92, pad=0.02)
    colorbar.set_label("paired t statistic for F(0.2) − F(0.05), clipped to ±12")
    fig.suptitle(
        "Conditional within-block t tests; black dots: Holm p<0.05 across 37 distances",
        fontsize=15,
    )
    fig.savefig(FIGURES / "primary_endpoint_ttest_heatmap.png", dpi=240)
    plt.close(fig)


def plot_experimental_unit_tests(experimental: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.2), constrained_layout=True)
    for ax, target in zip(axes.flat, TARGET_DISTANCES_NM, strict=True):
        positions = {
            ("map1_3_refresh_affected", "palindrome_mean"): 0.85,
            ("map1_3_refresh_affected", "linear_center_aligned"): 1.15,
            ("map3_5_no_refresh", "palindrome_mean"): 1.85,
            ("map3_5_no_refresh", "linear_center_aligned"): 2.15,
        }
        for (test, preprocessing), position in positions.items():
            selected = experimental[
                (experimental.test == test)
                & (experimental.preprocessing == preprocessing)
                & (experimental.baseline_method == "linear_drift_corrected")
                & (experimental.distance_nm == target)
                & (experimental.low_nominal_speed_um_per_s == 0.05)
                & (experimental.high_nominal_speed_um_per_s == 0.2)
            ]
            if len(selected) != 1:
                raise AssertionError("experimental-unit plot selection failed")
            item = selected.iloc[0]
            difference = np.asarray(
                [
                    float(value)
                    for value in item.block_differences_high_minus_low_pN.split(";")
                ]
            )
            ax.scatter(
                np.full(difference.size, position),
                difference,
                s=32,
                color=PREP_COLORS[preprocessing],
                alpha=0.75,
            )
            if np.isfinite(item.mean_difference_pN):
                ax.errorbar(
                    position,
                    item.mean_difference_pN,
                    yerr=[
                        [item.mean_difference_pN - item.mean_95CI_low_pN],
                        [item.mean_95CI_high_pN - item.mean_difference_pN],
                    ],
                    fmt="D",
                    color=PREP_COLORS[preprocessing],
                    capsize=4,
                    ms=5,
                    zorder=4,
                )
                ax.text(
                    position,
                    item.mean_95CI_high_pN,
                    f"p={p_text(float(item.paired_t_p_two_sided))}",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                )
        ax.axhline(0.0, color="black", lw=1.0)
        ax.set_xticks((1.0, 2.0), ("Test A\nmap1–3", "Test B\nmap3–5"))
        ax.set_title(f"D={target:.0f} nm")
        ax.set_ylabel("Block spatial-median ΔF, 0.2−0.05 (pN)")
        ax.grid(alpha=0.22, axis="y")
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="D",
            color=PREP_COLORS[prep],
            label=prep.replace("_", " "),
        )
        for prep in PREPROCESSING
    ]
    fig.legend(
        handles=handles,
        loc="outside lower center",
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    fig.suptitle(
        "Experimental-unit paired t tests: three block summaries per complete endpoint contrast",
        fontsize=15,
    )
    fig.savefig(FIGURES / "experimental_unit_endpoint_ttests.png", dpi=240)
    plt.close(fig)


def plot_pair_timing(timing: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.8), constrained_layout=True)
    for block in range(1, 6):
        selected = timing[timing.block == block].sort_values(
            "nominal_speed_um_per_s"
        )
        axes[0].plot(
            selected.nominal_speed_um_per_s,
            selected.pair_center_offset_from_block_mean_min,
            marker="o",
            label=f"map{block}",
        )
        axes[1].plot(
            selected.nominal_speed_um_per_s,
            selected.pair_half_span_min,
            marker="o",
            label=f"map{block}",
        )
    axes[0].axhline(0.0, color="black", lw=0.8)
    axes[0].set_ylabel("Pair-center offset from block mean (min)")
    axes[1].set_ylabel("Early/late half-span (min)")
    for ax in axes:
        ax.set_xlabel("Nominal approach speed (µm/s)")
        ax.set_xticks(SPEEDS)
        ax.grid(alpha=0.25)
    axes[0].set_title("Clock-time centers are close but not identical")
    axes[1].set_title("Palindrome history spans differ by speed/order")
    axes[1].legend(ncol=2, fontsize=8)
    fig.savefig(FIGURES / "palindrome_pair_clock_timing.png", dpi=240)
    plt.close(fig)


def selected_summary(
    pixel: pd.DataFrame,
    tile: pd.DataFrame,
    experimental: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict] = []
    for preprocessing in PREPROCESSING:
        for baseline in BASELINES:
            for target in TARGET_DISTANCES_NM:
                block_test = experimental[
                    (experimental.test == "map3_5_no_refresh")
                    & (experimental.preprocessing == preprocessing)
                    & (experimental.baseline_method == baseline)
                    & (experimental.distance_nm == target)
                    & (experimental.low_nominal_speed_um_per_s == 0.05)
                    & (experimental.high_nominal_speed_um_per_s == 0.2)
                ].iloc[0]
                for block in (3, 4, 5):
                    pix = pixel[
                        (pixel.preprocessing == preprocessing)
                        & (pixel.baseline_method == baseline)
                        & (pixel.block == block)
                        & (pixel.distance_nm == target)
                        & (pixel.low_nominal_speed_um_per_s == 0.05)
                        & (pixel.high_nominal_speed_um_per_s == 0.2)
                    ].iloc[0]
                    til = tile[
                        (tile.preprocessing == preprocessing)
                        & (tile.baseline_method == baseline)
                        & (tile.block == block)
                        & (tile.distance_nm == target)
                        & (tile.low_nominal_speed_um_per_s == 0.05)
                        & (tile.high_nominal_speed_um_per_s == 0.2)
                    ].iloc[0]
                    rows.append(
                        {
                            "test": "map3_5_no_refresh",
                            "preprocessing": preprocessing,
                            "baseline_method": baseline,
                            "distance_nm": target,
                            "block": block,
                            "contrast": "0.2_minus_0.05_um_per_s",
                            "pixel_mean_difference_pN": pix.mean_difference_pN,
                            "pixel_paired_t_p": pix.paired_t_p_two_sided,
                            "pixel_holm_p_across_37_distances": pix[
                                "paired_t_p_two_sided_holm_across_37_distances"
                            ],
                            "tile_mean_difference_pN": til.mean_difference_pN,
                            "tile_paired_t_p": til.paired_t_p_two_sided,
                            "tile_holm_p_across_37_distances": til[
                                "paired_t_p_two_sided_holm_across_37_distances"
                            ],
                            "experimental_block_count": block_test.sample_count,
                            "experimental_mean_of_block_median_differences_pN": block_test.mean_difference_pN,
                            "experimental_paired_t_p": block_test.paired_t_p_two_sided,
                            "experimental_holm_p": block_test[
                                "paired_t_p_two_sided_holm_across_valid_contrasts_and_37_distances"
                            ],
                            "interpretation": "pixel/tile tests are conditional within measured area; block test is the experimental-unit inference",
                        }
                    )
    return pd.DataFrame(rows)


def build_report(
    timing: pd.DataFrame,
    alignment: pd.DataFrame,
    pixel: pd.DataFrame,
    tile: pd.DataFrame,
    experimental: pd.DataFrame,
) -> str:
    valid_block_tests = experimental[
        np.isfinite(experimental.paired_t_p_two_sided)
    ]
    adjusted_block_field = (
        "paired_t_p_two_sided_holm_across_valid_contrasts_and_37_distances"
    )
    test_b_primary = valid_block_tests[
        (valid_block_tests.test == "map3_5_no_refresh")
        & (valid_block_tests.preprocessing == "palindrome_mean")
        & (valid_block_tests.baseline_method == "linear_drift_corrected")
        & (valid_block_tests.low_nominal_speed_um_per_s == 0.05)
        & (valid_block_tests.high_nominal_speed_um_per_s == 0.2)
    ]
    test_b_minimum = test_b_primary.loc[
        test_b_primary.paired_t_p_two_sided.idxmin()
    ]
    lines = [
        "# 27-08-26 回文平均后的速度分布与paired t-test",
        "",
        "## 直接结论",
        "",
        "回文平均后，若把64个pixel当样本，部分单独block和距离的速度分布确实给出很小的paired t-test p-value；但空间tile检验显著减少，而且map3/map4/map5的速度差会变号。以block作为真正实验重复时，Test B在20、50、100、200 nm均不显著；在完整20–200 nm扫描中最小raw p也只有 "
        f"`{test_b_minimum.paired_t_p_two_sided:.4f}`（{test_b_minimum.distance_nm:.0f} nm），Holm后为 `{test_b_minimum[adjusted_block_field]:.4f}`。全部{len(valid_block_tests)}个可计算的n=3 block-level sensitivity tests中，没有一个Holm-adjusted p<0.05。因此当前证据是 **conditional within-map separation存在，但replicated velocity separation未建立**。",
        "",
        "## 问题与统计单位",
        "",
        "对每个block、速度和同一physical pixel先计算 `F_sym=(F_early+F_late)/2`。如果时间漂移在early/late之间近似线性，且各速度pair具有共同中心时刻，该平均会消除一阶时间项。t-test使用同位置差值 `F_sym(v_high)-F_sym(v_low)`，因此是paired而不是Welch independent-sample test。",
        "",
        "结果分三层：",
        "",
        "1. 64-pixel paired t-test：回答本次已测区域内同位置值是否系统变化，但pixels空间相关，会低估standard error。",
        "2. 16个2×2 spatial-tile paired t-test：先对局部tile内差值取均值，作为较保守的区域内分离检验。",
        "3. block-level paired t-test：每个block/速度先取64-pixel spatial median，再以三个block作为n=3实验重复。这一层才对应跨map的速度推断，但自由度只有2。",
        "",
        "paired t-test检验的是mean paired difference是否为零，不等价于检验两个完整distribution的所有形状都不同。",
        "",
        "## 回文在clock time上的对称程度",
        "",
        "| block | complete speeds | pair-center spread (min) | largest early/late half-span (min) |",
        "|---:|---:|---:|---:|",
    ]
    for block in range(1, 6):
        selected = timing[timing.block == block]
        speeds = ", ".join(
            f"{value:g}" for value in selected.nominal_speed_um_per_s
        )
        lines.append(
            f"| map{block} | {speeds} | {selected.block_pair_center_spread_min.iloc[0]:.3f} | {selected.pair_half_span_min.max():.2f} |"
        )
    lines.extend(
        [
            "",
            "map3的三种速度pair center相差约4 min，因此原始回文平均并不严格对应同一clock time。`linear_center_aligned` sensitivity branch用所有同block early/late差值估计一个共同线性drift，再把每个F_sym平移到该block的平均pair center；它只能检验一阶中心错位，不能消除已观察到的非线性relaxation。",
            "",
            "## Primary Test B：0.2与0.05 µm/s",
            "",
            "以下采用far-linear baseline和未经额外对齐的palindrome mean。ΔF为0.2−0.05 µm/s。pixel/tile的Holm p在每个block和contrast内跨20–200 nm的37个距离校正；最后一列以map3/map4/map5三个block spatial median作paired t-test。",
            "",
            "| D | map3 Δmean / pixel Holm / tile Holm | map4 | map5 | n=3 block-median mean Δ / raw p / Holm p |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for target in TARGET_DISTANCES_NM:
        block_cells = []
        for block in (3, 4, 5):
            pix = pixel[
                (pixel.preprocessing == "palindrome_mean")
                & (pixel.baseline_method == "linear_drift_corrected")
                & (pixel.block == block)
                & (pixel.distance_nm == target)
                & (pixel.low_nominal_speed_um_per_s == 0.05)
                & (pixel.high_nominal_speed_um_per_s == 0.2)
            ].iloc[0]
            til = tile[
                (tile.preprocessing == "palindrome_mean")
                & (tile.baseline_method == "linear_drift_corrected")
                & (tile.block == block)
                & (tile.distance_nm == target)
                & (tile.low_nominal_speed_um_per_s == 0.05)
                & (tile.high_nominal_speed_um_per_s == 0.2)
            ].iloc[0]
            block_cells.append(
                f"{pix.mean_difference_pN:+.2f} pN / {p_text(pix.paired_t_p_two_sided_holm_across_37_distances)} / {p_text(til.paired_t_p_two_sided_holm_across_37_distances)}"
            )
        unit = experimental[
            (experimental.test == "map3_5_no_refresh")
            & (experimental.preprocessing == "palindrome_mean")
            & (experimental.baseline_method == "linear_drift_corrected")
            & (experimental.distance_nm == target)
            & (experimental.low_nominal_speed_um_per_s == 0.05)
            & (experimental.high_nominal_speed_um_per_s == 0.2)
        ].iloc[0]
        lines.append(
            f"| {target:.0f} nm | {block_cells[0]} | {block_cells[1]} | {block_cells[2]} | {unit.mean_difference_pN:+.2f} pN / {p_text(unit.paired_t_p_two_sided)} / {p_text(unit.paired_t_p_two_sided_holm_across_valid_contrasts_and_37_distances)} |"
        )

    lines.extend(
        [
            "",
            "逐block pixel/tile检验即使显著，也只说明该block、该区域的paired distributions具有非零mean shift。真正的block-level结果由三个Δ值组成；如果它们变号或离散很大，64 pixels不能增加独立速度重复数。n=3 exact sign test即使3/3同号，two-sided最小p也为0.25。",
            "",
            "Test A的far-linear palindrome endpoint在20 nm给出raw block-level p=0.0103，但它跨越换液干预，而且在预先保留的全距离/全contrast Holm family中p=1.000；far-constant对应raw p=0.0045、Holm p=0.499。它最多是下一轮实验的候选距离，不能作为当前速度因果结论。",
            "",
            "## 显著距离范围",
            "",
            "下表只看Test B的0.2−0.05 endpoint、far-linear baseline。`+`表示high speed force更大，`−`相反；Holm correction跨37个距离。",
            "",
            "| preprocessing | block | 64-pixel paired t | 16-tile paired t |",
            "|---|---:|---|---|",
        ]
    )
    for preprocessing in PREPROCESSING:
        for block in (3, 4, 5):
            pix = pixel[
                (pixel.preprocessing == preprocessing)
                & (pixel.baseline_method == "linear_drift_corrected")
                & (pixel.block == block)
                & (pixel.low_nominal_speed_um_per_s == 0.05)
                & (pixel.high_nominal_speed_um_per_s == 0.2)
            ]
            til = tile[
                (tile.preprocessing == preprocessing)
                & (tile.baseline_method == "linear_drift_corrected")
                & (tile.block == block)
                & (tile.low_nominal_speed_um_per_s == 0.05)
                & (tile.high_nominal_speed_um_per_s == 0.2)
            ]
            lines.append(
                f"| {preprocessing} | map{block} | {significant_segments(pix, 'paired_t_p_two_sided_holm_across_37_distances')} | {significant_segments(til, 'paired_t_p_two_sided_holm_across_37_distances')} |"
            )

    lines.extend(
        [
            "",
            "## 线性中心时刻对齐和baseline敏感性",
            "",
        ]
    )
    for target in TARGET_DISTANCES_NM:
        selected = alignment[
            (alignment.baseline_method == "linear_drift_corrected")
            & (alignment.distance_nm == target)
            & alignment.block.isin((3, 4, 5))
        ]
        lines.append(
            f"- {target:.0f} nm：map3–5的median absolute clock-center adjustment为 "
            + ", ".join(
                f"map{int(row.block)} {row.absolute_center_adjustment_median_pN:.2f} pN"
                for row in selected.itertuples(index=False)
            )
            + "。"
        )
    lines.extend(
        [
            "",
            "原始palindrome mean、linear-center alignment、far-linear和far-constant四种组合都保留在CSV。若显著区间或ΔF方向随这些合理处理改变，则不能称为robust velocity separation。特别是center alignment仍假设每个block存在一个共同线性drift；当前chronology显示relaxation具有曲率，而且该branch的t-test没有传播drift-estimation uncertainty，所以它是敏感性分析而不是修正后的ground truth。",
            "",
            "## 结论",
            "",
            "- 回文平均确实降低了同速度early/late的一阶时间偏差，并允许在同一block内做严格same-pixel paired comparison。",
            "- pixel-level paired t-test回答的是已测空间区域内是否分离；它不能把64个pixel变成64次独立实验。2×2 tile结果更保守，但仍不是独立map replication。",
            "- 速度因果的主统计单位是block。Test B的0.05–0.2 endpoint只有n=3；涉及0.1 µm/s的Test B block-level比较只有map3/map4两组，未执行正式t-test。",
            "- 因pair center并非完全同时、relaxation又明显非线性，回文平均不能保证彻底消除时间效应。是否可主张velocity separation必须以n=3 block-level一致性和处理敏感性为准，而不能只看很小的pixel p-value。",
            "",
            "## 输出",
            "",
            "- `palindrome_pixel_paired_ttests.csv`: 64-pixel paired t-tests and multiplicity corrections.",
            "- `palindrome_spatial_tile_paired_ttests.csv`: sixteen 2x2-tile paired t-tests.",
            "- `palindrome_experimental_unit_ttests.csv`: n=3 block-level tests and n<3 exclusions.",
            "- `pair_center_timing.csv`, `linear_center_alignment_QC.csv`: clock-time symmetry and alignment sensitivity.",
            "- `palindrome_symmetric_pixel_forces.npz`: all F_sym pixel distributions used by the tests.",
            "- `selected_target_summary.csv`, `figures/`: target-distance tables and visualizations.",
            "- `provenance.json`, `artifact_manifest.sha256`: upstream identities, definitions, software and hashes.",
        ]
    )
    return "\n".join(lines) + "\n"


def write_manifest(paths: list[Path]) -> None:
    lines = []
    for path in sorted(paths):
        lines.append(f"{sha256_file(path)}  {path.relative_to(ROOT)}")
    (RESULTS / "artifact_manifest.sha256").write_text("\n".join(lines) + "\n")


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    (
        distance,
        sources,
        arrays,
        map_frame,
        pair_frame,
        coordinates,
        upstream_provenance,
    ) = load_inputs()
    force, pair_metadata, timing, alignment = build_palindrome_forces(
        distance, sources, arrays, map_frame, pair_frame
    )
    pixel, tile = distribution_tests(
        distance, force, pair_metadata, coordinates
    )
    experimental = experimental_unit_tests(distance, force)
    selected = selected_summary(pixel, tile, experimental)

    timing.to_csv(RESULTS / "pair_center_timing.csv", index=False)
    alignment.to_csv(RESULTS / "linear_center_alignment_QC.csv", index=False)
    pixel.to_csv(RESULTS / "palindrome_pixel_paired_ttests.csv", index=False)
    tile.to_csv(RESULTS / "palindrome_spatial_tile_paired_ttests.csv", index=False)
    experimental.to_csv(
        RESULTS / "palindrome_experimental_unit_ttests.csv", index=False
    )
    selected.to_csv(RESULTS / "selected_target_summary.csv", index=False)
    save_symmetric_npz(distance, force)

    plot_force_distributions(distance, force)
    plot_endpoint_differences(distance, force)
    plot_ttest_heatmap(distance, pixel, tile)
    plot_experimental_unit_tests(experimental)
    plot_pair_timing(timing)

    report = build_report(timing, alignment, pixel, tile, experimental)
    (RESULTS / "REPORT.md").write_text(report)

    upstream_files = [
        UPSTREAM / "pixel_force_curves.npz",
        UPSTREAM / "map_inventory_QC.csv",
        UPSTREAM / "pair_inventory.csv",
        UPSTREAM / "pixel_QC.csv",
        UPSTREAM / "provenance.json",
    ]
    provenance = {
        "analysis": "same-pixel palindrome-mean velocity distribution paired t-tests",
        "temperature_C": upstream_provenance["temperature_C"],
        "map_environment": upstream_provenance["map_environment"],
        "geometry": upstream_provenance["geometry"],
        "same_region": upstream_provenance["same_region"],
        "distance_range_nm": [DISTANCE_MIN_NM, DISTANCE_MAX_NM],
        "distance_step_nm": 5.0,
        "target_distances_nm": TARGET_DISTANCES_NM,
        "primary_transformation": "F_sym=(F_early+F_late)/2 at the same physical pixel",
        "sensitivity_transformation": "common within-block linear drift inferred from early/late differences and F_sym shifted to mean pair-center clock time",
        "force_baselines": BASELINES,
        "statistical_levels": {
            "pixel": "64 same-position paired observations; spatial pseudoreplication caveat",
            "tile": "sixteen 2x2-pixel mean paired differences; conditional spatial inference",
            "block": "three block-level spatial medians; experimental-unit paired t-test",
        },
        "multiplicity": {
            "within_contrast": "Holm and BH across 37 distances",
            "within_block_family": "Holm and BH across all speed contrasts and 37 distances",
            "experimental_unit": "Holm and BH across all valid contrasts and 37 distances within test/baseline/preprocessing",
        },
        "missing_data": upstream_provenance["missing_data"],
        "time_coordinate": upstream_provenance["time_coordinate"],
        "upstream_input_hashes": {
            str(path.relative_to(ROOT)): sha256_file(path) for path in upstream_files
        },
        "upstream_raw_input_hashes": upstream_provenance["input_hashes"],
        "code_hash": sha256_file(Path(__file__).resolve()),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "randomness": "none",
        "claim_boundary": "paired t tests assess mean shifts; pixel/tile results are conditional within-area evidence and only block-level n=3 represents map replication",
    }
    (RESULTS / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n"
    )

    artifacts = [
        Path(__file__).resolve(),
        *[path for path in RESULTS.glob("*") if path.is_file() and path.name != "artifact_manifest.sha256"],
        *[path for path in FIGURES.glob("*.png")],
    ]
    write_manifest(artifacts)
    print(f"Wrote {RESULTS}")
    print(f"pixel tests={len(pixel)}, tile tests={len(tile)}, block tests={len(experimental)}")


if __name__ == "__main__":
    main()

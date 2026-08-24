#!/usr/bin/env python3
"""Conditional distribution-separation analysis for the three 0 wt% maps.

The three velocity maps cover the same nominal 16 x 16 raster but were
acquired sequentially (2, then 1, then 4 um/s).  This analysis therefore uses
same-coordinate paired force differences and spatial block resampling.  Its
confidence statements describe separation among these recorded maps over the
measured area; they are not independent-map confidence intervals for a causal
velocity effect.

For each faster-minus-slower comparison and each 5 nm separation bin from
25--250 nm, the script reports:

* the paired median force difference;
* the paired dominance fraction P(F_fast > F_slow);
* pointwise 95% spatial-block bootstrap intervals;
* a 95% studentized simultaneous band across all three comparisons and all
  distances; and
* an exact 4 x 4-block sign-flip max-T test, Holm-adjusted across the three
  pairwise comparisons.

Only map pixels enter the inferential calculation.  The five stand-alone
0 wt% force curves remain hard-contact sensitivity anchors but are not balanced
across speeds and have no 4 um/s counterpart, so treating them as replicate
map pixels would be pseudoreplication.
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

import analyze_velocity_joint_fit as joint
import analyze_velocity_systematics as events
import fit_glycerol_surface_forces as base


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "20-08-26" / "0"
RESULTS_DIR = ROOT / "analysis" / "velocity_joint_fit_results"
FIGURES_DIR = RESULTS_DIR / "figures"
SENSITIVITY_CSV = ROOT / "analysis" / "surface_force_results" / "sensitivity_by_source.csv"

DISTANCE_MIN_NM = 25.0
DISTANCE_MAX_NM = 250.0
SELECTED_DISTANCES_NM = (25.0, 50.0, 100.0, 150.0, 200.0, 250.0)
SPEEDS_UM_PER_S = (1, 2, 4)
COMPARISONS = ((1, 2), (2, 4), (1, 4))
PRIMARY_BLOCK_SIZE_PIXELS = 4
BOOTSTRAP_REPLICATES = 10_000
RANDOM_SEED = 20260824

DETAIL_CSV = RESULTS_DIR / "zero_wt_distribution_separation.csv"
SELECTED_CSV = RESULTS_DIR / "zero_wt_distribution_separation_selected.csv"
GLOBAL_CSV = RESULTS_DIR / "zero_wt_distribution_separation_global.csv"
REPORT_MD = RESULTS_DIR / "ZERO_WT_DISTRIBUTION_SEPARATION.md"
PROVENANCE_JSON = RESULTS_DIR / "zero_wt_distribution_separation_provenance.json"
FIGURE_STEM = FIGURES_DIR / "zero_wt_distribution_separation_confidence"


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


def curve_key(curve: base.RawCurve | events.BranchCurve, order: int) -> int:
    return int(curve.point_index) if curve.point_index is not None else int(order)


def calibrate_zero_wt_sources(sources: list[base.SourceData]) -> float:
    """Rebuild 0 wt% contact anchors and verify the saved consensus InvOLS."""

    for source in sources:
        values: list[float] = []
        for curve in source.curves:
            try:
                curve.far_field_fit = base.fit_far_field_drift(
                    curve.measured_height_m,
                    curve.deflection_V,
                )
            except (ValueError, np.linalg.LinAlgError):
                curve.far_field_fit = None
            curve.contact_fit = base.find_contact_fit(
                curve.measured_height_m,
                curve.deflection_V,
                curve.far_field_fit,
            )
            if curve.contact_fit is not None:
                values.append(curve.contact_fit.sensitivity_m_per_V)
        array = np.asarray(values, dtype=float)
        if array.size:
            median = float(np.median(array))
            mad = base.robust_mad(array)
            tolerance = max(3.5 * mad if np.isfinite(mad) else 0.0, 5e-9)
            retained = array[np.abs(array - median) <= tolerance]
        else:
            retained = array
        required = 1 if source.source_type == "force" else max(
            10,
            int(math.ceil(0.03 * len(source.curves))),
        )
        source.sensitivity_valid_curves = int(retained.size)
        if retained.size >= required:
            source.sensitivity_anchor_m_per_V = float(np.median(retained))
            source.sensitivity_anchor_mad_m_per_V = base.robust_mad(retained)

    anchors = np.asarray(
        [
            source.sensitivity_anchor_m_per_V
            for source in sources
            if np.isfinite(source.sensitivity_anchor_m_per_V)
        ],
        dtype=float,
    )
    if anchors.size < 2:
        raise RuntimeError("fewer than two reconstructed 0 wt% sensitivity anchors")
    center = float(np.median(anchors))
    retained_anchors = anchors[np.abs(anchors - center) <= 15e-9]
    reconstructed_consensus = float(np.median(retained_anchors))

    with SENSITIVITY_CSV.open(newline="", encoding="utf-8") as handle:
        saved_values = {
            float(row["sensitivity_used_nm_per_V"]) * 1e-9
            for row in csv.DictReader(handle)
            if int(row["concentration_wt_percent"]) == 0
        }
    if len(saved_values) != 1:
        raise RuntimeError(f"expected one saved 0 wt% sensitivity, got {saved_values}")
    saved_consensus = saved_values.pop()
    if not np.isclose(
        reconstructed_consensus,
        saved_consensus,
        rtol=0.0,
        atol=1e-15,
    ):
        raise RuntimeError(
            "reconstructed and saved 0 wt% sensitivities disagree: "
            f"{reconstructed_consensus * 1e9:.9g} vs "
            f"{saved_consensus * 1e9:.9g} nm/V"
        )
    for source in sources:
        source.sensitivity_used_m_per_V = saved_consensus
        source.sensitivity_method = "concentration_consensus_raw_contacts"
    return saved_consensus


def load_zero_wt_maps() -> tuple[list[base.SourceData], dict[int, joint.MapData]]:
    """Decode only 0 wt% raw sources using the established reconstruction."""

    paths = sorted(DATA_DIR.glob("*.jpk-force"))
    paths += sorted(DATA_DIR.glob("*.jpk-force-map"))
    sources = [base.load_source(path, 0) for path in paths]
    sources.sort(key=lambda item: item.timestamp)
    calibrate_zero_wt_sources(sources)

    maps_by_speed: dict[int, joint.MapData] = {}
    for source in sources:
        if source.source_type != "map":
            continue
        raw_approaches, _ = events.load_branch(source.path, "extend")
        retracts, _ = events.load_branch(source.path, "retract")
        raw_by_key = {
            curve_key(curve, order): curve
            for order, curve in enumerate(raw_approaches)
        }
        retract_by_key = {
            curve_key(curve, order): curve
            for order, curve in enumerate(retracts)
        }
        points: dict[tuple[int, int], joint.PointCurve] = {}
        for order, approach in enumerate(source.curves):
            key = curve_key(approach, order)
            retract = retract_by_key.get(key)
            try:
                event_record = events.analyze_pair(
                    approach,
                    retract,
                    source.sensitivity_used_m_per_V,
                )
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                event_record = events._nan_record()
            point = joint.prepare_point_curve(
                source,
                approach,
                raw_by_key.get(key),
                retract,
                order,
                event_record,
            )
            if point is not None:
                points[(point.row, point.column)] = point
        if not points:
            raise RuntimeError(f"no valid points reconstructed for {source.path}")

        terminal = joint.safe_median(
            point.approach_terminal_load_nN for point in points.values()
        )
        measured_speed = joint.safe_median(
            point.approach_speed_um_per_s for point in points.values()
        )
        label = int(round(joint.speed_label(measured_speed)))
        if label in maps_by_speed:
            raise RuntimeError(f"more than one 0 wt% map at {label} um/s")
        maps_by_speed[label] = joint.MapData(
            source=str(source.path.relative_to(ROOT)),
            concentration_wt_percent=0,
            timestamp=source.timestamp.isoformat(),
            speed_um_per_s=measured_speed,
            speed_label_um_per_s=float(label),
            load_regime=joint.classify_load_regime(0, terminal),
            primary_included=True,
            sensitivity_nm_per_V=source.sensitivity_used_m_per_V * 1e9,
            contact_response_center=joint.safe_median(
                point.approach_contact_response_ratio for point in points.values()
            ),
            contact_response_mad=base.robust_mad(
                np.asarray(
                    [point.approach_contact_response_ratio for point in points.values()]
                )
            ),
            contact_invOLS_center_nm_per_V=joint.safe_median(
                point.approach_contact_invOLS_nm_per_V for point in points.values()
            ),
            terminal_load_center_nN=terminal,
            far_slope_center_pN_per_100nm=joint.safe_median(
                point.approach_far_slope_pN_per_100nm for point in points.values()
            ),
            points=points,
        )

    if set(maps_by_speed) != set(SPEEDS_UM_PER_S):
        raise RuntimeError(
            f"expected speeds {SPEEDS_UM_PER_S}, got {sorted(maps_by_speed)}"
        )
    return sources, maps_by_speed


def holm_adjust(p_values: list[float]) -> list[float]:
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    adjusted = np.empty(values.shape, dtype=float)
    running = 0.0
    count = values.size
    for rank, index in enumerate(order):
        candidate = min(1.0, (count - rank) * values[index])
        running = max(running, candidate)
        adjusted[index] = running
    return [float(value) for value in adjusted]


def block_members(
    coordinates: list[tuple[int, int]], block_size: int
) -> tuple[np.ndarray, list[np.ndarray]]:
    labels = np.asarray(
        [
            (row // block_size) * 100 + column // block_size
            for row, column in coordinates
        ],
        dtype=np.int64,
    )
    unique = np.unique(labels)
    members = [np.flatnonzero(labels == label) for label in unique]
    if any(index.size == 0 for index in members):
        raise AssertionError("empty spatial block")
    return unique, members


def paired_dominance(difference: np.ndarray, axis: int) -> np.ndarray:
    return np.mean(difference > 0.0, axis=axis) + 0.5 * np.mean(
        difference == 0.0,
        axis=axis,
    )


def bootstrap_functional_statistics(
    forces: dict[int, np.ndarray],
    coordinates: list[tuple[int, int]],
    rng: np.random.Generator,
) -> dict[str, np.ndarray | float | int]:
    """Joint spatial-block bootstrap for all pairwise functional statistics."""

    _, members = block_members(coordinates, PRIMARY_BLOCK_SIZE_PIXELS)
    pair_count = len(COMPARISONS)
    distance_count = next(iter(forces.values())).shape[1]
    estimate_difference = np.empty((pair_count, distance_count), dtype=float)
    estimate_dominance = np.empty_like(estimate_difference)
    boot_difference = np.empty(
        (BOOTSTRAP_REPLICATES, pair_count, distance_count), dtype=float
    )
    boot_dominance = np.empty_like(boot_difference)

    paired_arrays: list[np.ndarray] = []
    for pair_index, (slow, fast) in enumerate(COMPARISONS):
        difference = forces[fast] - forces[slow]
        paired_arrays.append(difference)
        estimate_difference[pair_index] = np.median(difference, axis=0)
        estimate_dominance[pair_index] = paired_dominance(difference, axis=0)

    for replicate in range(BOOTSTRAP_REPLICATES):
        selected_blocks = rng.integers(0, len(members), size=len(members))
        selected_pixels = np.concatenate([members[index] for index in selected_blocks])
        for pair_index, difference in enumerate(paired_arrays):
            sample = difference[selected_pixels]
            boot_difference[replicate, pair_index] = np.median(sample, axis=0)
            boot_dominance[replicate, pair_index] = paired_dominance(sample, axis=0)

    difference_pointwise_low, difference_pointwise_high = np.quantile(
        boot_difference,
        (0.025, 0.975),
        axis=0,
    )
    dominance_pointwise_low, dominance_pointwise_high = np.quantile(
        boot_dominance,
        (0.025, 0.975),
        axis=0,
    )

    difference_se = np.std(boot_difference, axis=0, ddof=1)
    difference_valid = difference_se > np.finfo(float).eps
    difference_studentized = np.zeros_like(boot_difference)
    np.divide(
        boot_difference - estimate_difference,
        difference_se,
        out=difference_studentized,
        where=difference_valid,
    )
    difference_max_t = np.max(np.abs(difference_studentized), axis=(1, 2))
    difference_critical = float(np.quantile(difference_max_t, 0.95))
    difference_simultaneous_low = estimate_difference - difference_critical * difference_se
    difference_simultaneous_high = estimate_difference + difference_critical * difference_se

    dominance_se = np.std(boot_dominance, axis=0, ddof=1)
    dominance_valid = dominance_se > np.finfo(float).eps
    dominance_studentized = np.zeros_like(boot_dominance)
    np.divide(
        boot_dominance - estimate_dominance,
        dominance_se,
        out=dominance_studentized,
        where=dominance_valid,
    )
    dominance_max_t = np.max(np.abs(dominance_studentized), axis=(1, 2))
    dominance_critical = float(np.quantile(dominance_max_t, 0.95))
    dominance_simultaneous_low = np.clip(
        estimate_dominance - dominance_critical * dominance_se,
        0.0,
        1.0,
    )
    dominance_simultaneous_high = np.clip(
        estimate_dominance + dominance_critical * dominance_se,
        0.0,
        1.0,
    )

    if not np.all(difference_pointwise_low <= difference_pointwise_high):
        raise AssertionError("reversed force-difference bootstrap interval")
    if not np.all(dominance_pointwise_low <= dominance_pointwise_high):
        raise AssertionError("reversed dominance bootstrap interval")

    return {
        "estimate_difference": estimate_difference,
        "estimate_dominance": estimate_dominance,
        "boot_difference": boot_difference,
        "difference_pointwise_low": difference_pointwise_low,
        "difference_pointwise_high": difference_pointwise_high,
        "difference_simultaneous_low": difference_simultaneous_low,
        "difference_simultaneous_high": difference_simultaneous_high,
        "difference_critical": difference_critical,
        "dominance_pointwise_low": dominance_pointwise_low,
        "dominance_pointwise_high": dominance_pointwise_high,
        "dominance_simultaneous_low": dominance_simultaneous_low,
        "dominance_simultaneous_high": dominance_simultaneous_high,
        "dominance_critical": dominance_critical,
        "block_count": len(members),
    }


def block_sign_flip_max_t_pvalue(
    difference: np.ndarray,
    coordinates: list[tuple[int, int]],
    block_size: int,
    rng: np.random.Generator,
) -> tuple[float, float, int, str]:
    """Two-sided functional max-T sign-flip test on block medians."""

    _, members = block_members(coordinates, block_size)
    block_values = np.asarray(
        [np.median(difference[index], axis=0) for index in members],
        dtype=float,
    )
    block_count, distance_count = block_values.shape
    observed_mean = np.mean(block_values, axis=0)
    observed_se = np.std(block_values, axis=0, ddof=1) / math.sqrt(block_count)
    valid = observed_se > np.finfo(float).eps
    observed_t = np.zeros(distance_count, dtype=float)
    np.divide(
        observed_mean,
        observed_se,
        out=observed_t,
        where=valid,
    )
    observed_max_t = float(np.max(np.abs(observed_t)))

    if block_count <= 20:
        permutation_count = 1 << block_count
        identifiers = np.arange(permutation_count, dtype=np.uint64)[:, None]
        bits = (
            identifiers
            >> np.arange(block_count, dtype=np.uint64)[None, :]
        ) & 1
        signs = 1.0 - 2.0 * bits.astype(float)
        method = "exact"
    else:
        permutation_count = 100_000
        signs = rng.choice(
            np.asarray((-1.0, 1.0)),
            size=(permutation_count, block_count),
        )
        signs[0] = 1.0
        signs[1] = -1.0
        method = "monte_carlo"

    # Compute standard deviations explicitly in bounded chunks.  The algebraic
    # shortcut sum(x^2)-n*mean(x)^2 catastrophically cancels when all block
    # functions are almost equal, and can incorrectly return an exact p=0.
    chunk_size = max(64, min(4096, 2_000_000 // (block_count * distance_count)))
    exceedances = 0
    for start in range(0, permutation_count, chunk_size):
        stop = min(start + chunk_size, permutation_count)
        permuted_values = (
            signs[start:stop, :, None] * block_values[None, :, :]
        )
        permuted_mean = np.mean(permuted_values, axis=1)
        permuted_se = np.std(permuted_values, axis=1, ddof=1) / math.sqrt(
            block_count
        )
        permuted_t = np.zeros_like(permuted_mean)
        np.divide(
            permuted_mean,
            permuted_se,
            out=permuted_t,
            where=permuted_se > np.finfo(float).eps,
        )
        permuted_max_t = np.max(np.abs(permuted_t), axis=1)
        exceedances += int(
            np.count_nonzero(permuted_max_t >= observed_max_t - 1e-12)
        )
    if method == "exact":
        p_value = exceedances / permutation_count
        if p_value <= 0.0:
            raise AssertionError("an exact sign-flip p-value cannot be zero")
    else:
        p_value = (exceedances + 1) / (permutation_count + 1)
    return float(p_value), observed_max_t, block_count, method


def contiguous_segments(
    distance_nm: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
) -> list[tuple[float, float, str]]:
    sign = np.where(lower > 0.0, 1, np.where(upper < 0.0, -1, 0))
    segments: list[tuple[float, float, str]] = []
    start: int | None = None
    current = 0
    for index, value in enumerate(sign):
        if value != current:
            if current != 0 and start is not None:
                segments.append(
                    (
                        float(distance_nm[start]),
                        float(distance_nm[index - 1]),
                        "positive" if current > 0 else "negative",
                    )
                )
            start = index if value != 0 else None
            current = int(value)
    if current != 0 and start is not None:
        segments.append(
            (
                float(distance_nm[start]),
                float(distance_nm[-1]),
                "positive" if current > 0 else "negative",
            )
        )
    return segments


def segment_text(segments: list[tuple[float, float, str]]) -> str:
    if not segments:
        return "none"
    return "; ".join(
        f"{start:.0f}-{stop:.0f} nm ({direction})"
        for start, stop, direction in segments
    )


def plot_results(rows: list[dict]) -> None:
    colors = {(1, 2): "#1f78b4", (2, 4): "#e66101", (1, 4): "#238b45"}
    figure, axes = plt.subplots(
        3,
        2,
        figsize=(12.5, 11.0),
        sharex=True,
        constrained_layout=False,
    )
    for pair_index, pair in enumerate(COMPARISONS):
        selected = [
            row
            for row in rows
            if int(row["slow_speed_um_per_s"]) == pair[0]
            and int(row["fast_speed_um_per_s"]) == pair[1]
        ]
        selected.sort(key=lambda row: float(row["distance_nm"]))
        distance = np.asarray([float(row["distance_nm"]) for row in selected])
        difference = np.asarray(
            [float(row["median_difference_fast_minus_slow_pN"]) for row in selected]
        )
        point_low = np.asarray(
            [float(row["difference_pointwise_ci95_low_pN"]) for row in selected]
        )
        point_high = np.asarray(
            [float(row["difference_pointwise_ci95_high_pN"]) for row in selected]
        )
        simultaneous_low = np.asarray(
            [float(row["difference_simultaneous_ci95_low_pN"]) for row in selected]
        )
        simultaneous_high = np.asarray(
            [float(row["difference_simultaneous_ci95_high_pN"]) for row in selected]
        )
        dominance = np.asarray(
            [float(row["paired_dominance_fast_gt_slow"]) for row in selected]
        )
        dominance_low = np.asarray(
            [float(row["dominance_pointwise_ci95_low"]) for row in selected]
        )
        dominance_high = np.asarray(
            [float(row["dominance_pointwise_ci95_high"]) for row in selected]
        )
        dominance_sim_low = np.asarray(
            [float(row["dominance_simultaneous_ci95_low"]) for row in selected]
        )
        dominance_sim_high = np.asarray(
            [float(row["dominance_simultaneous_ci95_high"]) for row in selected]
        )

        color = colors[pair]
        force_axis = axes[pair_index, 0]
        probability_axis = axes[pair_index, 1]
        force_axis.fill_between(
            distance,
            simultaneous_low,
            simultaneous_high,
            color=color,
            alpha=0.15,
            label="95% simultaneous band",
        )
        force_axis.fill_between(
            distance,
            point_low,
            point_high,
            color=color,
            alpha=0.32,
            label="95% pointwise CI",
        )
        force_axis.plot(distance, difference, color=color, linewidth=2.0)
        force_axis.axhline(0.0, color="black", linewidth=0.8, linestyle="--")
        force_axis.set_ylabel(
            f"F({pair[1]}) - F({pair[0]})\n(pN)"
        )
        force_axis.grid(alpha=0.22)

        probability_axis.fill_between(
            distance,
            dominance_sim_low,
            dominance_sim_high,
            color=color,
            alpha=0.15,
            label="95% simultaneous band",
        )
        probability_axis.fill_between(
            distance,
            dominance_low,
            dominance_high,
            color=color,
            alpha=0.32,
            label="95% pointwise CI",
        )
        probability_axis.plot(distance, dominance, color=color, linewidth=2.0)
        probability_axis.axhline(0.5, color="black", linewidth=0.8, linestyle="--")
        probability_axis.set_ylim(-0.03, 1.03)
        probability_axis.set_ylabel(
            f"P[F({pair[1]}) > F({pair[0]})]"
        )
        probability_axis.grid(alpha=0.22)

    axes[0, 0].set_title("Paired median force difference")
    axes[0, 1].set_title("Same-location dominance fraction")
    axes[0, 0].legend(fontsize=8, loc="best")
    axes[0, 1].legend(fontsize=8, loc="best")
    for axis in axes[-1]:
        axis.set_xlabel("Sphere-plane separation, D (nm)")
    figure.suptitle(
        "0 wt% map-distribution separation: spatially paired conditional confidence",
        fontsize=14,
        y=0.985,
    )
    figure.text(
        0.5,
        0.012,
        (
            "4x4-pixel block bootstrap; simultaneous bands cover all three pairs "
            "and all 25-250 nm bins.  One sequential map per speed: speed is "
            "confounded with time/order."
        ),
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="0.30",
    )
    figure.subplots_adjust(
        left=0.105,
        right=0.985,
        bottom=0.075,
        top=0.935,
        hspace=0.18,
        wspace=0.20,
    )
    figure.savefig(FIGURE_STEM.with_suffix(".png"), dpi=300)
    figure.savefig(FIGURE_STEM.with_suffix(".svg"))
    plt.close(figure)


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)

    print("[1/5] decode 0 wt% sources and reconstruct line-corrected map forces")
    sources, maps = load_zero_wt_maps()
    common_coordinates = sorted(
        set.intersection(*(set(maps[speed].points) for speed in SPEEDS_UM_PER_S))
    )
    if len(common_coordinates) < 240:
        raise RuntimeError(
            f"only {len(common_coordinates)} complete spatial coordinates"
        )

    distance_mask = (
        (joint.BIN_CENTERS_NM >= DISTANCE_MIN_NM)
        & (joint.BIN_CENTERS_NM <= DISTANCE_MAX_NM)
    )
    distance_nm = joint.BIN_CENTERS_NM[distance_mask]
    forces: dict[int, np.ndarray] = {}
    for speed in SPEEDS_UM_PER_S:
        values = np.asarray(
            [
                maps[speed].points[coordinate].force_linear_corrected_pN[
                    distance_mask
                ]
                for coordinate in common_coordinates
            ],
            dtype=float,
        )
        forces[speed] = values

    complete = np.ones(len(common_coordinates), dtype=bool)
    for speed in SPEEDS_UM_PER_S:
        complete &= np.all(np.isfinite(forces[speed]), axis=1)
    if np.count_nonzero(complete) < len(common_coordinates):
        common_coordinates = [
            coordinate
            for coordinate, keep in zip(common_coordinates, complete)
            if keep
        ]
        forces = {speed: values[complete] for speed, values in forces.items()}
    paired_pixel_count = len(common_coordinates)
    if paired_pixel_count < 240:
        raise RuntimeError(
            f"only {paired_pixel_count} finite paired pixels across full D window"
        )
    if any(values.shape != (paired_pixel_count, distance_nm.size) for values in forces.values()):
        raise AssertionError("force matrix shape mismatch")

    print("[2/5] joint spatial-block bootstrap")
    bootstrap = bootstrap_functional_statistics(forces, common_coordinates, rng)

    print("[3/5] exact block sign-flip max-T tests and sensitivity")
    global_rows: list[dict] = []
    main_p_values: list[float] = []
    main_test_records: list[tuple[float, float, int, str]] = []
    sensitivity_records: list[dict[int, tuple[float, float, int, str]]] = []
    for slow, fast in COMPARISONS:
        difference = forces[fast] - forces[slow]
        sensitivity: dict[int, tuple[float, float, int, str]] = {}
        for block_size in (2, 4, 8):
            sensitivity[block_size] = block_sign_flip_max_t_pvalue(
                difference,
                common_coordinates,
                block_size,
                rng,
            )
        sensitivity_records.append(sensitivity)
        main_record = sensitivity[PRIMARY_BLOCK_SIZE_PIXELS]
        main_test_records.append(main_record)
        main_p_values.append(main_record[0])
    adjusted_p_values = holm_adjust(main_p_values)

    print("[4/5] result tables and figure")
    detail_rows: list[dict] = []
    selected_rows: list[dict] = []
    difference_estimate = np.asarray(bootstrap["estimate_difference"])
    dominance_estimate = np.asarray(bootstrap["estimate_dominance"])
    boot_difference = np.asarray(bootstrap["boot_difference"])
    for pair_index, (slow, fast) in enumerate(COMPARISONS):
        difference = forces[fast] - forces[slow]
        slow_q25, slow_median, slow_q75 = np.percentile(
            forces[slow],
            (25.0, 50.0, 75.0),
            axis=0,
        )
        fast_q25, fast_median, fast_q75 = np.percentile(
            forces[fast],
            (25.0, 50.0, 75.0),
            axis=0,
        )
        difference_mad = np.asarray(
            [base.robust_mad(difference[:, index]) for index in range(distance_nm.size)]
        )
        robust_effect = np.divide(
            difference_estimate[pair_index],
            difference_mad,
            out=np.full(distance_nm.shape, np.nan, dtype=float),
            where=difference_mad > np.finfo(float).eps,
        )
        sign_stability_positive = np.mean(
            boot_difference[:, pair_index] > 0.0,
            axis=0,
        )
        sign_stability_negative = np.mean(
            boot_difference[:, pair_index] < 0.0,
            axis=0,
        )
        for index, distance in enumerate(distance_nm):
            row = {
                "concentration_wt_percent": 0,
                "slow_speed_um_per_s": slow,
                "fast_speed_um_per_s": fast,
                "comparison": f"F_{fast}_minus_F_{slow}",
                "slow_map_timestamp": maps[slow].timestamp,
                "fast_map_timestamp": maps[fast].timestamp,
                "distance_nm": float(distance),
                "paired_pixel_count": paired_pixel_count,
                "spatial_block_size_pixels": PRIMARY_BLOCK_SIZE_PIXELS,
                "spatial_block_count": int(bootstrap["block_count"]),
                "slow_force_q25_pN": float(slow_q25[index]),
                "slow_force_median_pN": float(slow_median[index]),
                "slow_force_q75_pN": float(slow_q75[index]),
                "fast_force_q25_pN": float(fast_q25[index]),
                "fast_force_median_pN": float(fast_median[index]),
                "fast_force_q75_pN": float(fast_q75[index]),
                "median_difference_fast_minus_slow_pN": float(
                    difference_estimate[pair_index, index]
                ),
                "difference_mad_pN": float(difference_mad[index]),
                "robust_effect_median_over_mad": float(robust_effect[index]),
                "difference_pointwise_ci95_low_pN": float(
                    np.asarray(bootstrap["difference_pointwise_low"])[pair_index, index]
                ),
                "difference_pointwise_ci95_high_pN": float(
                    np.asarray(bootstrap["difference_pointwise_high"])[pair_index, index]
                ),
                "difference_simultaneous_ci95_low_pN": float(
                    np.asarray(bootstrap["difference_simultaneous_low"])[pair_index, index]
                ),
                "difference_simultaneous_ci95_high_pN": float(
                    np.asarray(bootstrap["difference_simultaneous_high"])[pair_index, index]
                ),
                "paired_dominance_fast_gt_slow": float(
                    dominance_estimate[pair_index, index]
                ),
                "dominance_pointwise_ci95_low": float(
                    np.asarray(bootstrap["dominance_pointwise_low"])[pair_index, index]
                ),
                "dominance_pointwise_ci95_high": float(
                    np.asarray(bootstrap["dominance_pointwise_high"])[pair_index, index]
                ),
                "dominance_simultaneous_ci95_low": float(
                    np.asarray(bootstrap["dominance_simultaneous_low"])[pair_index, index]
                ),
                "dominance_simultaneous_ci95_high": float(
                    np.asarray(bootstrap["dominance_simultaneous_high"])[pair_index, index]
                ),
                "bootstrap_median_positive_fraction": float(
                    sign_stability_positive[index]
                ),
                "bootstrap_median_negative_fraction": float(
                    sign_stability_negative[index]
                ),
                "simultaneous_force_difference_excludes_zero": bool(
                    np.asarray(bootstrap["difference_simultaneous_low"])[pair_index, index]
                    > 0.0
                    or np.asarray(bootstrap["difference_simultaneous_high"])[pair_index, index]
                    < 0.0
                ),
                "inference_scope": (
                    "conditional_within_area_spatial_blocks; one_sequential_map_"
                    "per_speed; speed_time_order_confounded"
                ),
            }
            detail_rows.append(row)
            if float(distance) in SELECTED_DISTANCES_NM:
                selected_rows.append(row.copy())

        simultaneous_low = np.asarray(bootstrap["difference_simultaneous_low"])[
            pair_index
        ]
        simultaneous_high = np.asarray(bootstrap["difference_simultaneous_high"])[
            pair_index
        ]
        segments = contiguous_segments(distance_nm, simultaneous_low, simultaneous_high)
        main_p, observed_max_t, main_blocks, main_method = main_test_records[pair_index]
        sensitivity = sensitivity_records[pair_index]
        global_rows.append(
            {
                "concentration_wt_percent": 0,
                "slow_speed_um_per_s": slow,
                "fast_speed_um_per_s": fast,
                "comparison": f"F_{fast}_minus_F_{slow}",
                "slow_map_timestamp": maps[slow].timestamp,
                "fast_map_timestamp": maps[fast].timestamp,
                "paired_pixel_count": paired_pixel_count,
                "primary_spatial_block_size_pixels": PRIMARY_BLOCK_SIZE_PIXELS,
                "primary_spatial_block_count": main_blocks,
                "functional_test_window_nm": f"{DISTANCE_MIN_NM:g}-{DISTANCE_MAX_NM:g}",
                "exact_block_sign_flip_maxT": observed_max_t,
                "exact_block_sign_flip_p_value": main_p,
                "holm_adjusted_p_value_across_three_pairs": adjusted_p_values[pair_index],
                "significant_at_0p05_conditional_spatial_level": bool(
                    adjusted_p_values[pair_index] < 0.05
                ),
                "fraction_distance_bins_simultaneous_ci_excludes_zero": float(
                    np.mean((simultaneous_low > 0.0) | (simultaneous_high < 0.0))
                ),
                "simultaneous_exclusion_segments": segment_text(segments),
                "block_2x2_count": sensitivity[2][2],
                "block_2x2_maxT_p_value": sensitivity[2][0],
                "block_2x2_test_method": sensitivity[2][3],
                "block_4x4_count": sensitivity[4][2],
                "block_4x4_maxT_p_value": sensitivity[4][0],
                "block_4x4_test_method": sensitivity[4][3],
                "block_8x8_count": sensitivity[8][2],
                "block_8x8_maxT_p_value": sensitivity[8][0],
                "block_8x8_test_method": sensitivity[8][3],
                "inference_scope": (
                    "conditional_within_area_spatial_blocks; not_independent_map_"
                    "replication_or_causal_velocity_confidence"
                ),
            }
        )

    write_csv(DETAIL_CSV, detail_rows)
    write_csv(SELECTED_CSV, selected_rows)
    write_csv(GLOBAL_CSV, global_rows)
    plot_results(detail_rows)

    print("[5/5] report and provenance")
    order_text = " -> ".join(
        f"{speed:g} um/s"
        for speed, _ in sorted(
            ((speed, maps[speed].timestamp) for speed in SPEEDS_UM_PER_S),
            key=lambda item: item[1],
        )
    )
    selected_lookup = {
        (
            int(row["slow_speed_um_per_s"]),
            int(row["fast_speed_um_per_s"]),
            float(row["distance_nm"]),
        ): row
        for row in selected_rows
    }
    report_lines = [
        "# 0 wt%三张map的分布分离置信度",
        "",
        "## 结论边界",
        "",
        (
            f"三张map的采集顺序为 `{order_text}`，每个速度只有一张map。"
            f"主分析使用三图共同有效的 `{paired_pixel_count}` 个物理点位和"
            f" `{int(bootstrap['block_count'])}` 个4x4空间block。"
        ),
        "",
        (
            "因此，下述95%区间是已测10x10 um区域内、以空间block为重采样单位的"
            "条件置信度；它不能解释为速度处理在重复实验中的95%置信度。速度与"
            "map时间/顺序完全混杂，真正的速度实验单位每组n=1。"
        ),
        "",
        "## 方法",
        "",
        (
            "在25-250 nm的每个5 nm距离bin，对相同物理pixel计算"
            " `F_fast - F_slow`。主效应是配对中位差；分布级效应是同点位"
            " `P(F_fast > F_slow)`。4x4 pixel block bootstrap共"
            f" {BOOTSTRAP_REPLICATES} 次。simultaneous band使用studentized max-T，"
            "一次覆盖三组比较和全部46个距离bin。"
        ),
        "",
        "## 全窗口检验",
        "",
        "| comparison | exact block max-T p | Holm p | simultaneous 95% band excludes 0 | 8x8-block sensitivity p |",
        "|:---|---:|---:|:---|---:|",
    ]
    for row in global_rows:
        report_lines.append(
            f"| {int(row['fast_speed_um_per_s'])}-{int(row['slow_speed_um_per_s'])} um/s "
            f"| {float(row['exact_block_sign_flip_p_value']):.5g} "
            f"| {float(row['holm_adjusted_p_value_across_three_pairs']):.5g} "
            f"| {row['simultaneous_exclusion_segments']} "
            f"| {float(row['block_8x8_maxT_p_value']):.5g} |"
        )
    report_lines.extend(
        [
            "",
            (
                "4x4结果量化区域内分离的一致性。8x8敏感性只有4个空间block，"
                "两侧exact sign-flip检验的最小可达p值受离散分辨率限制；它用于"
                "显示结论对更长空间相关尺度的依赖。"
            ),
            "",
            "## 选定距离",
            "",
            "数值为fast-minus-slow配对中位差；括号依次为pointwise 95% CI和全家族simultaneous 95% CI。",
            "",
        ]
    )
    for slow, fast in COMPARISONS:
        report_lines.extend(
            [
                f"### {fast}-{slow} um/s",
                "",
                "| D (nm) | median difference (pN) | pointwise 95% CI (pN) | simultaneous 95% CI (pN) | P(fast > slow) |",
                "|---:|---:|:---|:---|---:|",
            ]
        )
        for distance in SELECTED_DISTANCES_NM:
            row = selected_lookup[(slow, fast, distance)]
            report_lines.append(
                f"| {distance:.0f} "
                f"| {float(row['median_difference_fast_minus_slow_pN']):.3g} "
                f"| [{float(row['difference_pointwise_ci95_low_pN']):.3g}, "
                f"{float(row['difference_pointwise_ci95_high_pN']):.3g}] "
                f"| [{float(row['difference_simultaneous_ci95_low_pN']):.3g}, "
                f"{float(row['difference_simultaneous_ci95_high_pN']):.3g}] "
                f"| {float(row['paired_dominance_fast_gt_slow']):.3f} |"
            )
        report_lines.append("")
    report_lines.extend(
        [
            "## 独立force curves为何不进入主置信度",
            "",
            (
                "0 wt%独立force curves在1/2/4 um/s分别只有3/2/0条，既不平衡也"
                "不能与map pixel同点位配对。它们继续参与raw hard-contact sensitivity"
                "共识，但若加入速度分布检验会把不同层级的观测混为独立重复。"
            ),
            "",
            "## 文件",
            "",
            "- `zero_wt_distribution_separation.csv`: 全部距离的效应、pointwise和simultaneous区间。",
            "- `zero_wt_distribution_separation_selected.csv`: 六个选定距离。",
            "- `zero_wt_distribution_separation_global.csv`: functional max-T检验和block-size敏感性。",
            "- `figures/zero_wt_distribution_separation_confidence.png`: 置信带图。",
            "",
        ]
    )
    REPORT_MD.write_text("\n".join(report_lines), encoding="utf-8")

    provenance = {
        "analysis_script": str(Path(__file__).relative_to(ROOT)),
        "analysis_script_sha256": sha256_file(Path(__file__)),
        "dependency_script_sha256": {
            "analyze_velocity_joint_fit.py": sha256_file(
                ROOT / "analysis" / "analyze_velocity_joint_fit.py"
            ),
            "analyze_velocity_systematics.py": sha256_file(
                ROOT / "analysis" / "analyze_velocity_systematics.py"
            ),
            "fit_glycerol_surface_forces.py": sha256_file(
                ROOT / "analysis" / "fit_glycerol_surface_forces.py"
            ),
        },
        "raw_sources": [
            {
                "path": str(source.path.relative_to(ROOT)),
                "sha256": source.sha256,
                "source_type": source.source_type,
                "timestamp": source.timestamp.isoformat(),
            }
            for source in sources
        ],
        "python": platform.python_version(),
        "numpy": np.__version__,
        "random_seed": RANDOM_SEED,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "primary_spatial_block_size_pixels": PRIMARY_BLOCK_SIZE_PIXELS,
        "paired_pixel_count": paired_pixel_count,
        "distance_window_nm": [DISTANCE_MIN_NM, DISTANCE_MAX_NM],
        "distance_step_nm": float(np.diff(distance_nm)[0]),
        "force_unit": "pN",
        "force_branch": "per-curve far-linear-drift-corrected approach",
        "sensitivity_nm_per_V": maps[1].sensitivity_nm_per_V,
        "spring_constant_N_per_m": base.SPRING_CONSTANT_N_PER_M,
        "map_acquisition_order_um_per_s": [
            speed
            for speed, _ in sorted(
                ((speed, maps[speed].timestamp) for speed in SPEEDS_UM_PER_S),
                key=lambda item: item[1],
            )
        ],
        "experimental_unit_limit": (
            "one sequential map per speed; spatial blocks quantify conditional "
            "within-area heterogeneity, not replicated velocity-treatment uncertainty"
        ),
    }
    PROVENANCE_JSON.write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(f"paired pixels: {paired_pixel_count}")
    print(f"spatial blocks: {int(bootstrap['block_count'])}")
    print(f"difference simultaneous max-T critical value: {bootstrap['difference_critical']:.5g}")
    print(f"dominance simultaneous max-T critical value: {bootstrap['dominance_critical']:.5g}")
    for row in global_rows:
        print(
            f"{int(row['fast_speed_um_per_s'])}-{int(row['slow_speed_um_per_s'])} um/s: "
            f"p={float(row['exact_block_sign_flip_p_value']):.5g}, "
            f"Holm p={float(row['holm_adjusted_p_value_across_three_pairs']):.5g}, "
            f"segments={row['simultaneous_exclusion_segments']}"
        )
    print(f"wrote {DETAIL_CSV}")
    print(f"wrote {SELECTED_CSV}")
    print(f"wrote {GLOBAL_CSV}")
    print(f"wrote {REPORT_MD}")
    print(f"wrote {FIGURE_STEM.with_suffix('.png')}")
    print(f"wrote {FIGURE_STEM.with_suffix('.svg')}")
    print(f"wrote {PROVENANCE_JSON}")


if __name__ == "__main__":
    main()

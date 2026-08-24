#!/usr/bin/env python3
"""Plot every measured 0 wt% approach F-D curve in speed-highlight figures.

The reconstruction deliberately uses the same line-corrected operational
branch as the reported numerical force maps: each approach trace has its
initial 20% robust far-field line removed before voltage-to-deflection and
force conversion.  Each output highlights one speed while retaining the other
two as a light-gray reference.  The highlighted pointwise median and 25th--
75th percentile band are evaluated on a common separation grid covered by
every curve in that speed group; no distance extrapolation is performed.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

import fit_glycerol_surface_forces as base


CONCENTRATION = 0
MIN_DISTANCE_NM = 20.0
MAX_DISTANCE_NM = 900.0
ZOOM_MAX_DISTANCE_NM = 250.0
COLORS = {
    1: "#4292c6",  # blue
    2: "#f16913",  # orange
    4: "#41ab5d",  # green
}
MEDIAN_COLORS = {
    1: "#084594",
    2: "#a63603",
    4: "#006d2c",
}
BACKGROUND_COLOR = "#bdbdbd"
SUMMARY_GRID_STEP_NM = 1.0
OUTPUT_DIR = base.ROOT / "analysis" / "velocity_joint_fit_results" / "figures"
OUTPUT_STEM_TEMPLATE = "zero_wt_fd_highlight_{speed}_um_s"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def point_key(source: str, point_index: int | None) -> tuple[str, int | None]:
    return source, point_index


def parse_point_index(value: str) -> int | None:
    return None if value == "" else int(value)


def pointwise_summary(
    values: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return a pointwise IQR and median with every curve contributing.

    Exact duplicate distance samples, if present, are averaged before linear
    interpolation.  The grid is the intersection of the measured distance
    domains, so ``np.interp`` never extrapolates.
    """

    lower_bound = max(float(distance_nm[0]) for distance_nm, _ in values)
    upper_bound = min(float(distance_nm[-1]) for distance_nm, _ in values)
    grid_start = math.ceil(lower_bound / SUMMARY_GRID_STEP_NM) * SUMMARY_GRID_STEP_NM
    grid_stop = math.floor(upper_bound / SUMMARY_GRID_STEP_NM) * SUMMARY_GRID_STEP_NM
    if grid_stop <= grid_start:
        raise RuntimeError(
            f"no shared distance support: {lower_bound:.6g} to {upper_bound:.6g} nm"
        )
    grid_nm = np.arange(
        grid_start,
        grid_stop + 0.5 * SUMMARY_GRID_STEP_NM,
        SUMMARY_GRID_STEP_NM,
        dtype=float,
    )
    interpolated = np.empty((len(values), grid_nm.size), dtype=float)

    for row_index, (distance_nm, force_nN) in enumerate(values):
        unique_distance, inverse = np.unique(distance_nm, return_inverse=True)
        if unique_distance.size < 2:
            raise RuntimeError("curve has fewer than two unique distance samples")
        duplicate_counts = np.bincount(inverse)
        unique_force = np.bincount(inverse, weights=force_nN) / duplicate_counts
        interpolated[row_index] = np.interp(
            grid_nm,
            unique_distance,
            unique_force,
        )

    if not np.all(np.isfinite(interpolated)):
        raise RuntimeError("non-finite value in common-grid force matrix")
    q25_nN, median_nN, q75_nN = np.percentile(
        interpolated,
        (25.0, 50.0, 75.0),
        axis=0,
    )
    if not (
        np.all(q25_nN <= median_nN) and np.all(median_nN <= q75_nN)
    ):
        raise RuntimeError("pointwise quantile ordering failed")
    return grid_nm, q25_nN, median_nN, q75_nN


def main() -> None:
    event_rows = read_csv(
        base.ROOT
        / "analysis"
        / "velocity_systematics_results"
        / "curve_event_metrics.csv"
    )
    far_rows = read_csv(
        base.ROOT
        / "analysis"
        / "surface_force_results"
        / "far_field_drift_by_curve.csv"
    )
    sensitivity_rows = read_csv(
        base.ROOT
        / "analysis"
        / "surface_force_results"
        / "sensitivity_by_source.csv"
    )

    events: dict[tuple[str, int | None], dict[str, str]] = {}
    for row in event_rows:
        if int(row["concentration_wt_percent"]) != CONCENTRATION:
            continue
        key = point_key(row["source"], parse_point_index(row["point_index"]))
        events[key] = row

    far_fits: dict[tuple[str, int | None], tuple[float, float]] = {}
    for row in far_rows:
        if int(row["concentration_wt_percent"]) != CONCENTRATION:
            continue
        if row["far_field_fit_valid"] != "True":
            continue
        key = point_key(row["source"], parse_point_index(row["point_index"]))
        far_fits[key] = (
            float(row["slope_V_per_m"]),
            float(row["far_field_intercept_V"]),
        )

    sensitivities = {
        float(row["sensitivity_used_nm_per_V"]) * 1e-9
        for row in sensitivity_rows
        if int(row["concentration_wt_percent"]) == CONCENTRATION
    }
    if len(sensitivities) != 1:
        raise RuntimeError(f"expected one 0 wt% consensus sensitivity, got {sensitivities}")
    sensitivity_m_per_V = sensitivities.pop()

    curves: dict[int, list[tuple[np.ndarray, np.ndarray]]] = {
        speed: [] for speed in COLORS
    }
    source_counts: dict[int, dict[str, int]] = {
        speed: {"map": 0, "force": 0} for speed in COLORS
    }
    endpoint_aligned = 0

    data_dir = base.DATA_ROOT / str(CONCENTRATION)
    paths = sorted(data_dir.glob("*.jpk-force")) + sorted(
        data_dir.glob("*.jpk-force-map")
    )
    for path in paths:
        source = base.load_source(path, CONCENTRATION)
        source_name = str(path.relative_to(base.ROOT))
        for curve in source.curves:
            key = point_key(source_name, curve.point_index)
            if key not in events or key not in far_fits:
                raise RuntimeError(f"missing saved reconstruction metadata for {key}")
            event = events[key]
            speed = int(round(float(event["approach_speed_um_per_s"])))
            if speed not in COLORS:
                raise RuntimeError(f"unexpected speed {speed} um/s for {key}")

            slope_V_per_m, intercept_V = far_fits[key]
            baseline_V = slope_V_per_m * curve.measured_height_m + intercept_V
            corrected_V = curve.deflection_V - baseline_V
            delta_m = sensitivity_m_per_V * corrected_V

            contact_height_um = float(event["approach_contact_height_um"])
            if math.isfinite(contact_height_um):
                contact_height_m = contact_height_um * 1e-6
            else:
                # Match the primary source-preparation fallback for the single
                # 0 wt% pixel without a resolved rigid-contact interval.
                terminal = curve.measured_height_m[-20:] + delta_m[-20:]
                contact_height_m = float(np.quantile(terminal, 0.95))
                endpoint_aligned += 1

            distance_nm = (
                curve.measured_height_m + delta_m - contact_height_m
            ) * 1e9
            force_nN = base.SPRING_CONSTANT_N_PER_M * delta_m * 1e9
            usable = (
                np.isfinite(distance_nm)
                & np.isfinite(force_nN)
                & (distance_nm >= MIN_DISTANCE_NM)
                & (distance_nm <= MAX_DISTANCE_NM)
            )
            if np.count_nonzero(usable) < 20:
                raise RuntimeError(f"too little non-contact support for {key}")
            order = np.argsort(distance_nm[usable])
            curves[speed].append(
                (distance_nm[usable][order], force_nN[usable][order])
            )
            source_counts[speed][source.source_type] += 1

    total = sum(len(values) for values in curves.values())
    if total != len(events):
        raise RuntimeError(f"plotted {total} curves but expected {len(events)}")

    summaries = {
        speed: pointwise_summary(curves[speed]) for speed in (1, 2, 4)
    }

    # A shared full-data scale preserves direct comparability among the three
    # figures and ensures that no measured curve is silently clipped.
    all_force = np.concatenate(
        [force for values in curves.values() for _, force in values]
    )
    lower, upper = float(np.min(all_force)), float(np.max(all_force))
    span = upper - lower
    ylim = (lower - 0.02 * span, upper + 0.02 * span)

    plt.rcParams.update(
        {
            "font.size": 10.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    written_paths: list[Path] = []
    for highlighted_speed in (1, 2, 4):
        figure, axes = plt.subplots(
            1,
            2,
            figsize=(12.4, 5.2),
            gridspec_kw={"width_ratios": [1.15, 1.0]},
            constrained_layout=False,
        )

        # Draw all measurements from the two comparison speeds first so that
        # they remain visible but do not compete with the highlighted group.
        for speed in (1, 2, 4):
            if speed == highlighted_speed:
                continue
            for distance_nm, force_nN in curves[speed]:
                for axis in axes:
                    axis.plot(
                        distance_nm,
                        force_nN,
                        color=BACKGROUND_COLOR,
                        alpha=0.055,
                        linewidth=0.45,
                        solid_capstyle="round",
                        rasterized=True,
                        zorder=1,
                    )

        grid_nm, q25_nN, median_nN, q75_nN = summaries[highlighted_speed]
        for axis in axes:
            axis.fill_between(
                grid_nm,
                q25_nN,
                q75_nN,
                color=COLORS[highlighted_speed],
                alpha=0.25,
                linewidth=0,
                zorder=2,
            )

        for distance_nm, force_nN in curves[highlighted_speed]:
            for axis in axes:
                axis.plot(
                    distance_nm,
                    force_nN,
                    color=COLORS[highlighted_speed],
                    alpha=0.13,
                    linewidth=0.58,
                    solid_capstyle="round",
                    rasterized=True,
                    zorder=3,
                )

        for axis in axes:
            axis.plot(
                grid_nm,
                median_nN,
                color=MEDIAN_COLORS[highlighted_speed],
                linewidth=2.4,
                solid_capstyle="round",
                zorder=5,
            )
            axis.axhline(
                0.0,
                color="#666666",
                linewidth=0.8,
                linestyle=(0, (3, 3)),
                zorder=0,
            )
            axis.set_xlabel("Sphere-plane separation, D (nm)")
            axis.set_ylim(*ylim)
            axis.grid(color="#d9d9d9", linewidth=0.45, alpha=0.55)

        axes[0].set_xlim(MIN_DISTANCE_NM, MAX_DISTANCE_NM)
        axes[1].set_xlim(MIN_DISTANCE_NM, ZOOM_MAX_DISTANCE_NM)
        axes[0].set_title("Full non-contact range")
        axes[1].set_title("Surface-force range")
        axes[0].set_ylabel("Measured force, F (nN)")

        highlighted_counts = source_counts[highlighted_speed]
        background_count = total - len(curves[highlighted_speed])
        legend_handles = [
            Line2D(
                [0],
                [0],
                color=COLORS[highlighted_speed],
                alpha=0.65,
                linewidth=1.4,
                label=(
                    f"{highlighted_speed} µm/s individual curves "
                    f"(n={len(curves[highlighted_speed])}: "
                    f"{highlighted_counts['map']} map + "
                    f"{highlighted_counts['force']} independent)"
                ),
            ),
            Line2D(
                [0],
                [0],
                color=BACKGROUND_COLOR,
                linewidth=1.4,
                label=f"Other two speeds (n={background_count})",
            ),
            Patch(
                facecolor=COLORS[highlighted_speed],
                alpha=0.25,
                edgecolor="none",
                label="Pointwise 25th–75th percentile",
            ),
            Line2D(
                [0],
                [0],
                color=MEDIAN_COLORS[highlighted_speed],
                linewidth=2.4,
                label="Pointwise median",
            ),
        ]
        axes[0].legend(
            handles=legend_handles,
            loc="upper right",
            frameon=True,
            facecolor="white",
            edgecolor="#d0d0d0",
            framealpha=0.94,
            fontsize=8.8,
        )

        figure.suptitle(
            (
                "0 wt% glycerol-water: all approach F-D curves; "
                f"{highlighted_speed} µm/s highlighted"
            ),
            fontsize=14,
            y=0.975,
        )
        figure.text(
            0.5,
            0.015,
            (
                "Gray curves are all measurements at the other two speeds.  "
                "Median/IQR use only the common measured D-domain; no extrapolation.\n"
                "Per-curve initial-20% linear baseline; "
                f"S={sensitivity_m_per_V * 1e9:.4f} nm/V, "
                f"k={base.SPRING_CONSTANT_N_PER_M:.7f} N/m; "
                f"endpoint-aligned={endpoint_aligned}."
            ),
            ha="center",
            va="bottom",
            fontsize=8.2,
            color="#4d4d4d",
            linespacing=1.35,
        )
        figure.subplots_adjust(
            left=0.075,
            right=0.985,
            bottom=0.16,
            top=0.87,
            wspace=0.20,
        )

        output_stem = OUTPUT_DIR / OUTPUT_STEM_TEMPLATE.format(
            speed=highlighted_speed
        )
        png_path = output_stem.with_suffix(".png")
        svg_path = output_stem.with_suffix(".svg")
        figure.savefig(png_path, dpi=300)
        figure.savefig(svg_path)
        plt.close(figure)
        written_paths.extend((png_path, svg_path))

    for speed in (1, 2, 4):
        counts = source_counts[speed]
        print(
            f"{speed} um/s: {len(curves[speed])} curves "
            f"({counts['map']} map, {counts['force']} independent)"
        )
    print(f"endpoint-aligned curves: {endpoint_aligned}")
    print(f"plotted force range: {lower:.6g} to {upper:.6g} nN")
    for speed in (1, 2, 4):
        grid_nm = summaries[speed][0]
        print(
            f"{speed} um/s summary: {grid_nm[0]:.3f} to "
            f"{grid_nm[-1]:.3f} nm, step={SUMMARY_GRID_STEP_NM:.3f} nm, "
            f"n={len(curves[speed])} at every grid point"
        )
    for path in written_paths:
        print(f"wrote {path}")


if __name__ == "__main__":
    main()

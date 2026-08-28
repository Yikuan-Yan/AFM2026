#!/usr/bin/env python3
"""Plot one fixed 8x8 map pixel through the 27-08-26 chronology."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import analyze_27_08_26_palindrome_pilot as pilot  # noqa: E402
import fit_glycerol_surface_forces as base  # noqa: E402


DATA_ROOT = ROOT / "27-08-26" / "0"
RESULTS = ROOT / "analysis" / "palindrome_27_08_26_pilot_results"
ROW = 3
COLUMN = 3
TARGET_DISTANCES_NM = (20.0, 50.0, 100.0, 200.0)
COLORS = {0.05: "#2a9d8f", 0.1: "#e9c46a", 0.2: "#e76f51"}


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def curve_at_physical_pixel(source: base.SourceData) -> base.RawCurve:
    matches = [
        curve
        for curve in source.curves
        if curve.point_index is not None
        and base.map_pixel_from_index(
            curve.point_index,
            int(source.map_grid_i),
            source.map_back_and_forth,
        )
        == (ROW, COLUMN)
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"{source.path}: expected exactly one curve at {(ROW, COLUMN)}, got {len(matches)}"
        )
    return matches[0]


def reconstruct_curve(
    source: base.SourceData,
    sensitivity_m_per_V: float,
    spring_N_per_m: float,
) -> tuple[base.RawCurve, np.ndarray, np.ndarray, np.ndarray]:
    curve = curve_at_physical_pixel(source)
    far = base.fit_far_field_drift(curve.measured_height_m, curve.deflection_V)
    baseline = base.baseline_voltage(curve.measured_height_m, curve.deflection_V, far)
    corrected_V = curve.deflection_V - baseline
    contact = pilot.terminal_contact_fit(
        curve.measured_height_m,
        corrected_V,
        pilot.MAP_CONTACT_SPAN_NM,
    )
    slope, intercept, _, _ = base.robust_line(
        curve.measured_height_m[contact["start"] : contact["stop"]],
        corrected_V[contact["start"] : contact["stop"]],
    )
    contact_height_m = -intercept / slope
    delta_m = sensitivity_m_per_V * corrected_V
    distance_nm = (
        curve.measured_height_m + delta_m - contact_height_m
    ) * 1e9
    force_line_pN = spring_N_per_m * delta_m * 1e12
    reference_V = float(np.median(curve.deflection_V[: far.n_points]))
    force_constant_pN = (
        spring_N_per_m
        * sensitivity_m_per_V
        * (curve.deflection_V - reference_V)
        * 1e12
    )
    precontact = np.arange(distance_nm.size) < int(contact["start"])
    line_binned = pilot.bin_median(distance_nm, force_line_pN, precontact)
    constant_binned = pilot.bin_median(distance_nm, force_constant_pN, precontact)
    return curve, line_binned, constant_binned, distance_nm


def main() -> None:
    provenance = json.loads((RESULTS / "provenance.json").read_text(encoding="utf-8"))
    sensitivity_m_per_V = float(provenance["water_InvOLS_nm_per_V"]) * 1e-9
    spring_N_per_m = float(provenance["D4_spring_constant_N_per_m"])
    source_paths = sorted(
        path
        for path in DATA_ROOT.glob("*.jpk-force-map")
        if pilot.map_identity(path)[0] in (1, 2)
    )
    sources = [base.load_source(path.resolve(), 0) for path in source_paths]
    sources.sort(key=lambda source: source.timestamp)
    if len(sources) != 12:
        raise RuntimeError(f"expected 12 current pilot maps, got {len(sources)}")

    curve_rows: list[dict] = []
    slice_rows: list[dict] = []
    plot_records: list[dict] = []
    for order, source in enumerate(sources, start=1):
        block, speed = pilot.map_identity(source.path)
        curve, line, constant, raw_distance = reconstruct_curve(
            source, sensitivity_m_per_V, spring_N_per_m
        )
        if source.map_grid_i != 8 or source.map_grid_j != 8:
            raise RuntimeError(f"unexpected grid in {source.path}")
        physical_point_index = int(curve.point_index)
        reconstructed_pixel = base.map_pixel_from_index(
            physical_point_index, 8, source.map_back_and_forth
        )
        if reconstructed_pixel != (ROW, COLUMN):
            raise RuntimeError("physical pixel mapping changed during reconstruction")
        record = {
            "order": order,
            "block": block,
            "speed": speed,
            "timestamp": source.timestamp,
            "source": str(source.path.relative_to(ROOT)),
            "point_index": physical_point_index,
            "line": line,
            "constant": constant,
        }
        plot_records.append(record)
        for index, distance in enumerate(pilot.BIN_CENTERS_NM):
            if not np.isfinite(line[index]) and not np.isfinite(constant[index]):
                continue
            curve_rows.append(
                {
                    "source": record["source"],
                    "timestamp": source.timestamp.isoformat(),
                    "acquisition_order": order,
                    "block": block,
                    "nominal_speed_um_per_s": speed,
                    "physical_row_zero_based": ROW,
                    "physical_column_zero_based": COLUMN,
                    "point_index_in_archive": physical_point_index,
                    "distance_nm": distance,
                    "force_linear_drift_corrected_pN": line[index],
                    "force_far_constant_referenced_pN": constant[index],
                }
            )
        for distance in TARGET_DISTANCES_NM:
            index = int(np.flatnonzero(pilot.BIN_CENTERS_NM == distance)[0])
            slice_rows.append(
                {
                    "source": record["source"],
                    "timestamp": source.timestamp.isoformat(),
                    "acquisition_order": order,
                    "block": block,
                    "nominal_speed_um_per_s": speed,
                    "physical_row_zero_based": ROW,
                    "physical_column_zero_based": COLUMN,
                    "point_index_in_archive": physical_point_index,
                    "distance_nm": distance,
                    "force_linear_drift_corrected_pN": line[index],
                    "force_far_constant_referenced_pN": constant[index],
                }
            )

    write_csv(RESULTS / "fixed_pixel_row3_col3_force_curves.csv", curve_rows)
    write_csv(RESULTS / "fixed_pixel_row3_col3_force_slices.csv", slice_rows)

    # Consistent axes make chronology-to-chronology amplitude changes visible.
    visible_values = np.concatenate(
        [
            np.asarray(record["line"])[
                (pilot.BIN_CENTERS_NM >= 20) & (pilot.BIN_CENTERS_NM <= 200)
            ]
            for record in plot_records
        ]
    )
    visible_values = visible_values[np.isfinite(visible_values)]
    y_min = min(-25.0, float(np.quantile(visible_values, 0.005)))
    y_max = float(np.quantile(visible_values, 0.995) * 1.08)
    fig, axes = plt.subplots(3, 4, figsize=(15.0, 10.0), sharex=True, sharey=True)
    for ax, record in zip(axes.flat, plot_records, strict=True):
        mask = (
            (pilot.BIN_CENTERS_NM >= 20)
            & (pilot.BIN_CENTERS_NM <= 200)
            & np.isfinite(record["line"])
        )
        ax.plot(
            pilot.BIN_CENTERS_NM[mask],
            np.asarray(record["line"])[mask],
            color=COLORS[record["speed"]],
            lw=1.8,
        )
        ax.axhline(0.0, color="0.35", lw=0.6)
        ax.set_xlim(20, 200)
        ax.set_ylim(y_min, y_max)
        ax.set_title(
            f"#{record['order']}  B{record['block']}  {record['speed']:g} µm/s\n"
            f"{record['timestamp'].strftime('%H:%M:%S')}"
        )
        ax.grid(alpha=0.2)
    for ax in axes[-1]:
        ax.set_xlabel("Separation D (nm)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Force (pN)")
    fig.suptitle(
        "Fixed physical pixel (row 3, column 3): chronological line-corrected F–D"
    )
    fig.tight_layout()
    fig.savefig(RESULTS / "figures" / "fixed_pixel_row3_col3_chronological_FD.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.0), sharex=True)
    for ax, distance in zip(axes.flat, TARGET_DISTANCES_NM, strict=True):
        rows = sorted(
            [row for row in slice_rows if row["distance_nm"] == distance],
            key=lambda row: row["acquisition_order"],
        )
        order = np.asarray([row["acquisition_order"] for row in rows])
        force = np.asarray([row["force_linear_drift_corrected_pN"] for row in rows])
        ax.plot(order, force, color="0.4", lw=1.0)
        ax.scatter(
            order,
            force,
            c=[COLORS[row["nominal_speed_um_per_s"]] for row in rows],
            s=46,
            zorder=3,
        )
        ax.axvline(6.5, color="0.25", ls="--", lw=0.8)
        ax.set_title(f"D = {distance:g} nm")
        ax.set_ylabel("Force (pN)")
        ax.grid(alpha=0.2)
    labels = [
        f"B{record['block']} {record['speed']:g}" for record in plot_records
    ]
    for ax in axes[-1]:
        ax.set_xticks(np.arange(1, 13), labels, rotation=50, ha="right")
        ax.set_xlabel("Chronological order (block, speed µm/s)")
    handles = [
        plt.Line2D([], [], marker="o", ls="", color=color, label=f"{speed:g} µm/s")
        for speed, color in COLORS.items()
    ]
    fig.suptitle("Fixed physical pixel: force at selected separations", y=0.995)
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
        ncol=3,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(RESULTS / "figures" / "fixed_pixel_row3_col3_force_slices.png", dpi=220)
    plt.close(fig)

    pilot.create_manifest(
        [
            Path(__file__).resolve(),
            (RESULTS / "provenance.json").resolve(),
            (RESULTS / "fixed_pixel_row3_col3_force_curves.csv").resolve(),
            (RESULTS / "fixed_pixel_row3_col3_force_slices.csv").resolve(),
            (
                RESULTS
                / "figures"
                / "fixed_pixel_row3_col3_chronological_FD.png"
            ).resolve(),
            (
                RESULTS / "figures" / "fixed_pixel_row3_col3_force_slices.png"
            ).resolve(),
        ],
        RESULTS / "fixed_pixel_row3_col3_manifest.sha256",
    )

    # This supplemental script modifies files inside the pilot result tree.
    # Refresh the directory-level manifest as well so it cannot retain hashes
    # from an earlier fixed-pixel reconstruction.
    pilot_artifacts = [
        path
        for path in RESULTS.rglob("*")
        if path.is_file() and path.name != "artifact_manifest.sha256"
    ]
    pilot_artifacts.append(
        (ROOT / "analysis" / "analyze_27_08_26_palindrome_pilot.py").resolve()
    )
    pilot.create_manifest(
        pilot_artifacts,
        RESULTS / "artifact_manifest.sha256",
    )

    print(
        f"Wrote fixed pixel ({ROW}, {COLUMN}); point indices = "
        + ",".join(str(record["point_index"]) for record in plot_records)
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Plot one 8x8 map's force at a fixed separation versus acquisition order."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import theilslopes


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "analysis" / "palindrome_27_08_26_full_results"
PIXEL_SLICES = RESULTS / "pixel_force_slices_20_50_100_200nm.csv"
MAP_INVENTORY = RESULTS / "map_inventory_QC.csv"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--map-order", type=int, default=1)
    parser.add_argument("--distance-nm", type=float, default=50.0)
    args = parser.parse_args()

    inventory = read_csv(MAP_INVENTORY)
    matches = [
        row for row in inventory if int(row["acquisition_order"]) == args.map_order
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one map at acquisition order {args.map_order}, got {len(matches)}"
        )
    map_row = matches[0]
    source = map_row["source"]

    force_rows = [
        row
        for row in read_csv(PIXEL_SLICES)
        if row["source"] == source
        and np.isclose(
            float(row["distance_nm"]),
            args.distance_nm,
            rtol=0.0,
            atol=1e-9,
        )
    ]
    force_rows.sort(key=lambda row: int(row["point_index"]))
    point_indices = np.asarray(
        [int(row["point_index"]) for row in force_rows], dtype=np.int64
    )
    expected_indices = np.arange(64, dtype=np.int64)
    if len(force_rows) != 64 or not np.array_equal(point_indices, expected_indices):
        raise RuntimeError(
            "the selected map/distance must contain point indices 0 through 63 exactly once"
        )

    force_pN = np.asarray(
        [float(row["force_linear_drift_corrected_pN"]) for row in force_rows],
        dtype=np.float64,
    )
    if force_pN.shape != (64,) or not np.all(np.isfinite(force_pN)):
        raise RuntimeError("selected force vector is not 64 finite float64 values")

    sample_order = point_indices + 1
    slope_pN_per_sample, intercept_pN, _, _ = theilslopes(
        force_pN, sample_order, alpha=0.95
    )
    trend_pN = slope_pN_per_sample * sample_order + intercept_pN
    delta_pN = float(force_pN[-1] - force_pN[0])
    first_to_last_minutes = (
        float(map_row["map_protocol_duration_estimate_s"]) * 63.0 / 64.0 / 60.0
    )

    fig, ax = plt.subplots(figsize=(11.0, 6.2))
    for raster_row in range(8):
        if raster_row % 2 == 0:
            ax.axvspan(
                8 * raster_row + 0.5,
                8 * (raster_row + 1) + 0.5,
                color="#457b9d",
                alpha=0.045,
                linewidth=0,
            )
    for boundary in range(8, 64, 8):
        ax.axvline(boundary + 0.5, color="0.78", linewidth=0.8)

    ax.plot(sample_order, force_pN, color="#8ecae6", linewidth=1.2, zorder=1)
    ax.scatter(
        sample_order,
        force_pN,
        s=34,
        color="#1261a0",
        edgecolor="white",
        linewidth=0.45,
        zorder=2,
        label="Individual force curves",
    )
    ax.plot(
        sample_order,
        trend_pN,
        color="#e76f51",
        linestyle="--",
        linewidth=2.0,
        label="Theil–Sen visual guide",
    )
    ax.axhline(
        float(np.median(force_pN)),
        color="0.30",
        linestyle=":",
        linewidth=1.2,
        label="64-point median",
    )
    ax.scatter(
        [1, 64],
        [force_pN[0], force_pN[-1]],
        s=78,
        color=["#2a9d8f", "#d62828"],
        edgecolor="white",
        linewidth=0.7,
        zorder=4,
    )

    ax.annotate(
        f"1st: {force_pN[0]:.1f} pN",
        (1, force_pN[0]),
        xytext=(8, 10),
        textcoords="offset points",
        fontsize=9,
    )
    ax.annotate(
        f"64th: {force_pN[-1]:.1f} pN",
        (64, force_pN[-1]),
        xytext=(-8, -18),
        textcoords="offset points",
        ha="right",
        fontsize=9,
    )

    ax.set_xlim(0.5, 64.5)
    ax.set_xticks(np.arange(1, 65, 4))
    ax.set_xlabel("Acquisition order within the 8×8 map")
    ax.set_ylabel(f"Force at D = {args.distance_nm:g} nm (pN)")
    ax.grid(axis="y", alpha=0.22)
    ax.legend(frameon=False, ncol=3, loc="upper right")
    ax.set_title(
        f"Map #{args.map_order}: force versus acquisition order\n"
        f"block {int(map_row['block'])}, "
        f"U = {float(map_row['nominal_speed_um_per_s']):g} µm/s, "
        f"ΔF₆₄₋₁ = {delta_pN:+.1f} pN over ≈{first_to_last_minutes:.2f} min"
    )
    fig.text(
        0.5,
        0.01,
        "Vertical bands mark successive 8-point raster rows; the robust guide still mixes time and spatial position.",
        ha="center",
        fontsize=9,
        color="0.35",
    )
    fig.tight_layout(rect=(0, 0.035, 1, 1))

    distance_label = f"{args.distance_nm:g}".replace(".", "p")
    output = (
        RESULTS
        / "figures"
        / f"map_order{args.map_order:02d}_force_{distance_label}nm_vs_acquisition_order.png"
    )
    fig.savefig(output, dpi=240)
    plt.close(fig)

    manifest = output.with_suffix(".sha256")
    bound_paths = [Path(__file__).resolve(), PIXEL_SLICES, MAP_INVENTORY, output]
    manifest.write_text(
        "\n".join(
            f"{sha256_file(path)}  {path.relative_to(ROOT)}" for path in bound_paths
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Wrote {output}")
    print(
        f"F1={force_pN[0]:.6f} pN, F64={force_pN[-1]:.6f} pN, "
        f"delta={delta_pN:+.6f} pN, Theil-Sen endpoint change="
        f"{slope_pN_per_sample * 63.0:+.6f} pN"
    )


if __name__ == "__main__":
    main()

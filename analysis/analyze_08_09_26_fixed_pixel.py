#!/usr/bin/env python3
"""Reconstruct the 08-09-26 water maps with user-corrected cantilever-2 k.

Reuse the reference plot's raw decode, hard-contact consensus, force operator,
and 5 nm bins. Palindrome groups describe speed order, not liquid changes.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
import csv
import json
from pathlib import Path
import platform
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import scipy
from zipfile import ZipFile

import analyze_27_08_26_palindrome_pilot as pilot
import fit_glycerol_surface_forces as base
import plot_27_08_26_fixed_pixel_chronology as reference


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw" / "keeper_6f30d59535114f89b568"
DATA = RAW / "08-09-26"
OUT = ROOT / "analysis" / "water_08_09_26_results"
FIG = OUT / "figures"
K = 0.2736
K_UNCERTAINTY = 0.0063
ROW, COLUMN = 3, 3  # Zero-based, as in the requested reference figure.
COLORS = {1.0: "#2a9d8f", 2.0: "#e9c46a", 4.0: "#e76f51"}
TARGETS = np.asarray([20.0, 50.0, 100.0, 200.0])
TARGET_INDICES = np.asarray([int(np.flatnonzero(pilot.BIN_CENTERS_NM == d)[0]) for d in TARGETS])


def timestamp(value: str) -> datetime:
    return datetime.strptime(value, "%Y-%m-%d %H:%M:%S.%f %z")


def read_derived_csv(path: Path) -> list[dict]:
    """Read this script's own typed CSV outputs for an explicit export retry."""
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for key, value in row.items():
            if value in ("True", "False"):
                row[key] = value == "True"
            elif re.fullmatch(r"-?\d+", value):
                row[key] = int(value)
            else:
                try:
                    row[key] = float(value)
                except ValueError:
                    pass
    return rows


def metadata(path: Path) -> dict:
    with ZipFile(path) as archive:
        header = base.parse_properties(archive.read("header.properties"))
        key = "force-scan-map.settings.force-settings."
        duration = float(header[key + "extend-scan-time"])
        nominal_speed = abs(float(header[key + "relative-z-start"]) - float(header[key + "relative-z-end"])) / duration * 1e6
        segments = []
        flags = []
        fixed_segment = None
        for name in archive.namelist():
            match = re.fullmatch(r"index/(\d+)/segments/(\d+)/segment-header.properties", name)
            if match is None:
                continue
            segment = base.parse_properties(archive.read(name))
            start = timestamp(segment["force-segment-header.time-stamp"])
            stop = start + timedelta(seconds=float(segment["force-segment-header.duration"]))
            segments.append((start, stop))
            if match.group(1) == "28" and match.group(2) == "0":
                fixed_segment = segment
            flags.extend(k for k, v in segment.items() if k.startswith("force-segment-header.force-scan-flags.") and v == "true" and any(word in k for word in ("aborted", "out-of-range", "limit-exceeded")))
        if fixed_segment is None:
            raise RuntimeError(f"Missing fixed approach segment: {path.name}")
        point = base.parse_properties(archive.read("index/28/header.properties"))
        g = "force-scan-map.position-pattern.grid."
        return {
            "source": path.relative_to(ROOT).as_posix(),
            "map_start_time": min(x[0] for x in segments).isoformat(),
            "map_end_time": max(x[1] for x in segments).isoformat(),
            "fixed_pixel_time": timestamp(fixed_segment["force-segment-header.time-stamp"]).isoformat(),
            "instrument_scan_number": int(header["force-scan-map.scan-number"]),
            "nominal_speed_um_per_s": float(nominal_speed),
            "grid_i": int(header[g + "ilength"]),
            "grid_j": int(header[g + "jlength"]),
            "field_u_um": float(header[g + "ulength"]) * 1e6,
            "field_v_um": float(header[g + "vlength"]) * 1e6,
            "grid_center_x_um": float(header[g + "xcenter"]) * 1e6,
            "grid_center_y_um": float(header[g + "ycenter"]) * 1e6,
            "grid_theta": float(header[g + "theta"]),
            "grid_reflect": header[g + "reflect"],
            "back_and_forth": header["force-scan-map.position-pattern.back-and-forth"],
            "fixed_pixel_x_um": float(point["force-scan-series.header.position.x"]) * 1e6,
            "fixed_pixel_y_um": float(point["force-scan-series.header.position.y"]) * 1e6,
            "instrument_failure_flag_count": len(flags),
            "segment_count": len(segments),
        }


def alternate_curve(source: base.SourceData, sensitivity: float, span_nm: float = 50.0) -> np.ndarray:
    curve = reference.curve_at_physical_pixel(source)
    corrected = curve.deflection_V - base.baseline_voltage(curve.measured_height_m, curve.deflection_V, curve.far_field_fit)
    fit = pilot.terminal_contact_fit(curve.measured_height_m, corrected, span_nm)
    contact_height = -fit["intercept_V"] / fit["slope_V_per_m"]
    displacement = sensitivity * corrected
    distance = (curve.measured_height_m + displacement - contact_height) * 1e9
    return pilot.bin_median(distance, K * displacement * 1e12, np.arange(distance.size) < fit["start"])


def measurement_times(records: list[dict], basis: str = "fixed_pixel_time") -> np.ndarray:
    """Matplotlib dates from real segment times, preserving unequal time gaps."""
    if basis == "map_midpoint_time":
        times = [datetime.fromisoformat(r["map_start_time"]) +
                 (datetime.fromisoformat(r["map_end_time"]) - datetime.fromisoformat(r["map_start_time"])) / 2
                 for r in records]
    else:
        times = [datetime.fromisoformat(r[basis]) for r in records]
    result = np.asarray(mdates.date2num(times), dtype=float)
    if np.any(np.diff(result) <= 0):
        raise RuntimeError("Plot acquisition times must be strictly increasing")
    return result


def interval_between_maps(records: list[dict], before: int) -> tuple[float, float]:
    """Actual idle interval after a one-based map number."""
    return (float(mdates.date2num(datetime.fromisoformat(records[before - 1]["map_end_time"]))),
            float(mdates.date2num(datetime.fromisoformat(records[before]["map_start_time"]))))


def axes_chronology(ax: plt.Axes, records: list[dict], basis: str = "fixed_pixel_time", show_blocks: bool = True) -> None:
    if show_blocks:
        for boundary in range(6, len(records), 6):
            begin, end = interval_between_maps(records, boundary)
            ax.axvline((begin + end) / 2, color="0.45", ls="--", lw=0.7)
    first = datetime.fromisoformat(records[0]["map_start_time"])
    last = datetime.fromisoformat(records[-1]["map_end_time"])
    ax.set_xlim(mdates.date2num(first - timedelta(seconds=90)), mdates.date2num(last + timedelta(seconds=90)))
    ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(0, 60, 10), tz=first.tzinfo))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=first.tzinfo))
    ax.tick_params(axis="x", labelsize=9, labelrotation=0)
    label = {"fixed_pixel_time": "Fixed-pixel approach start time", "map_midpoint_time": "Map acquisition midpoint time", "map_start_time": "Map start time"}[basis]
    offset = first.strftime("%z")
    ax.set_xlabel(f"{label} (UTC{offset[:3]}:{offset[3:]})")
    ax.grid(alpha=0.2)


def force_figure(records: list[dict], matrices: np.ndarray, name: str, title: str, comparison: np.ndarray | None = None,
                 time_basis: str = "fixed_pixel_time") -> None:
    fig, axes = plt.subplots(2, 2, figsize=(18, 9), sharex=True)
    times = measurement_times(records, time_basis)
    colors = [COLORS[r["nominal_speed_um_per_s"]] for r in records]
    for ax, j, distance in zip(axes.flat, TARGET_INDICES, TARGETS, strict=True):
        ax.plot(times, matrices[:, j], color="0.4", lw=1)
        ax.scatter(times, matrices[:, j], c=colors, s=37, zorder=3)
        if comparison is not None:
            ax.plot(times, comparison[:, j], color="#4b67a1", ls=":", marker="x", lw=1, markersize=4, label="Far-constant baseline")
        ax.set_title(f"D = {distance:g} nm", fontsize=14)
        ax.set_ylabel("Force (pN)")
        missing = int(np.count_nonzero(~np.isfinite(matrices[:, j])))
        if missing:
            ax.text(0.98, 0.04, f"{missing}/{len(records)} bins without precontact samples", transform=ax.transAxes, ha="right", fontsize=8, color="0.35")
            ax.plot(times[~np.isfinite(matrices[:, j])], np.full(missing, 0.025), transform=ax.get_xaxis_transform(),
                    ls="", marker="x", ms=6, color="0.35", clip_on=False)
        if time_basis == "fixed_pixel_time" and distance == 100:
            i = next(i for i, r in enumerate(records) if r["acquisition_order"] == 25)
            ax.annotate("#25 · 1 µm/s", xy=(times[i], matrices[i, j]), xytext=(10, -18), textcoords="offset points",
                        fontsize=9, arrowprops={"arrowstyle": "-", "color": "0.4", "lw": 0.8})
        axes_chronology(ax, records, time_basis)
    handles = [plt.Line2D([], [], marker="o", ls="", color=color, label=f"{speed:g} µm/s") for speed, color in COLORS.items()]
    if comparison is not None:
        handles += [plt.Line2D([], [], ls=":", marker="x", color="#4b67a1", label="Far-constant baseline")]
    fig.suptitle(title, fontsize=16, y=0.985)
    aggregation = "Within-curve 5 nm bin medians" if time_basis == "fixed_pixel_time" else "Median across available pixels after 5 nm binning"
    fig.text(0.5, 0.945, f"08 Sep 2026 · water · cantilever 2 · k = 0.2736 ± 0.0063 N/m · {aggregation}", ha="center", fontsize=10)
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.93), ncol=len(handles), frameon=False)
    fig.text(0.5, 0.015, "Dashed lines: speed-palindrome boundaries. Bottom × marks: empty distance bins. Force: far-linear baseline unless indicated.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.035, 1, 0.885))
    fig.savefig(FIG / name, dpi=200)
    if name == "fixed_pixel_row3_col3_force_slices.png":
        fig.savefig(FIG / name.replace(".png", ".svg"))
    plt.close(fig)


def diagnostic_figure(records: list[dict]) -> None:
    panels = [
        ("map_contact_invOLS_median_nm_per_V", "Hard-contact InvOLS (nm/V)"),
        ("far_slope_median_pN_per_100nm", "Far-field slope (pN / 100 nm)"),
        ("terminal_load_median_nN", "Terminal contact load (nN)"),
        ("retract_pull_off_force_median_nN", "Retract pull-off (nN)"),
        ("contact_height_median_um", "Contact scanner height (µm)"),
        ("idle_before_map_s", "Idle time before map (s)"),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(18, 11), sharex=True)
    for ax, (key, label) in zip(axes.flat, panels, strict=True):
        values = [r[key] for r in records]
        time_basis = "map_start_time" if key == "idle_before_map_s" else "map_midpoint_time"
        times = measurement_times(records, time_basis)
        ax.plot(times, values, color="0.4", lw=1)
        ax.scatter(times, values, c=[COLORS[r["nominal_speed_um_per_s"]] for r in records], s=30, zorder=3)
        ax.set_ylabel(label)
        axes_chronology(ax, records, time_basis)
    fig.suptitle("08 Sep 2026: independent map and timing diagnostics", fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(FIG / "map_history_diagnostics.png", dpi=180)
    plt.close(fig)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    download = json.loads((RAW / "download_manifest.json").read_text(encoding="utf-8"))
    paths = sorted(DATA.glob("*.jpk-force-map"))
    if len(paths) != 36 or download["total_files"] != 36:
        raise RuntimeError("Expected the complete 36-map share snapshot")
    expected_hashes = {r["local_path"]: r["sha256"] for r in download["files"]}
    sources = []
    records = []
    for order, path in enumerate(paths, 1):
        source = base.load_source(path.resolve(), 0)
        m = metadata(path)
        if source.sha256 != expected_hashes[path.relative_to(ROOT).as_posix()]:
            raise RuntimeError(f"Input hash mismatch: {path}")
        if source.map_grid_i != 8 or source.map_grid_j != 8 or len(source.curves) != 64 or source.skipped_curves:
            raise RuntimeError(f"Incomplete 8x8 approach map: {path}")
        if sorted(c.point_index for c in source.curves) != list(range(64)) or m["segment_count"] != 128 or m["instrument_failure_flag_count"]:
            raise RuntimeError(f"Invalid curve identity or instrument flag: {path}")
        # Filename time is generated during acquisition (usually after the first
        # approach), with variable software delay. Use segment times for gaps.
        filename_offset = (source.timestamp - datetime.fromisoformat(m["map_start_time"]).replace(tzinfo=None)).total_seconds()
        acquisition_duration = (datetime.fromisoformat(m["map_end_time"]) - datetime.fromisoformat(m["map_start_time"])).total_seconds()
        if not 0 <= filename_offset <= acquisition_duration:
            raise RuntimeError("Filename timestamp falls outside actual acquisition")
        m["filename_time_minus_first_segment_s"] = filename_offset
        m.update(acquisition_order=order, block=(order - 1) // 6 + 1,
                 stored_spring_constant_N_per_m=source.stored_spring_constant_N_per_m,
                 spring_constant_used_N_per_m=K, spring_constant_supplied_uncertainty_N_per_m=K_UNCERTAINTY,
                 stored_InvOLS_nm_per_V=source.stored_sensitivity_m_per_V * 1e9)
        m["idle_before_map_s"] = 0.0 if not records else (datetime.fromisoformat(m["map_start_time"]) - datetime.fromisoformat(records[-1]["map_end_time"])).total_seconds()
        if m["idle_before_map_s"] < 0:
            raise RuntimeError("Overlapping map acquisition times")
        records.append(m)
        sources.append(source)
        if order % 6 == 0:
            print(f"Decoded and checked {order}/36 maps", flush=True)
    geometry = ["grid_i", "grid_j", "field_u_um", "field_v_um", "grid_center_x_um", "grid_center_y_um", "grid_theta", "grid_reflect", "back_and_forth"]
    if any(any(r[key] != records[0][key] for key in geometry) for r in records):
        raise RuntimeError("Map geometry changes: cannot claim a fixed physical grid pixel")
    if any(np.ptp([r[f"fixed_pixel_{axis}_um"] for r in records]) > 0.001 for axis in ("x", "y")):
        raise RuntimeError("Fixed pixel moves by more than 1 nm in scanner coordinates")
    for start in range(0, 36, 6):
        speeds = [r["nominal_speed_um_per_s"] for r in records[start:start + 6]]
        if speeds != speeds[::-1] or set(speeds) != set(COLORS):
            raise RuntimeError(f"Unexpected speed palindrome at map {start + 1}")
    # The legacy routines require identity from filenames. For these unnamed maps,
    # supply identity from verified scan settings only within this Python process.
    lookup = {source.path: (r["block"], r["nominal_speed_um_per_s"]) for source, r in zip(sources, records, strict=True)}
    pilot.map_identity = lambda path: lookup[path.resolve()]
    contact_rows, sensitivity = pilot.set_map_contact_fits(sources)
    print(f"Global water InvOLS = {sensitivity * 1e9:.8f} nm/V; corrected k = {K} N/m", flush=True)
    # Record the final within-map MAD retention separately from the first QC gate.
    for source in sources:
        subset = [r for r in contact_rows if Path(ROOT / r["source"]) == source.path]
        valid = np.asarray([r["sensitivity_nm_per_V"] for r in subset if r["valid_for_consensus"]])
        center, tolerance = np.median(valid), max(3.5 * pilot.robust_mad(valid), 0.5)
        for row in subset:
            row["retained_for_global_consensus"] = bool(row["valid_for_consensus"] and abs(row["sensitivity_nm_per_V"] - center) <= tolerance)
    retained = [r for r in contact_rows if r["retained_for_global_consensus"]]
    if not np.isclose(np.median([r["sensitivity_nm_per_V"] for r in retained]) * 1e-9, sensitivity, rtol=1e-12):
        raise RuntimeError("Contact-consensus identity failed")
    alternatives = {span: np.median([r[f"sensitivity_{span}nm_nm_per_V"] for r in retained]) * 1e-9 for span in (40, 60)}
    reused_reconstruction = "--reuse-reconstruction" in sys.argv
    if reused_reconstruction:
        map_rows = read_derived_csv(OUT / "map_inventory_QC.csv")
        pixel_rows = read_derived_csv(OUT / "pixel_QC.csv")
        map_force_rows = read_derived_csv(OUT / "map_force_curves.csv")
        if len(map_rows) != len(sources) or len(pixel_rows) != 36 * 64:
            raise RuntimeError("Incomplete saved reconstruction")
        for source, row in zip(sources, map_rows, strict=True):
            if (row["sha256"] != source.sha256 or row["spring_constant_used_N_per_m"] != K
                or row["spring_constant_supplied_uncertainty_N_per_m"] != K_UNCERTAINTY
                or not np.isclose(row["global_water_InvOLS_used_nm_per_V"] * 1e-9, sensitivity, rtol=1e-12)):
                raise RuntimeError("Saved inputs or calibration differ; run without --reuse-reconstruction")
        with np.load(OUT / "all_pixel_force_curves.npz", allow_pickle=False) as saved:
            if not np.array_equal(saved["distance_nm"], pilot.BIN_CENTERS_NM) or saved["line_pN"].shape != (36, 64, 100) or saved["constant_pN"].shape != (36, 64, 100):
                raise RuntimeError("Saved reconstruction axes differ")
            matrices = {str(s.path.relative_to(ROOT)) + "|" + method: saved[method + "_pN"][i].copy()
                        for i, s in enumerate(sources) for method in ("line", "constant")}
        print("Reused completed pixel reconstruction after input-hash/calibration checks; rechecking every fixed-pixel curve from raw data.", flush=True)
    else:
        print("Reconstructing all approach/retract pairs and force curves...", flush=True)
        pixel_rows, map_rows, map_force_rows, matrices = [], [], [], {}
        for i, source in enumerate(sources, 1):
            curve_part, map_part, force_part, matrix_part = pilot.map_analysis(K, [source])
            for row in curve_part + map_part + force_part:
                row["acquisition_order"] = i
            pixel_rows.extend(curve_part)
            map_rows.extend(map_part)
            map_force_rows.extend(force_part)
            matrices.update(matrix_part)
            if i % 6 == 0:
                print(f"Reconstructed {i}/36 maps, {i * 64} approach/retract pairs", flush=True)
    for m, r in zip(map_rows, records, strict=True):
        m.update(r)
    # Save the expensive reconstruction before selected-bin and figure checks.
    for filename, rows in {"water_contact_sensitivity_curves.csv": contact_rows, "map_inventory_QC.csv": map_rows,
                           "pixel_QC.csv": pixel_rows, "map_force_curves.csv": map_force_rows}.items():
        pilot.write_csv(OUT / filename, rows)
    np.savez_compressed(OUT / "all_pixel_force_curves.npz", distance_nm=pilot.BIN_CENTERS_NM,
                        acquisition_order=np.arange(1, 37), point_index=np.arange(64),
                        line_pN=np.stack([matrices[str(s.path.relative_to(ROOT)) + "|line"] for s in sources]),
                        constant_pN=np.stack([matrices[str(s.path.relative_to(ROOT)) + "|constant"] for s in sources]))
    fixed_rows, fixed_curve_rows, pixel_slice_rows = [], [], []
    fixed_line, fixed_constant, map_line, map_constant = [], [], [], []
    max_reference_error = 0.0
    for source, record in zip(sources, map_rows, strict=True):
        curve, line, constant, raw_distance = reference.reconstruct_curve(source, sensitivity, K)
        if curve.point_index != 28 or curve.contact_fit.r2 < pilot.MAP_CONTACT_R2_MIN:
            raise RuntimeError("Invalid fixed pixel or hard contact")
        key = str(source.path.relative_to(ROOT))
        matrix, constant_matrix = matrices[key + "|line"], matrices[key + "|constant"]
        if not np.allclose(line, matrix[28], equal_nan=True, rtol=1e-12, atol=1e-9):
            raise RuntimeError("Reference and all-pixel reconstructions disagree")
        max_reference_error = max(max_reference_error, float(np.nanmax(abs(line - matrix[28]))))
        fixed_line.append(line)
        fixed_constant.append(constant)
        map_line.append(np.nanmedian(matrix, axis=0))
        map_constant.append(np.nanmedian(constant_matrix, axis=0))
        local_line = alternate_curve(source, source.sensitivity_anchor_m_per_V)
        lines_by_span = {span: alternate_curve(source, alternatives[span], span) for span in (40, 60)}
        for j, distance in enumerate(pilot.BIN_CENTERS_NM):
            row = {"source": record["source"], "acquisition_order": record["acquisition_order"], "block": record["block"],
                   "map_start_time": record["map_start_time"], "fixed_pixel_time": record["fixed_pixel_time"],
                   "nominal_speed_um_per_s": record["nominal_speed_um_per_s"],
                   "physical_row_zero_based": ROW, "physical_column_zero_based": COLUMN, "point_index_in_archive": 28,
                   "distance_nm": float(distance), "force_linear_drift_corrected_pN": line[j],
                   "force_far_constant_referenced_pN": constant[j], "force_local_map_InvOLS_pN": local_line[j],
                   "force_contact_40nm_pN": lines_by_span[40][j], "force_contact_60nm_pN": lines_by_span[60][j],
                   "force_k_only_uncertainty_pN": abs(line[j]) * K_UNCERTAINTY / K,
                   "raw_samples_in_bin": int(np.count_nonzero((np.arange(raw_distance.size) < curve.contact_fit.start) & (raw_distance >= distance - 2.5) & (raw_distance < distance + 2.5)))}
            fixed_curve_rows.append(row)
            if distance in TARGETS:
                if np.isfinite(line[j]) != (row["raw_samples_in_bin"] > 0):
                    raise RuntimeError("Finite-bin value and raw sample count disagree")
                fixed_rows.append(row)
        for point in range(64):
            for j, distance in zip(TARGET_INDICES, TARGETS, strict=True):
                pr, pc = base.map_pixel_from_index(point, 8, source.map_back_and_forth)
                pixel_slice_rows.append({"acquisition_order": record["acquisition_order"], "block": record["block"], "nominal_speed_um_per_s": record["nominal_speed_um_per_s"], "point_index": point, "row": pr, "column": pc, "distance_nm": distance, "force_linear_drift_corrected_pN": matrix[point, j], "force_far_constant_referenced_pN": constant_matrix[point, j]})
    fixed_line, fixed_constant, map_line, map_constant = map(np.asarray, (fixed_line, fixed_constant, map_line, map_constant))
    # Adjacent differences are descriptive paired-pixel changes, not independent
    # replicate t tests and not proof that the operation was liquid replacement.
    changes = []
    for i in range(1, len(sources)):
        before, after = sources[i - 1], sources[i]
        for method in ("line", "constant"):
            difference = matrices[str(after.path.relative_to(ROOT)) + "|" + method] - matrices[str(before.path.relative_to(ROOT)) + "|" + method]
            for j, distance in zip(TARGET_INDICES, TARGETS, strict=True):
                values = difference[:, j]
                values = values[np.isfinite(values)]
                changes.append({"before_map": i, "after_map": i + 1, "distance_nm": distance, "baseline": method,
                                "before_speed_um_per_s": map_rows[i - 1]["nominal_speed_um_per_s"], "after_speed_um_per_s": map_rows[i]["nominal_speed_um_per_s"],
                                "idle_s": map_rows[i]["idle_before_map_s"], "paired_pixels": len(values),
                                "median_paired_change_pN": float(np.median(values)), "q25_change_pN": float(np.quantile(values, 0.25)),
                                "q75_change_pN": float(np.quantile(values, 0.75)), "positive_pixel_fraction": float(np.mean(values > 0)),
                                "fixed_pixel_change_pN": difference[28, j]})
    outputs = {"water_contact_sensitivity_curves.csv": contact_rows, "map_inventory_QC.csv": map_rows, "pixel_QC.csv": pixel_rows,
               "map_force_curves.csv": map_force_rows, "fixed_pixel_row3_col3_force_curves.csv": fixed_curve_rows,
               "fixed_pixel_row3_col3_force_slices.csv": fixed_rows, "pixel_force_slices.csv": pixel_slice_rows,
               "adjacent_map_changes.csv": changes}
    for filename, rows in outputs.items():
        pilot.write_csv(OUT / filename, rows)
    np.savez_compressed(OUT / "all_pixel_force_curves.npz", distance_nm=pilot.BIN_CENTERS_NM,
                        acquisition_order=np.arange(1, 37), point_index=np.arange(64),
                        line_pN=np.stack([matrices[str(s.path.relative_to(ROOT)) + "|line"] for s in sources]),
                        constant_pN=np.stack([matrices[str(s.path.relative_to(ROOT)) + "|constant"] for s in sources]))
    force_figure(map_rows, fixed_line, "fixed_pixel_row3_col3_force_slices.png", "Fixed pixel (3,3): force versus acquisition time")
    force_figure(map_rows, fixed_line, "fixed_pixel_row3_col3_baseline_sensitivity.png", "Fixed pixel (3,3): baseline sensitivity versus acquisition time", fixed_constant)
    force_figure(map_rows, map_line, "map_median_force_slices.png", "Map medians: force versus acquisition time", time_basis="map_midpoint_time")
    diagnostic_figure(map_rows)
    provenance = {
        "analysis": "08-09-26 water fixed physical pixel and map history diagnostics", "created_utc": datetime.now(timezone.utc).isoformat(),
        "sample_medium": "water", "sample_medium_authority": "user in current task", "cantilever": "cantilever2",
        "spring_constant_N_per_m": K, "supplied_k_uncertainty_N_per_m": K_UNCERTAINTY, "k_uncertainty_statistical_definition": "not specified by user; do not assign a confidence level",
        "spring_constant_authority": "user explicitly corrected wrong cantilever3 k in all file headers", "stored_wrong_k_N_per_m": sorted({s.stored_spring_constant_N_per_m for s in sources}),
        "header_k_to_corrected_k_force_ratio": K / sources[0].stored_spring_constant_N_per_m,
        "water_InvOLS_nm_per_V": sensitivity * 1e9, "water_InvOLS_method": "same reference global median after per-map hard-contact QC and MAD retention",
        "water_InvOLS_retained_contacts": len(retained), "total_contacts": len(contact_rows),
        "contact_span_nm": 50, "contact_sensitivity_checks_nm_per_V": {str(k): v * 1e9 for k, v in alternatives.items()},
        "contact_QC": {"r2_min": 0.995, "InvOLS_nm_per_V_range": [35, 90], "minimum_eligible_per_map": 56},
        "far_baseline": "per-curve robust line on initial max(80, ceil(0.2*N)) samples; subtract slope and intercept",
        "secondary_baseline": "subtract median voltage in same far window; separation coordinate still uses primary line baseline",
        "separation": "D_nm = (measured_height_m + InvOLS_m_per_V * V_corrected - contact_height_m) * 1e9",
        "force": "F_pN = user_k_N_per_m * InvOLS_m_per_V * V_corrected * 1e12",
        "binning": "median of precontact raw samples in [D-2.5,D+2.5) nm; no interpolation or extrapolation",
        "pixel": {"row_zero_based": ROW, "column_zero_based": COLUMN, "archive_index": 28, "x_um": map_rows[0]["fixed_pixel_x_um"], "y_um": map_rows[0]["fixed_pixel_y_um"]},
        "pixel_location_scope": "identical commanded scanner coordinates and map geometry; physical sample drift not measured",
        "block_definition": "six consecutive maps with verified palindromic nominal speeds; not inferred refresh events",
        "chronology": "map order sorted by timestamped filenames and checked using actual segment timestamps (+0200); idle time excludes actual acquisition duration",
        "figure_time_axes": {"fixed_pixel_figures": "actual point-28 approach segment start time, fixed_pixel_time; each approach lasts at most 1 s",
                             "map_aggregate_figures": "midpoint between actual map acquisition start and end",
                             "idle_time_diagnostic": "actual map start time", "timezone": "instrument UTC+02:00",
                             "palindrome_dividers": "midpoints of actual idle intervals between six-map groups; no operation-time claim"},
        "figure_aggregation": {"fixed_pixel_figures": "median only within each 5 nm bin of the single point-28 curve; no spatial averaging",
                               "map_median_figure": "median across available pixels of their individual 5 nm bin medians"},
        "validation": {"input_hashes_and_archive_CRC": "passed", "maps": 36, "approach_curves": sum(len(s.curves) for s in sources),
                       "reused_completed_reconstruction": reused_reconstruction,
                       "approach_parser_skips": sum(s.skipped_curves for s in sources), "retract_parser_skips": sum(r["retract_skipped_curves"] for r in map_rows),
                       "reference_fixed_curve_max_abs_difference_pN": max_reference_error, "fixed_slices_requested": len(fixed_rows),
                       "fixed_slices_observed": int(sum(np.isfinite(r["force_linear_drift_corrected_pN"]) for r in fixed_rows)),
                       "fixed_slices_missing": [{"map": r["acquisition_order"], "distance_nm": r["distance_nm"]} for r in fixed_rows if not np.isfinite(r["force_linear_drift_corrected_pN"])],
                       "fixed_slice_min_raw_samples": min(r["raw_samples_in_bin"] for r in fixed_rows)},
        "download": download, "software": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "code_sha256": {p.relative_to(ROOT).as_posix(): pilot.sha256_file(p) for p in [Path(__file__).resolve(), ROOT / "analysis" / "download_08_09_26_keeper.py", Path(base.__file__).resolve(), Path(pilot.__file__).resolve(), Path(reference.__file__).resolve(), ROOT / "analysis" / "analyze_velocity_systematics.py"]},
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    artifacts = [p for p in OUT.rglob("*") if p.is_file() and p.name != "artifact_manifest.sha256"]
    artifacts += [Path(__file__).resolve(), ROOT / "analysis" / "download_08_09_26_keeper.py"]
    pilot.create_manifest(artifacts, OUT / "artifact_manifest.sha256")
    print(f"Wrote {OUT}", flush=True)
    print(json.dumps(provenance["validation"], indent=2), flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Model-free time/history analysis of the 10-09-26 water force maps.

The experiment contains two distinct 18-map, 8x8, 2 um x 2 um locations.
Each location repeats the same three balanced six-map speed palindromes.  The
two locations are reconstructed with one independently calibrated D4 spring
constant and one global in-water hard-contact InvOLS, but all time and
same-pixel inferences are made within location.  The intervening survey maps
are retained in the inventory and excluded from the primary analysis.

Primary observables are measured force and endpoint-normalized force.  PB
parameters are intentionally left to a separate, secondary script.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import sys
from zipfile import ZipFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import analyze_08_09_26_fixed_pixel as chronology  # noqa: E402
import analyze_27_08_26_palindrome_pilot as pilot  # noqa: E402
import analyze_31_08_26_advanced as advanced  # noqa: E402
import fit_glycerol_surface_forces as base  # noqa: E402
import plot_27_08_26_fixed_pixel_chronology as reference  # noqa: E402


DATA = ROOT / "10-09-26"
RESULTS = ROOT / "analysis" / "time_dependence_10_09_26_results"
FIGURES = RESULTS / "figures"

TEMPERATURE_C = 25.6
PROBE_RADIUS_M = 4.546848945303745e-6
SPRING_CONSTANT_N_PER_M = 0.16482877909505173
SPRING_CONSTANT_REPEATABILITY_SD_N_PER_M = 0.01728282630942144
EMBEDDED_SPRING_CONSTANT_N_PER_M = 0.17713312641473047
EMBEDDED_INVOLS_NM_PER_V = 71.00506298341986
TARGET_DISTANCES_NM = np.asarray((20.0, 50.0, 100.0, 200.0))
TARGET_INDICES = np.asarray(
    [int(np.flatnonzero(np.isclose(pilot.BIN_CENTERS_NM, value))[0]) for value in TARGET_DISTANCES_NM]
)
SPEEDS_UM_PER_S = (1.0, 2.0, 4.0)
PALINDROME_BLOCKS = (
    (2.0, 1.0, 4.0, 4.0, 1.0, 2.0),
    (1.0, 4.0, 2.0, 2.0, 4.0, 1.0),
    (4.0, 2.0, 1.0, 1.0, 2.0, 4.0),
)
LOCATION_SCAN_RANGES = {"A": range(5, 23), "B": range(25, 43)}
COLORS = {1.0: "#2a9d8f", 2.0: "#e9c46a", 4.0: "#e76f51"}
RANDOM_SEED = 20260910
NORMALIZATION_MIN_DENOMINATOR_PN = 50.0


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
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict]:
    """Read a derived CSV with the repository's typed scalar conversion."""
    return chronology.read_derived_csv(path)


def finite(values) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    return array[np.isfinite(array)]


def quantile(values, q: float) -> float:
    selected = finite(values)
    return float(np.quantile(selected, q)) if selected.size else float("nan")


def timestamp(value: str) -> datetime:
    return datetime.fromisoformat(value)


def map_midpoint(record: dict) -> datetime:
    start = timestamp(record["map_start_time"])
    stop = timestamp(record["map_end_time"])
    return start + (stop - start) / 2


def raw_archive_metadata(path: Path) -> dict:
    """Inventory metadata only; it does not decode force-channel arrays."""
    with ZipFile(path) as archive:
        header = base.parse_properties(archive.read("header.properties"))
        key = "force-scan-map.settings.force-settings."
        grid = "force-scan-map.position-pattern.grid."
        duration = float(header[key + "extend-scan-time"])
        speed = abs(float(header[key + "relative-z-start"]) - float(header[key + "relative-z-end"])) / duration * 1e6
        segment_names = [
            name for name in archive.namelist()
            if name.startswith("index/") and name.endswith("segment-header.properties")
        ]
        crc_ok = archive.testzip() is None
        return {
            "source": path.relative_to(ROOT).as_posix(),
            "filename": path.name,
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "zip_crc_ok": crc_ok,
            "instrument_scan_number": int(header["force-scan-map.scan-number"]),
            "nominal_speed_um_per_s": speed,
            "grid_i": int(header[grid + "ilength"]),
            "grid_j": int(header[grid + "jlength"]),
            "field_u_um": float(header[grid + "ulength"]) * 1e6,
            "field_v_um": float(header[grid + "vlength"]) * 1e6,
            "grid_center_x_um": float(header[grid + "xcenter"]) * 1e6,
            "grid_center_y_um": float(header[grid + "ycenter"]) * 1e6,
            "grid_theta": float(header[grid + "theta"]),
            "grid_reflect": header[grid + "reflect"],
            "back_and_forth": header["force-scan-map.position-pattern.back-and-forth"],
            "segment_header_count": len(segment_names),
        }


def classify(path: Path, scan: int) -> str:
    if path.parent.name == "tests":
        return "test_excluded"
    if scan in LOCATION_SCAN_RANGES["A"]:
        return "location_A_primary"
    if scan in LOCATION_SCAN_RANGES["B"]:
        return "location_B_primary"
    if scan in (23, 24):
        return "survey_reposition_excluded"
    return "unclassified"


def build_inventory() -> tuple[list[dict], list[dict]]:
    paths = sorted(DATA.rglob("*.jpk-force-map"))
    if len(paths) != 42:
        raise RuntimeError(f"Expected 42 downloaded archives, found {len(paths)}")
    rows = []
    for path in paths:
        row = raw_archive_metadata(path)
        row["analysis_role"] = classify(path, row["instrument_scan_number"])
        rows.append(row)
    if not all(row["zip_crc_ok"] for row in rows):
        raise RuntimeError("At least one raw JPK archive failed ZIP CRC")
    scans = sorted(row["instrument_scan_number"] for row in rows)
    if scans != list(range(1, 43)):
        raise RuntimeError(f"Instrument scan sequence is not 1..42: {scans}")
    primary = [row for row in rows if row["analysis_role"].endswith("_primary")]
    if len(primary) != 36:
        raise RuntimeError("Expected two 18-map primary series")
    return rows, sorted(primary, key=lambda row: row["instrument_scan_number"])


def prepare_sources(primary_inventory: list[dict]) -> tuple[list[base.SourceData], list[dict]]:
    sources: list[base.SourceData] = []
    records: list[dict] = []
    for global_order, inventory in enumerate(primary_inventory, 1):
        path = ROOT / inventory["source"]
        source = base.load_source(path.resolve(), 0)
        if source.sha256 != inventory["sha256"]:
            raise RuntimeError(f"Raw input hash changed during load: {path}")
        if source.map_grid_i != 8 or source.map_grid_j != 8:
            raise RuntimeError(f"Primary map is not 8x8: {path}")
        point_indices = sorted(int(curve.point_index) for curve in source.curves if curve.point_index is not None)
        expected_indices = list(range(44)) if inventory["instrument_scan_number"] == 36 else list(range(64))
        if source.skipped_curves or point_indices != expected_indices:
            raise RuntimeError(
                f"Unexpected approach decode at scan {inventory['instrument_scan_number']}: "
                f"points={point_indices}, skipped={source.skipped_curves}"
            )
        metadata = chronology.metadata(path)
        scan = int(metadata["instrument_scan_number"])
        location = "A" if scan in LOCATION_SCAN_RANGES["A"] else "B"
        location_order = scan - (4 if location == "A" else 24)
        block = (location_order - 1) // 6 + 1
        metadata.update(
            {
                "location": location,
                "location_order": location_order,
                "block_in_location": block,
                "global_primary_order": global_order,
                "source": inventory["source"],
                "sha256": inventory["sha256"],
                "stored_spring_constant_N_per_m": source.stored_spring_constant_N_per_m,
                "stored_InvOLS_nm_per_V": source.stored_sensitivity_m_per_V * 1e9,
                "saved_approach_curves": len(source.curves),
                "missing_point_indices": ";".join(
                    str(point) for point in sorted(set(range(64)) - set(point_indices))
                ),
                "spring_constant_used_N_per_m": SPRING_CONSTANT_N_PER_M,
                "spring_constant_repeatability_sd_N_per_m": SPRING_CONSTANT_REPEATABILITY_SD_N_PER_M,
            }
        )
        metadata["map_midpoint_time"] = map_midpoint(metadata).isoformat()
        records.append(metadata)
        sources.append(source)

    for location in LOCATION_SCAN_RANGES:
        subset = [row for row in records if row["location"] == location]
        if len(subset) != 18 or [row["location_order"] for row in subset] != list(range(1, 19)):
            raise RuntimeError(f"Location {location} is not a complete 18-map sequence")
        geometry_fields = (
            "grid_i", "grid_j", "field_u_um", "field_v_um", "grid_center_x_um",
            "grid_center_y_um", "grid_theta", "grid_reflect", "back_and_forth",
        )
        for row in subset[1:]:
            for field in geometry_fields:
                if row[field] != subset[0][field]:
                    raise RuntimeError(f"Location {location} geometry changed at {field}")
        if any(
            np.ptp([row[f"fixed_pixel_{axis}_um"] for row in subset]) > 0.001
            for axis in ("x", "y")
        ):
            raise RuntimeError(f"Location {location} fixed pixel moved by more than 1 nm in scanner coordinates")
        for block_index, expected in enumerate(PALINDROME_BLOCKS, 1):
            observed = tuple(
                row["nominal_speed_um_per_s"]
                for row in subset[(block_index - 1) * 6 : block_index * 6]
            )
            if not np.allclose(observed, expected, rtol=0.0, atol=1e-9):
                raise RuntimeError(f"Location {location} block {block_index}: {observed} != {expected}")
        first_midpoint = map_midpoint(subset[0])
        for index, row in enumerate(subset):
            row["relative_time_min"] = (map_midpoint(row) - first_midpoint).total_seconds() / 60.0
            row["idle_before_map_s"] = (
                float("nan") if index == 0 else
                (timestamp(row["map_start_time"]) - timestamp(subset[index - 1]["map_end_time"])).total_seconds()
            )

    embedded_k = {round(source.stored_spring_constant_N_per_m, 14) for source in sources}
    embedded_s = {source.stored_sensitivity_m_per_V * 1e9 for source in sources}
    if (
        embedded_k != {round(EMBEDDED_SPRING_CONSTANT_N_PER_M, 14)}
        or len(embedded_s) != 1
        or not np.isclose(next(iter(embedded_s)), EMBEDDED_INVOLS_NM_PER_V, rtol=1e-12)
    ):
        raise RuntimeError(f"Unexpected embedded calibration fingerprint: k={embedded_k}, InvOLS={embedded_s}")
    return sources, records


def set_map_contact_fits_allowing_documented_partial_map(
    sources: list[base.SourceData], identity: dict[Path, tuple[int, float]]
) -> tuple[list[dict], float]:
    """Repository contact-consensus method with one explicit 44/64-map exception."""
    rows: list[dict] = []
    retained_all: list[float] = []
    for source in sources:
        block, speed = identity[source.path]
        valid_values: list[float] = []
        for curve in source.curves:
            curve.far_field_fit = base.fit_far_field_drift(
                curve.measured_height_m, curve.deflection_V
            )
            corrected = curve.deflection_V - base.baseline_voltage(
                curve.measured_height_m, curve.deflection_V, curve.far_field_fit
            )
            primary = pilot.terminal_contact_fit(
                curve.measured_height_m, corrected, pilot.MAP_CONTACT_SPAN_NM
            )
            checks = {
                span: pilot.terminal_contact_fit(curve.measured_height_m, corrected, span)
                for span in pilot.MAP_CONTACT_CHECKS_NM
            }
            curve.contact_fit = base.ContactFit(
                sensitivity_m_per_V=primary["sensitivity_m_per_V"],
                start=primary["start"],
                stop=primary["stop"],
                slope_V_per_m=primary["slope_V_per_m"],
                intercept_V=primary["intercept_V"],
                r2=primary["r2"],
            )
            valid = bool(
                primary["r2"] >= pilot.MAP_CONTACT_R2_MIN
                and 35e-9 <= primary["sensitivity_m_per_V"] <= 90e-9
            )
            if valid:
                valid_values.append(primary["sensitivity_m_per_V"])
            rows.append(
                {
                    "source": source.path.relative_to(ROOT).as_posix(),
                    "timestamp": source.timestamp.isoformat(),
                    "block": block,
                    "nominal_speed_um_per_s": speed,
                    "point_index": curve.point_index,
                    "sensitivity_nm_per_V": primary["sensitivity_m_per_V"] * 1e9,
                    "sensitivity_40nm_nm_per_V": checks[40.0]["sensitivity_m_per_V"] * 1e9,
                    "sensitivity_60nm_nm_per_V": checks[60.0]["sensitivity_m_per_V"] * 1e9,
                    "contact_fit_R2": primary["r2"],
                    "contact_points": primary["points"],
                    "valid_for_consensus": valid,
                }
            )
        values = np.asarray(valid_values, dtype=np.float64)
        minimum_valid = min(56, int(math.ceil(0.875 * len(source.curves))))
        if values.size < minimum_valid:
            raise RuntimeError(
                f"Too few valid contacts in {source.path}: {values.size}/{len(source.curves)}; "
                f"required {minimum_valid}"
            )
        center = float(np.median(values))
        mad = pilot.robust_mad(values)
        tolerance = max(3.5 * mad, 0.5e-9)
        retained = values[np.abs(values - center) <= tolerance]
        retained_all.extend(retained.tolist())
        source.sensitivity_anchor_m_per_V = float(np.median(retained))
        source.sensitivity_anchor_mad_m_per_V = pilot.robust_mad(retained)
        source.sensitivity_valid_curves = int(retained.size)
    consensus = float(np.median(np.asarray(retained_all)))
    for source in sources:
        source.sensitivity_used_m_per_V = consensus
        source.sensitivity_method = "all_primary_water_map_hard_contacts_global_median_one_documented_partial_map"
    return rows, consensus


def attach_calibration(
    sources: list[base.SourceData], records: list[dict]
) -> tuple[list[dict], float, dict[str, float]]:
    record_by_path = {Path(ROOT / row["source"]).resolve(): row for row in records}
    identity = {
        source.path: (
            (0 if record_by_path[source.path]["location"] == "A" else 3)
            + record_by_path[source.path]["block_in_location"],
            record_by_path[source.path]["nominal_speed_um_per_s"],
        )
        for source in sources
    }
    pilot.map_identity = lambda path: identity[path.resolve()]
    contact_rows, sensitivity = set_map_contact_fits_allowing_documented_partial_map(
        sources, identity
    )
    retained_values = []
    for source in sources:
        relative = source.path.relative_to(ROOT).as_posix()
        subset = [row for row in contact_rows if row["source"] == relative]
        valid = np.asarray(
            [row["sensitivity_nm_per_V"] for row in subset if row["valid_for_consensus"]],
            dtype=np.float64,
        )
        center = float(np.median(valid))
        tolerance = max(3.5 * pilot.robust_mad(valid), 0.5)
        for row in subset:
            keep = bool(
                row["valid_for_consensus"]
                and abs(float(row["sensitivity_nm_per_V"]) - center) <= tolerance
            )
            row["retained_for_global_consensus"] = keep
            row.update(
                {
                    "location": record_by_path[source.path]["location"],
                    "location_order": record_by_path[source.path]["location_order"],
                    "instrument_scan_number": record_by_path[source.path]["instrument_scan_number"],
                }
            )
            if keep:
                retained_values.append(float(row["sensitivity_nm_per_V"]))
    if not np.isclose(np.median(retained_values) * 1e-9, sensitivity, rtol=1e-12, atol=0.0):
        raise RuntimeError("Reconstructed global InvOLS consensus does not match pilot implementation")
    alternatives = {
        str(span): float(np.median([
            row[f"sensitivity_{int(span)}nm_nm_per_V"] for row in contact_rows
            if row["retained_for_global_consensus"]
        ]))
        for span in pilot.MAP_CONTACT_CHECKS_NM
    }
    return contact_rows, sensitivity, alternatives


def reconstruct(
    sources: list[base.SourceData], records: list[dict]
) -> tuple[list[dict], list[dict], list[dict], dict[str, np.ndarray], float]:
    curve_rows, map_rows, force_rows, matrices = pilot.map_analysis(
        SPRING_CONSTANT_N_PER_M, sources
    )
    record_by_source = {row["source"]: row for row in records}
    shared_fields = (
        "location", "location_order", "block_in_location", "global_primary_order",
        "instrument_scan_number", "map_start_time", "map_end_time", "map_midpoint_time",
        "fixed_pixel_time", "relative_time_min", "idle_before_map_s",
        "spring_constant_used_N_per_m", "spring_constant_repeatability_sd_N_per_m",
        "stored_spring_constant_N_per_m", "stored_InvOLS_nm_per_V",
    )
    for collection in (curve_rows, map_rows, force_rows):
        for row in collection:
            metadata = record_by_source[row["source"]]
            row.update({field: metadata[field] for field in shared_fields})
            row["acquisition_order"] = metadata["global_primary_order"]
            row["block"] = metadata["block_in_location"]

    # The reused pilot summary takes a finite median of candidate snap distances
    # even when the detector flag is false.  Preserve that diagnostic, then make
    # the public map statistic explicitly conditional on detected events, matching
    # pixel_QC and the PB guard-window definition.
    for map_row in map_rows:
        related = [row for row in curve_rows if row["source"] == map_row["source"]]
        detected = [
            float(row["approach_snap_distance_nm"])
            for row in related
            if row["approach_snap_detected"]
            and np.isfinite(float(row["approach_snap_distance_nm"]))
        ]
        map_row["approach_snap_candidate_unconditional_median_nm"] = map_row[
            "approach_snap_distance_median_nm"
        ]
        map_row["approach_snap_distance_median_nm"] = quantile(detected, 0.5)
        map_row["approach_snap_distance_max_nm"] = max(detected, default=float("nan"))
        map_row["approach_snap_detected_count"] = len(detected)

    maximum_reference_error = 0.0
    for source in sources:
        key = source.path.relative_to(ROOT).as_posix()
        curve, line, constant, _ = reference.reconstruct_curve(
            source, source.sensitivity_used_m_per_V, SPRING_CONSTANT_N_PER_M
        )
        if base.map_pixel_from_index(
            int(curve.point_index), int(source.map_grid_i), source.map_back_and_forth
        ) != (3, 3):
            raise RuntimeError("Reference fixed pixel identity changed")
        for method, independent in (("line", line), ("constant", constant)):
            saved = matrices[key + "|" + method][int(curve.point_index)]
            if not np.allclose(independent, saved, equal_nan=True, rtol=1e-12, atol=1e-9):
                raise RuntimeError(f"Independent fixed-pixel reconstruction mismatch: {key} {method}")
            finite_mask = np.isfinite(independent) & np.isfinite(saved)
            if np.any(finite_mask):
                maximum_reference_error = max(
                    maximum_reference_error,
                    float(np.max(np.abs(independent[finite_mask] - saved[finite_mask]))),
                )
    return curve_rows, map_rows, force_rows, matrices, maximum_reference_error


def selected_force_tables(
    sources: list[base.SourceData], records: list[dict], matrices: dict[str, np.ndarray]
) -> tuple[list[dict], list[dict]]:
    pixel_rows: list[dict] = []
    map_rows: list[dict] = []
    for source, record in zip(sources, records, strict=True):
        key = source.path.relative_to(ROOT).as_posix()
        for baseline, matrix in (("linear", matrices[key + "|line"]), ("constant", matrices[key + "|constant"])):
            for target_index, distance in zip(TARGET_INDICES, TARGET_DISTANCES_NM, strict=True):
                values = matrix[:, target_index]
                selected = finite(values)
                map_rows.append(
                    {
                        "source": key,
                        "location": record["location"],
                        "location_order": record["location_order"],
                        "block_in_location": record["block_in_location"],
                        "instrument_scan_number": record["instrument_scan_number"],
                        "map_midpoint_time": record["map_midpoint_time"],
                        "relative_time_min": record["relative_time_min"],
                        "nominal_speed_um_per_s": record["nominal_speed_um_per_s"],
                        "baseline": baseline,
                        "distance_nm": distance,
                        "force_median_pN": quantile(selected, 0.5),
                        "force_q25_pN": quantile(selected, 0.25),
                        "force_q75_pN": quantile(selected, 0.75),
                        "finite_pixels": selected.size,
                    }
                )
                for point, value in enumerate(values):
                    row, column = base.map_pixel_from_index(point, 8, source.map_back_and_forth)
                    pixel_rows.append(
                        {
                            "source": key,
                            "location": record["location"],
                            "location_order": record["location_order"],
                            "block_in_location": record["block_in_location"],
                            "instrument_scan_number": record["instrument_scan_number"],
                            "map_midpoint_time": record["map_midpoint_time"],
                            "relative_time_min": record["relative_time_min"],
                            "nominal_speed_um_per_s": record["nominal_speed_um_per_s"],
                            "baseline": baseline,
                            "distance_nm": distance,
                            "point_index": point,
                            "row": row,
                            "column": column,
                            "force_pN": float(value),
                        }
                    )
    return pixel_rows, map_rows


def exact_row_statistics(
    difference: np.ndarray, rows: np.ndarray, rng: np.random.Generator
) -> dict:
    valid = np.isfinite(difference)
    values = difference[valid]
    selected_rows = rows[valid]
    low, high = advanced.row_cluster_bootstrap_ci(values, selected_rows, rng)
    wilcoxon_p = float("nan")
    if values.size and np.any(values != 0.0):
        wilcoxon_p = float(stats.wilcoxon(values).pvalue)
    return {
        "paired_pixels": int(values.size),
        "difference_median": quantile(values, 0.5),
        "difference_q25": quantile(values, 0.25),
        "difference_q75": quantile(values, 0.75),
        "negative_pixel_fraction": float(np.mean(values < 0)) if values.size else float("nan"),
        "positive_pixel_fraction": float(np.mean(values > 0)) if values.size else float("nan"),
        "naive_pixel_wilcoxon_p": wilcoxon_p,
        "exact_physical_row_signflip_p": advanced.exact_row_signflip_p(values, selected_rows),
        "row_cluster_bootstrap_median_95CI_low": low,
        "row_cluster_bootstrap_median_95CI_high": high,
    }


def palindrome_analysis(
    sources: list[base.SourceData], records: list[dict], matrices: dict[str, np.ndarray], rng: np.random.Generator
) -> tuple[list[dict], list[dict]]:
    pair_rows: list[dict] = []
    aggregate_rows: list[dict] = []
    physical_rows = np.asarray(
        [base.map_pixel_from_index(point, 8, sources[0].map_back_and_forth)[0] for point in range(64)]
    )
    arrays = {
        (record["location"], record["location_order"], baseline): matrices[
            source.path.relative_to(ROOT).as_posix() + "|" + key
        ]
        for source, record in zip(sources, records, strict=True)
        for baseline, key in (("linear", "line"), ("constant", "constant"))
    }
    record_lookup = {(row["location"], row["location_order"]): row for row in records}
    for location in LOCATION_SCAN_RANGES:
        for block in range(1, 4):
            block_start = (block - 1) * 6 + 1
            for pair_depth, (early_offset, late_offset) in enumerate(((0, 5), (1, 4), (2, 3)), 1):
                early_order = block_start + early_offset
                late_order = block_start + late_offset
                early_record = record_lookup[(location, early_order)]
                late_record = record_lookup[(location, late_order)]
                speed = float(early_record["nominal_speed_um_per_s"])
                if not np.isclose(speed, late_record["nominal_speed_um_per_s"]):
                    raise RuntimeError("Palindrome pair speed identity failed")
                pair_center = map_midpoint(early_record) + (map_midpoint(late_record) - map_midpoint(early_record)) / 2
                separation_min = (map_midpoint(late_record) - map_midpoint(early_record)).total_seconds() / 60.0
                for baseline in ("linear", "constant"):
                    early = arrays[(location, early_order, baseline)]
                    late = arrays[(location, late_order, baseline)]
                    for target_index, distance in zip(TARGET_INDICES, TARGET_DISTANCES_NM, strict=True):
                        difference = late[:, target_index] - early[:, target_index]
                        row = {
                            "location": location,
                            "block_in_location": block,
                            "pair_depth": pair_depth,
                            "speed_um_per_s": speed,
                            "early_location_order": early_order,
                            "late_location_order": late_order,
                            "early_scan_number": early_record["instrument_scan_number"],
                            "late_scan_number": late_record["instrument_scan_number"],
                            "pair_center_time": pair_center.isoformat(),
                            "pair_center_relative_min": (
                                pair_center - map_midpoint(record_lookup[(location, 1)])
                            ).total_seconds() / 60.0,
                            "early_late_separation_min": separation_min,
                            "baseline": baseline,
                            "distance_nm": distance,
                        }
                        row.update(exact_row_statistics(difference, physical_rows, rng))
                        row["history_half_difference_median_pN"] = row["difference_median"] / 2.0
                        pair_rows.append(row)

        for baseline in ("linear", "constant"):
            for distance in TARGET_DISTANCES_NM:
                subset = [
                    row for row in pair_rows
                    if row["location"] == location and row["baseline"] == baseline
                    and row["distance_nm"] == distance
                ]
                medians = np.asarray([row["difference_median"] for row in subset])
                aggregate_rows.append(
                    {
                        "location": location,
                        "baseline": baseline,
                        "distance_nm": distance,
                        "palindrome_pairs": len(subset),
                        "pair_median_difference_median_pN": quantile(medians, 0.5),
                        "pair_median_difference_q25_pN": quantile(medians, 0.25),
                        "pair_median_difference_q75_pN": quantile(medians, 0.75),
                        "negative_pair_medians": int(np.sum(medians < 0)),
                        "positive_pair_medians": int(np.sum(medians > 0)),
                    }
                )
                for speed in SPEEDS_UM_PER_S:
                    speed_subset = [row for row in subset if row["speed_um_per_s"] == speed]
                    speed_medians = np.asarray([row["difference_median"] for row in speed_subset])
                    aggregate_rows.append(
                        {
                            "location": location,
                            "baseline": baseline,
                            "distance_nm": distance,
                            "speed_um_per_s": speed,
                            "palindrome_pairs": len(speed_subset),
                            "pair_median_difference_mean_pN": float(np.mean(speed_medians)),
                            "pair_median_difference_sd_pN": float(np.std(speed_medians, ddof=1)),
                            "pair_median_difference_values_pN": ";".join(f"{value:.9g}" for value in speed_medians),
                            "negative_pair_medians": int(np.sum(speed_medians < 0)),
                            "positive_pair_medians": int(np.sum(speed_medians > 0)),
                        }
                    )
    primary_p = [
        row["exact_physical_row_signflip_p"] for row in pair_rows if row["baseline"] == "linear"
    ]
    adjusted = advanced.benjamini_hochberg(primary_p)
    for row, q in zip(
        [row for row in pair_rows if row["baseline"] == "linear"], adjusted, strict=True
    ):
        row["exact_row_signflip_BH_q_across_primary_pair_tests"] = q
    return pair_rows, aggregate_rows


def same_speed_endpoint_analysis(
    sources: list[base.SourceData], records: list[dict], matrices: dict[str, np.ndarray], rng: np.random.Generator
) -> list[dict]:
    rows: list[dict] = []
    physical_rows = np.asarray(
        [base.map_pixel_from_index(point, 8, sources[0].map_back_and_forth)[0] for point in range(64)]
    )
    for location in LOCATION_SCAN_RANGES:
        subset = [(source, record) for source, record in zip(sources, records, strict=True) if record["location"] == location]
        for speed in SPEEDS_UM_PER_S:
            matches = [(source, record) for source, record in subset if np.isclose(record["nominal_speed_um_per_s"], speed)]
            early_source, early_record = matches[0]
            late_source, late_record = matches[-1]
            for baseline, key in (("linear", "line"), ("constant", "constant")):
                early = matrices[early_source.path.relative_to(ROOT).as_posix() + "|" + key]
                late = matrices[late_source.path.relative_to(ROOT).as_posix() + "|" + key]
                for target_index, distance in zip(TARGET_INDICES, TARGET_DISTANCES_NM, strict=True):
                    difference = late[:, target_index] - early[:, target_index]
                    row = {
                        "location": location,
                        "speed_um_per_s": speed,
                        "early_location_order": early_record["location_order"],
                        "late_location_order": late_record["location_order"],
                        "early_scan_number": early_record["instrument_scan_number"],
                        "late_scan_number": late_record["instrument_scan_number"],
                        "elapsed_min": (
                            map_midpoint(late_record) - map_midpoint(early_record)
                        ).total_seconds() / 60.0,
                        "baseline": baseline,
                        "distance_nm": distance,
                    }
                    row.update(exact_row_statistics(difference, physical_rows, rng))
                    rows.append(row)
    adjusted = advanced.benjamini_hochberg(
        [row["exact_physical_row_signflip_p"] for row in rows if row["baseline"] == "linear"]
    )
    for row, q in zip([row for row in rows if row["baseline"] == "linear"], adjusted, strict=True):
        row["exact_row_signflip_BH_q_across_primary_endpoint_tests"] = q
    return rows


def ols_hc3(x: np.ndarray, y: np.ndarray, contrast: np.ndarray) -> dict:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    contrast = np.asarray(contrast, dtype=np.float64)
    if x.ndim != 2 or y.shape != (x.shape[0],) or contrast.shape != (x.shape[1],):
        raise ValueError("OLS shape contract failed")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise FloatingPointError("OLS inputs contain NaN or Inf")
    rank = int(np.linalg.matrix_rank(x))
    if rank != x.shape[1]:
        raise RuntimeError("Time/speed design matrix is rank deficient")
    xtx_inverse = np.linalg.inv(x.T @ x)
    beta = xtx_inverse @ x.T @ y
    residual = y - x @ beta
    leverage = np.sum((x @ xtx_inverse) * x, axis=1)
    if np.any(leverage >= 1.0 - 1e-12):
        raise RuntimeError("HC3 leverage is singular")
    adjusted = residual / (1.0 - leverage)
    covariance = xtx_inverse @ (x.T @ ((adjusted**2)[:, None] * x)) @ xtx_inverse
    estimate = float(contrast @ beta)
    variance = float(contrast @ covariance @ contrast)
    se = math.sqrt(max(0.0, variance))
    degrees = x.shape[0] - x.shape[1]
    critical = float(stats.t.ppf(0.975, degrees))
    p = float(2 * stats.t.sf(abs(estimate / se), degrees)) if se > 0 else float("nan")
    sse = float(np.sum(residual**2))
    sst = float(np.sum((y - np.mean(y)) ** 2))
    residual_lag1_correlation = float("nan")
    if residual.size > 2 and np.std(residual[:-1]) > 0 and np.std(residual[1:]) > 0:
        residual_lag1_correlation = float(np.corrcoef(residual[:-1], residual[1:])[0, 1])

    def hac_summary(max_lag: int) -> dict:
        """Small-sample-corrected Newey-West sensitivity for ordered maps."""
        meat = x.T @ ((residual**2)[:, None] * x)
        for lag in range(1, max_lag + 1):
            weight = 1.0 - lag / (max_lag + 1.0)
            cross = np.zeros((x.shape[1], x.shape[1]), dtype=np.float64)
            for index in range(lag, x.shape[0]):
                cross += (
                    residual[index] * residual[index - lag]
                    * np.outer(x[index], x[index - lag])
                )
            meat += weight * (cross + cross.T)
        hac_covariance = (
            x.shape[0] / degrees * xtx_inverse @ meat @ xtx_inverse
        )
        hac_variance = float(contrast @ hac_covariance @ contrast)
        hac_se = math.sqrt(max(0.0, hac_variance))
        hac_p = (
            float(2 * stats.t.sf(abs(estimate / hac_se), degrees))
            if hac_se > 0 else float("nan")
        )
        return {
            f"hac{max_lag}_se": hac_se,
            f"hac{max_lag}_ci_low": estimate - critical * hac_se,
            f"hac{max_lag}_ci_high": estimate + critical * hac_se,
            f"hac{max_lag}_p": hac_p,
        }

    output = {
        "estimate": estimate,
        "hc3_se": se,
        "hc3_ci_low": estimate - critical * se,
        "hc3_ci_high": estimate + critical * se,
        "hc3_p": p,
        "n": x.shape[0],
        "parameters": x.shape[1],
        "rank": rank,
        "condition_number": float(np.linalg.cond(x)),
        "r2": 1.0 - sse / sst if sst > 0 else float("nan"),
        "rmse": math.sqrt(sse / x.shape[0]),
        "max_leverage": float(np.max(leverage)),
        "residual_lag1_correlation": residual_lag1_correlation,
        "beta": beta,
    }
    output.update(hac_summary(1))
    output.update(hac_summary(2))
    return output


def time_models(
    records: list[dict], map_qc: list[dict], map_slices: list[dict],
    normalized_summary: list[dict]
) -> list[dict]:
    rows: list[dict] = []
    for location in LOCATION_SCAN_RANGES:
        metadata = [row for row in records if row["location"] == location]
        times = np.asarray([row["relative_time_min"] for row in metadata], dtype=np.float64)
        time_fraction = (times - times[0]) / (times[-1] - times[0])
        speeds = np.asarray([
            next(
                row["gap_speed_20_200nm_median_um_per_s"] for row in map_qc
                if row["location"] == location and row["location_order"] == meta["location_order"]
            )
            for meta in metadata
        ], dtype=np.float64)
        speed_centered = speeds - np.mean(speeds)
        for baseline in ("linear", "constant"):
            for distance in TARGET_DISTANCES_NM:
                selected = sorted(
                    [row for row in map_slices if row["location"] == location and row["baseline"] == baseline and row["distance_nm"] == distance],
                    key=lambda row: row["location_order"],
                )
                values = np.asarray([row["force_median_pN"] for row in selected], dtype=np.float64)
                for model in ("linear_time_speed", "quadratic_time_speed"):
                    if model == "linear_time_speed":
                        design = np.column_stack([np.ones(18), time_fraction, speed_centered])
                        contrast = np.asarray([0.0, 1.0, 0.0])
                    else:
                        design = np.column_stack([np.ones(18), time_fraction, time_fraction**2, speed_centered])
                        contrast = np.asarray([0.0, 1.0, 1.0, 0.0])
                    fit = ols_hc3(design, values, contrast)
                    rows.append(
                        {
                            "location": location,
                            "metric": "force",
                            "baseline": baseline,
                            "distance_nm": distance,
                            "time_model": model,
                            "full_interval_time_change_at_mean_speed_pN": fit["estimate"],
                            "time_change_hc3_se_pN": fit["hc3_se"],
                            "time_change_hc3_95CI_low_pN": fit["hc3_ci_low"],
                            "time_change_hc3_95CI_high_pN": fit["hc3_ci_high"],
                            "time_change_hc3_p": fit["hc3_p"],
                            "time_change_hac1_se_pN": fit["hac1_se"],
                            "time_change_hac1_95CI_low_pN": fit["hac1_ci_low"],
                            "time_change_hac1_95CI_high_pN": fit["hac1_ci_high"],
                            "time_change_hac1_p": fit["hac1_p"],
                            "time_change_hac2_se_pN": fit["hac2_se"],
                            "time_change_hac2_95CI_low_pN": fit["hac2_ci_low"],
                            "time_change_hac2_95CI_high_pN": fit["hac2_ci_high"],
                            "time_change_hac2_p": fit["hac2_p"],
                            "actual_time_span_min": times[-1] - times[0],
                            "actual_gap_speed_coefficient_pN_per_um_s": float(fit["beta"][-1]),
                            "n_maps": fit["n"],
                            "design_rank": fit["rank"],
                            "design_parameters": fit["parameters"],
                            "design_condition_number": fit["condition_number"],
                            "r2": fit["r2"],
                            "rmse_pN": fit["rmse"],
                            "max_leverage": fit["max_leverage"],
                            "residual_lag1_correlation": fit["residual_lag1_correlation"],
                            "claim_status": "chronology_association_adjusted_for_measured_gap_speed",
                        }
                    )
        for distance in (50.0, 100.0):
            selected = sorted(
                [row for row in normalized_summary if row["location"] == location and row["distance_nm"] == distance],
                key=lambda row: row["location_order"],
            )
            values = np.asarray([row["normalized_force_median"] for row in selected], dtype=np.float64)
            design = np.column_stack([np.ones(18), time_fraction, speed_centered])
            fit = ols_hc3(design, values, np.asarray([0.0, 1.0, 0.0]))
            rows.append(
                {
                    "location": location,
                    "metric": "endpoint_normalized_force",
                    "baseline": "linear",
                    "distance_nm": distance,
                    "time_model": "linear_time_speed",
                    "full_interval_normalized_change_at_mean_speed": fit["estimate"],
                    "normalized_change_hc3_se": fit["hc3_se"],
                    "normalized_change_hc3_95CI_low": fit["hc3_ci_low"],
                    "normalized_change_hc3_95CI_high": fit["hc3_ci_high"],
                    "time_change_hc3_p": fit["hc3_p"],
                    "normalized_change_hac1_se": fit["hac1_se"],
                    "normalized_change_hac1_95CI_low": fit["hac1_ci_low"],
                    "normalized_change_hac1_95CI_high": fit["hac1_ci_high"],
                    "normalized_change_hac1_p": fit["hac1_p"],
                    "normalized_change_hac2_se": fit["hac2_se"],
                    "normalized_change_hac2_95CI_low": fit["hac2_ci_low"],
                    "normalized_change_hac2_95CI_high": fit["hac2_ci_high"],
                    "normalized_change_hac2_p": fit["hac2_p"],
                    "actual_time_span_min": times[-1] - times[0],
                    "actual_gap_speed_coefficient_per_um_s": float(fit["beta"][-1]),
                    "n_maps": fit["n"],
                    "design_rank": fit["rank"],
                    "design_parameters": fit["parameters"],
                    "design_condition_number": fit["condition_number"],
                    "r2": fit["r2"],
                    "rmse_pN": fit["rmse"],
                    "max_leverage": fit["max_leverage"],
                    "residual_lag1_correlation": fit["residual_lag1_correlation"],
                    "claim_status": "model_free_shape_chronology_adjusted_for_measured_gap_speed",
                }
            )
    return rows


def normalized_shape(
    sources: list[base.SourceData], records: list[dict], matrices: dict[str, np.ndarray], curve_rows: list[dict]
) -> tuple[list[dict], list[dict]]:
    noise_lookup = {
        (row["source"], int(row["point_index"])): float(row["far_noise_pN"])
        for row in curve_rows
    }
    selected_rows: list[dict] = []
    summary_rows: list[dict] = []
    near_index = int(np.flatnonzero(np.isclose(pilot.BIN_CENTERS_NM, 20.0))[0])
    far_index = int(np.flatnonzero(np.isclose(pilot.BIN_CENTERS_NM, 200.0))[0])
    for source, record in zip(sources, records, strict=True):
        key = source.path.relative_to(ROOT).as_posix()
        force = matrices[key + "|line"]
        denominator = force[:, near_index] - force[:, far_index]
        threshold = np.asarray([
            (
                max(3.0 * noise_lookup[(key, point)], NORMALIZATION_MIN_DENOMINATOR_PN)
                if (key, point) in noise_lookup
                else float("inf")
            )
            for point in range(64)
        ])
        valid = np.isfinite(denominator) & (denominator > threshold)
        normalized = np.full_like(force, np.nan)
        np.divide(
            force - force[:, far_index][:, None],
            denominator[:, None],
            out=normalized,
            where=valid[:, None],
        )
        for target_index, distance in zip(TARGET_INDICES, TARGET_DISTANCES_NM, strict=True):
            values = normalized[:, target_index]
            summary_rows.append(
                {
                    "source": key,
                    "location": record["location"],
                    "location_order": record["location_order"],
                    "relative_time_min": record["relative_time_min"],
                    "nominal_speed_um_per_s": record["nominal_speed_um_per_s"],
                    "distance_nm": distance,
                    "valid_pixels": int(np.count_nonzero(np.isfinite(values))),
                    "normalized_force_median": quantile(values, 0.5),
                    "normalized_force_q25": quantile(values, 0.25),
                    "normalized_force_q75": quantile(values, 0.75),
                }
            )
            for point, value in enumerate(values):
                row, column = base.map_pixel_from_index(point, 8, source.map_back_and_forth)
                selected_rows.append(
                    {
                        "source": key,
                        "location": record["location"],
                        "location_order": record["location_order"],
                        "relative_time_min": record["relative_time_min"],
                        "nominal_speed_um_per_s": record["nominal_speed_um_per_s"],
                        "point_index": point,
                        "row": row,
                        "column": column,
                        "distance_nm": distance,
                        "normalization_valid": bool(valid[point]),
                        "normalization_denominator_pN": float(denominator[point]),
                        "normalization_threshold_pN": float(threshold[point]),
                        "normalized_force": float(value),
                    }
                )
    return selected_rows, summary_rows


def qc_associations(map_rows: list[dict], map_slices: list[dict], normalized_summary: list[dict]) -> list[dict]:
    rows: list[dict] = []
    for location in LOCATION_SCAN_RANGES:
        metadata = sorted([row for row in map_rows if row["location"] == location], key=lambda row: row["location_order"])
        far_slope = np.asarray([row["far_slope_median_pN_per_100nm"] for row in metadata])
        for metric, source_rows, value_key in (
            (
                "absolute_force_50nm",
                sorted([row for row in map_slices if row["location"] == location and row["baseline"] == "linear" and row["distance_nm"] == 50.0], key=lambda row: row["location_order"]),
                "force_median_pN",
            ),
            (
                "endpoint_normalized_force_50nm",
                sorted([row for row in normalized_summary if row["location"] == location and row["distance_nm"] == 50.0], key=lambda row: row["location_order"]),
                "normalized_force_median",
            ),
        ):
            values = np.asarray([row[value_key] for row in source_rows])
            rho, p, n = advanced.stable_spearman(values, far_slope)
            rows.append(
                {
                    "location": location,
                    "metric": metric,
                    "qc_metric": "far_slope_median_pN_per_100nm",
                    "spearman_rho": rho,
                    "spearman_p": p,
                    "n_maps": n,
                    "interpretation": "association_not_causal_correction",
                }
            )
    for row, q in zip(rows, advanced.benjamini_hochberg([row["spearman_p"] for row in rows]), strict=True):
        row["BH_q_across_four_planned_associations"] = q
    return rows


def plot_chronology(records: list[dict], map_slices: list[dict]) -> None:
    figure, axes = plt.subplots(2, 4, figsize=(18, 8.5), sharex="row")
    for row_index, location in enumerate(("A", "B")):
        metadata = [row for row in records if row["location"] == location]
        times = np.asarray([row["relative_time_min"] for row in metadata])
        speeds = [row["nominal_speed_um_per_s"] for row in metadata]
        for column, distance in enumerate(TARGET_DISTANCES_NM):
            selected = sorted(
                [row for row in map_slices if row["location"] == location and row["baseline"] == "linear" and row["distance_nm"] == distance],
                key=lambda row: row["location_order"],
            )
            values = np.asarray([row["force_median_pN"] for row in selected])
            ax = axes[row_index, column]
            ax.plot(times, values, color="0.45", lw=1.0)
            ax.scatter(times, values, c=[COLORS[speed] for speed in speeds], s=38, zorder=3)
            for boundary in (6, 12):
                ax.axvline((times[boundary - 1] + times[boundary]) / 2, color="0.6", ls="--", lw=0.8)
            ax.axhline(0.0, color="0.7", lw=0.7)
            ax.grid(alpha=0.2)
            ax.set_title(f"D = {distance:g} nm")
            if column == 0:
                ax.set_ylabel(f"Location {location}\nMap-median force (pN)")
            if row_index == 1:
                ax.set_xlabel("Elapsed time within location (min)")
    handles = [plt.Line2D([], [], marker="o", ls="", color=COLORS[speed], label=f"{speed:g} µm/s") for speed in SPEEDS_UM_PER_S]
    figure.legend(handles=handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.965))
    figure.suptitle("10 Sep 2026: measured force chronology at two distinct locations", y=0.995)
    figure.text(0.5, 0.015, "Far-linear baseline; each point is the median of available 8×8 pixel-bin medians. Dashed lines delimit balanced palindrome blocks.", ha="center", fontsize=9)
    figure.tight_layout(rect=(0, 0.035, 1, 0.94))
    figure.savefig(FIGURES / "force_time_chronology.png", dpi=220)
    plt.close(figure)


def plot_palindrome(pair_rows: list[dict]) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(13.5, 8.5), sharex=False)
    for row_index, location in enumerate(("A", "B")):
        for column, distance in enumerate((20.0, 50.0)):
            ax = axes[row_index, column]
            selected = [row for row in pair_rows if row["location"] == location and row["baseline"] == "linear" and row["distance_nm"] == distance]
            for row in selected:
                y = row["difference_median"]
                low = y - row["difference_q25"]
                high = row["difference_q75"] - y
                ax.errorbar(row["pair_center_relative_min"], y, yerr=[[low], [high]], fmt="o", color=COLORS[row["speed_um_per_s"]], alpha=0.9, capsize=2)
            ax.axhline(0.0, color="0.45", lw=0.8)
            ax.grid(alpha=0.2)
            ax.set_title(f"Location {location}, D = {distance:g} nm")
            ax.set_xlabel("Pair-center elapsed time (min)")
            ax.set_ylabel("Late − early paired force (pN)")
    handles = [plt.Line2D([], [], marker="o", ls="", color=COLORS[speed], label=f"{speed:g} µm/s") for speed in SPEEDS_UM_PER_S]
    figure.legend(handles=handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.955))
    figure.suptitle("Same-speed palindrome history residuals", y=0.995)
    figure.tight_layout(rect=(0, 0, 1, 0.90))
    figure.savefig(FIGURES / "palindrome_history_residuals.png", dpi=220)
    plt.close(figure)


def plot_qc(records: list[dict], map_rows: list[dict]) -> None:
    metrics = (
        ("map_contact_invOLS_median_nm_per_V", "Local water InvOLS (nm/V)"),
        ("far_slope_median_pN_per_100nm", "Far slope (pN/100 nm)"),
        ("terminal_load_median_nN", "Terminal load (nN)"),
        ("approach_snap_distance_median_nm", "Detected snap distance (nm)"),
        ("retract_pull_off_force_median_nN", "Retract pull-off (nN)"),
        ("gap_speed_20_200nm_median_um_per_s", "Measured gap speed (µm/s)"),
    )
    figure, axes = plt.subplots(2, 3, figsize=(17, 8.5))
    for ax, (field, label) in zip(axes.flat, metrics, strict=True):
        for location, marker in (("A", "o"), ("B", "s")):
            subset = sorted([row for row in map_rows if row["location"] == location], key=lambda row: row["location_order"])
            times = [row["relative_time_min"] for row in subset]
            values = [row[field] for row in subset]
            ax.plot(times, values, marker=marker, ms=3.5, lw=1.2, label=f"Location {location}")
        ax.set_xlabel("Elapsed time within location (min)")
        ax.set_ylabel(label)
        ax.grid(alpha=0.2)
    axes[0, 0].legend(frameon=False)
    figure.suptitle("Calibration, baseline, contact and speed diagnostics")
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    figure.savefig(FIGURES / "qc_chronology.png", dpi=220)
    plt.close(figure)


def plot_normalized(normalized_summary: list[dict]) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(13.5, 8.5), sharex="row")
    for row_index, location in enumerate(("A", "B")):
        for column, distance in enumerate((50.0, 100.0)):
            selected = sorted([row for row in normalized_summary if row["location"] == location and row["distance_nm"] == distance], key=lambda row: row["location_order"])
            times = np.asarray([row["relative_time_min"] for row in selected])
            values = np.asarray([row["normalized_force_median"] for row in selected])
            q25 = np.asarray([row["normalized_force_q25"] for row in selected])
            q75 = np.asarray([row["normalized_force_q75"] for row in selected])
            speeds = [row["nominal_speed_um_per_s"] for row in selected]
            ax = axes[row_index, column]
            ax.fill_between(times, q25, q75, color="#457b9d", alpha=0.15)
            ax.plot(times, values, color="0.45", lw=1.0)
            ax.scatter(times, values, c=[COLORS[speed] for speed in speeds], s=38, zorder=3)
            ax.grid(alpha=0.2)
            ax.set_title(f"Location {location}, D = {distance:g} nm")
            ax.set_xlabel("Elapsed time within location (min)")
            ax.set_ylabel("Endpoint-normalized force")
    figure.suptitle("Model-free curve-shape chronology: (F(D)−F(200))/(F(20)−F(200))")
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(FIGURES / "normalized_shape_chronology.png", dpi=220)
    plt.close(figure)


def snap_safe_sensitivity_from_saved() -> tuple[list[dict], list[dict]]:
    """Re-evaluate headline summaries after a per-pixel five-nm snap guard."""
    records = read_csv(RESULTS / "map_inventory_QC.csv")
    pixels = read_csv(RESULTS / "pixel_QC.csv")
    normalized = read_csv(RESULTS / "normalized_force_pixels.csv")
    with np.load(RESULTS / "all_pixel_force_curves.npz", allow_pickle=False) as archive:
        distance = np.asarray(archive["distance_nm"], dtype=np.float64)
        force = np.asarray(archive["line_pN"], dtype=np.float64)
        orders = np.asarray(archive["global_primary_order"], dtype=np.int64)
    if force.shape != (36, 64, 100) or not np.array_equal(orders, np.arange(1, 37)):
        raise RuntimeError("Snap sensitivity received unexpected force-array axes")

    pixel_lookup = {
        (int(row["acquisition_order"]), int(row["point_index"])): row
        for row in pixels
    }
    normalized_lookup = {
        (row["source"], int(row["point_index"]), float(row["distance_nm"])): row
        for row in normalized
    }

    def is_safe(order: int, point: int, target_nm: float) -> bool:
        row = pixel_lookup.get((order, point))
        if row is None:
            return False
        if not row["approach_snap_detected"]:
            return True
        return bool(
            target_nm - pilot.BIN_WIDTH_NM / 2.0
            >= float(row["approach_snap_distance_nm"]) + 5.0
        )

    model_rows: list[dict] = []
    pair_rows: list[dict] = []
    for location in ("A", "B"):
        metadata = sorted(
            [row for row in records if row["location"] == location],
            key=lambda row: row["location_order"],
        )
        times = np.asarray([row["relative_time_min"] for row in metadata], dtype=np.float64)
        speeds = np.asarray(
            [row["gap_speed_20_200nm_median_um_per_s"] for row in metadata],
            dtype=np.float64,
        )
        design = np.column_stack(
            [np.ones(len(metadata)), times / times[-1], speeds - np.mean(speeds)]
        )
        for target in TARGET_DISTANCES_NM:
            column = int(np.flatnonzero(np.isclose(distance, target))[0])
            medians: list[float] = []
            counts: list[int] = []
            for row in metadata:
                order = int(row["acquisition_order"])
                valid = np.asarray(
                    [is_safe(order, point, target) for point in range(64)], dtype=bool
                )
                values = force[order - 1, :, column]
                selected = values[valid & np.isfinite(values)]
                if selected.size < 32:
                    raise RuntimeError("Snap-safe force median has fewer than 32 pixels")
                medians.append(float(np.median(selected)))
                counts.append(int(selected.size))
            fit = ols_hc3(design, np.asarray(medians), np.asarray([0.0, 1.0, 0.0]))
            model_rows.append(
                snap_model_row(location, "force", "pN", target, fit, counts)
            )

        for target in (50.0, 100.0):
            medians = []
            counts = []
            for row in metadata:
                order = int(row["acquisition_order"])
                values = []
                for point in range(64):
                    item = normalized_lookup.get((row["source"], point, target))
                    if item is None or not item["normalization_valid"]:
                        continue
                    value = float(item["normalized_force"])
                    if is_safe(order, point, 20.0) and np.isfinite(value):
                        values.append(value)
                if len(values) < 32:
                    raise RuntimeError("Snap-safe normalized median has fewer than 32 pixels")
                medians.append(float(np.median(values)))
                counts.append(len(values))
            fit = ols_hc3(design, np.asarray(medians), np.asarray([0.0, 1.0, 0.0]))
            model_rows.append(
                snap_model_row(
                    location, "endpoint_normalized_force", "dimensionless", target,
                    fit, counts,
                )
            )

        order_lookup = {int(row["location_order"]): row for row in metadata}
        for target in TARGET_DISTANCES_NM:
            column = int(np.flatnonzero(np.isclose(distance, target))[0])
            pair_medians: list[float] = []
            pair_counts: list[int] = []
            for block in range(3):
                start = block * 6 + 1
                for early_order, late_order in (
                    (start, start + 5), (start + 1, start + 4), (start + 2, start + 3)
                ):
                    early = int(order_lookup[early_order]["acquisition_order"])
                    late = int(order_lookup[late_order]["acquisition_order"])
                    valid = np.asarray(
                        [
                            is_safe(early, point, target)
                            and is_safe(late, point, target)
                            for point in range(64)
                        ],
                        dtype=bool,
                    )
                    difference = force[late - 1, :, column] - force[early - 1, :, column]
                    selected = difference[valid & np.isfinite(difference)]
                    if selected.size < 24:
                        raise RuntimeError("Snap-safe palindrome has fewer than 24 paired pixels")
                    pair_medians.append(float(np.median(selected)))
                    pair_counts.append(int(selected.size))
            pair_rows.append(
                {
                    "location": location,
                    "distance_nm": target,
                    "pairs": len(pair_medians),
                    "negative_pair_medians": int(np.sum(np.asarray(pair_medians) < 0)),
                    "pair_median_difference_median_pN": quantile(pair_medians, 0.5),
                    "minimum_snap_safe_paired_pixels": min(pair_counts),
                    "maximum_snap_safe_paired_pixels": max(pair_counts),
                    "guard": "bin_lower_edge_at_least_5nm_beyond_detected_snap_in_both_maps",
                }
            )
    return model_rows, pair_rows


def snap_model_row(
    location: str, metric: str, unit: str, distance_nm: float,
    fit: dict, counts: list[int],
) -> dict:
    return {
        "location": location,
        "metric": metric,
        "unit": unit,
        "distance_nm": distance_nm,
        "full_interval_change_at_mean_speed": fit["estimate"],
        "HC3_95CI_low": fit["hc3_ci_low"],
        "HC3_95CI_high": fit["hc3_ci_high"],
        "HC3_p": fit["hc3_p"],
        "HAC2_95CI_low": fit["hac2_ci_low"],
        "HAC2_95CI_high": fit["hac2_ci_high"],
        "HAC2_p": fit["hac2_p"],
        "residual_lag1_correlation": fit["residual_lag1_correlation"],
        "minimum_snap_safe_pixels_per_map": min(counts),
        "maximum_snap_safe_pixels_per_map": max(counts),
        "design_rank": fit["rank"],
        "design_parameters": fit["parameters"],
        "design_condition_number": fit["condition_number"],
        "guard": "bin_lower_edge_at_least_5nm_beyond_detected_snap",
    }


def create_manifest() -> None:
    files = [
        path for path in RESULTS.rglob("*")
        if path.is_file()
        and path.name != "artifact_manifest.sha256"
        and "apparent_pb" not in path.relative_to(RESULTS).parts
    ]
    lines = [f"{sha256_file(path)}  {path.relative_to(RESULTS).as_posix()}" for path in sorted(files)]
    (RESULTS / "artifact_manifest.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


def render_report(
    sensitivity_nm_per_v: float,
    alternatives: dict[str, float],
    records: list[dict],
    pair_aggregate: list[dict],
    endpoints: list[dict],
    models: list[dict],
    associations: list[dict],
    snap_models: list[dict],
    snap_pairs: list[dict],
    validation: dict,
) -> None:
    def aggregate(location: str, distance: float) -> dict:
        return next(
            row for row in pair_aggregate
            if row["location"] == location
            and row["baseline"] == "linear"
            and row["distance_nm"] == distance
            and ("speed_um_per_s" not in row or row["speed_um_per_s"] == "")
        )

    def linear_model(location: str, distance: float) -> dict:
        return next(row for row in models if row["location"] == location and row["metric"] == "force" and row["baseline"] == "linear" and row["distance_nm"] == distance and row["time_model"] == "linear_time_speed")

    def normalized_model(location: str, distance: float) -> dict:
        return next(
            row for row in models
            if row["location"] == location
            and row["metric"] == "endpoint_normalized_force"
            and row["distance_nm"] == distance
        )

    def snap_model(location: str, metric: str, distance: float) -> dict:
        return next(
            row for row in snap_models
            if row["location"] == location
            and row["metric"] == metric
            and row["distance_nm"] == distance
        )

    def snap_pair(location: str, distance: float) -> dict:
        return next(
            row for row in snap_pairs
            if row["location"] == location and row["distance_nm"] == distance
        )

    lines = [
        "# 10-09-26 pure-water AFM: two-location time/history analysis",
        "",
        "## Direct result",
        "",
        "The primary dataset consists of two separate 18-map balanced palindrome series. Time/history is estimated within each fixed scanner grid; the intervening survey/reposition maps are not treated as a continuation of the same physical pixels.",
        "",
        f"All primary forces use the independently air-calibrated D4 spring constant `k={SPRING_CONSTANT_N_PER_M:.9f} N/m` (repeatability SD `{SPRING_CONSTANT_REPEATABILITY_SD_N_PER_M:.9f} N/m`) and the current-map hard-contact consensus `InvOLS={sensitivity_nm_per_v:.6f} nm/V`. The uniform embedded `k={EMBEDDED_SPRING_CONSTANT_N_PER_M:.9f} N/m` and embedded `InvOLS={EMBEDDED_INVOLS_NM_PER_V:.6f} nm/V` are recorded but not used.",
        "",
        "| Location | span / min | D / nm | linear time+speed full-span change / pN | HC3 95% CI / pN | HAC(2) 95% CI / pN | 9 palindrome pair medians below zero | median late−early / pN |",
        "|---|---:|---:|---:|:---|:---|---:|---:|",
    ]
    for location in ("A", "B"):
        span = max(row["relative_time_min"] for row in records if row["location"] == location)
        for distance in TARGET_DISTANCES_NM:
            model = linear_model(location, distance)
            history = aggregate(location, distance)
            lines.append(
                f"| {location} | {span:.2f} | {distance:g} | {model['full_interval_time_change_at_mean_speed_pN']:+.2f} | [{model['time_change_hc3_95CI_low_pN']:+.2f}, {model['time_change_hc3_95CI_high_pN']:+.2f}] | [{model['time_change_hac2_95CI_low_pN']:+.2f}, {model['time_change_hac2_95CI_high_pN']:+.2f}] | {history['negative_pair_medians']}/9 | {history['pair_median_difference_median_pN']:+.2f} |"
            )
    lines += [
        "",
        "The regression is a chronological association adjusted for measured gap speed, not a kinetic-law or causal speed estimate. Palindrome late−early differences are same-speed contrasts, but pairs within a location share one evolving experimental history.",
        "",
        "### Per-pixel snap-safe sensitivity",
        "",
        "The table below removes a pixel at a target distance unless the full 5 nm bin begins at least 5 nm beyond its detected snap-in. This is a conservative sensitivity, not a replacement for the repository's pre-hard-contact binning.",
        "",
        "| Location | D / nm | snap-safe full-span change / pN | HC3 95% CI / pN | HAC(2) 95% CI / pN | safe pixels/map | negative palindrome pairs | snap-safe pair median / pN |",
        "|---|---:|---:|:---|:---|:---|---:|---:|",
    ]
    for location in ("A", "B"):
        for distance in TARGET_DISTANCES_NM:
            model = snap_model(location, "force", distance)
            history = snap_pair(location, distance)
            lines.append(
                f"| {location} | {distance:g} | "
                f"{model['full_interval_change_at_mean_speed']:+.2f} | "
                f"[{model['HC3_95CI_low']:+.2f}, {model['HC3_95CI_high']:+.2f}] | "
                f"[{model['HAC2_95CI_low']:+.2f}, {model['HAC2_95CI_high']:+.2f}] | "
                f"{model['minimum_snap_safe_pixels_per_map']}–"
                f"{model['maximum_snap_safe_pixels_per_map']} | "
                f"{history['negative_pair_medians']}/9 | "
                f"{history['pair_median_difference_median_pN']:+.2f} |"
            )
    lines += [
        "",
        "The 20 nm decrease remains under this guard. Location B's snap-safe 20 nm magnitude is smaller, so the unguarded value should not be interpreted as an all-pixel noncontact force. The 200 nm conclusion remains baseline-limited rather than snap-limited.",
        "The far-constant baseline branch retains negative HC3 intervals at A:20–100 nm and B:50–100 nm. A:200 nm and B:20 nm are therefore not promoted to baseline-robust effects; B:200 nm is inconclusive in both branches.",
        "",
        "![Force chronology](figures/force_time_chronology.png)",
        "",
        "![Palindrome history](figures/palindrome_history_residuals.png)",
        "",
        "## Experimental structure",
        "",
        "- Location A: instrument scans 5–22, fixed 8×8 grid, 2×2 µm, center (-0.967317, 3.962016) µm.",
        "- Scans 23–24 change field and center and are retained only as survey/reposition records.",
        "- Location B: instrument scans 25–42, fixed 8×8 grid, 2×2 µm, center (-3.176984, 3.524092) µm.",
        "- Each location uses [2,1,4,4,1,2], [1,4,2,2,4,1], [4,2,1,1,2,4] µm/s. Each speed is early/late paired once per block and occupies every palindrome depth once.",
        "- Location B scan 36 stopped after point indices 0–43: 44/64 approach/retract pairs are present and points 44–63 are explicitly missing. No values are imputed; all affected medians and pairs use their recorded finite-pixel count.",
        "- Commanded coordinates repeat exactly within each location, but sample drift is not independently measured. Locations A and B are not same-pixel replicates.",
        "",
        "## Reconstruction and calibration",
        "",
        "Each approach branch is independently corrected by a robust line fitted to its initial max(80,20%) samples. Hard contact uses a contiguous terminal 50 nm window. Separation and force are",
        "",
        "`D = h + InvOLS * V_corrected - h_contact`,  `F = k * InvOLS * V_corrected`.",
        "",
        "Absolute pN values and their regression intervals condition on the fixed k and InvOLS. The supplied k repeatability SD is not propagated into those intervals; it scales absolute force, while a common multiplicative scale cancels from endpoint normalization.",
        "",
        "Only raw samples preceding the terminal hard-contact fit window are binned in 5 nm intervals; empty bins remain missing. The separate snap-safe sensitivity removes bins too close to a detected jump-to-contact. A far-constant branch is preserved. The global hard-contact consensus uses map-level validity and MAD filtering, with 40/60 nm alternatives of " + f"`{alternatives['40.0']:.6f}/{alternatives['60.0']:.6f} nm/V`.",
        "",
        "![QC chronology](figures/qc_chronology.png)",
        "",
        "## Same-speed full-interval contrasts",
        "",
        "Each row compares the earliest and latest map acquired at the same nominal speed within one location. Pixels are paired by commanded grid coordinate; row-cluster intervals retain the eight physical rows as clusters.",
        "",
        "| Location | speed / µm/s | D / nm | elapsed / min | median paired change / pN | row-cluster 95% CI / pN | row sign-flip p |",
        "|---|---:|---:|---:|---:|:---|---:|",
    ]
    for row in endpoints:
        if row["baseline"] != "linear":
            continue
        lines.append(
            f"| {row['location']} | {row['speed_um_per_s']:g} | {row['distance_nm']:g} | {row['elapsed_min']:.2f} | {row['difference_median']:+.2f} | [{row['row_cluster_bootstrap_median_95CI_low']:+.2f}, {row['row_cluster_bootstrap_median_95CI_high']:+.2f}] | {row['exact_physical_row_signflip_p']:.5g} |"
        )
    lines += [
        "",
        "## Model-free curve shape and baseline coupling",
        "",
        "Endpoint normalization is `(F(D)-F(200))/(F(20)-F(200))`. A pixel is retained only when `F(20)-F(200) > max(3*far-noise, 50 pN)`. This removes constant offsets and separates curve-shape evolution from a uniform force scale, but it is not a Debye-length estimator.",
        "",
        "| Location | D / nm | full-span normalized change | HC3 95% CI | HAC(2) 95% CI | HC3 p | HAC(2) p | residual lag-1 r |",
        "|---|---:|---:|:---|:---|---:|---:|---:|",
    ]
    for location in ("A", "B"):
        for distance in (50.0, 100.0):
            model = normalized_model(location, distance)
            lines.append(
                f"| {location} | {distance:g} | {model['full_interval_normalized_change_at_mean_speed']:+.5f} | [{model['normalized_change_hc3_95CI_low']:+.5f}, {model['normalized_change_hc3_95CI_high']:+.5f}] | [{model['normalized_change_hac2_95CI_low']:+.5f}, {model['normalized_change_hac2_95CI_high']:+.5f}] | {model['time_change_hc3_p']:.4g} | {model['normalized_change_hac2_p']:.4g} | {model['residual_lag1_correlation']:+.3f} |"
            )
    lines += [
        "",
        "The same endpoint-normalized model after requiring the 20 nm denominator bin to pass the per-pixel snap guard is:",
        "",
        "| Location | D / nm | snap-safe normalized change | HC3 95% CI | HAC(2) 95% CI | safe pixels/map |",
        "|---|---:|---:|:---|:---|:---|",
    ]
    for location in ("A", "B"):
        for distance in (50.0, 100.0):
            model = snap_model(location, "endpoint_normalized_force", distance)
            lines.append(
                f"| {location} | {distance:g} | "
                f"{model['full_interval_change_at_mean_speed']:+.5f} | "
                f"[{model['HC3_95CI_low']:+.5f}, {model['HC3_95CI_high']:+.5f}] | "
                f"[{model['HAC2_95CI_low']:+.5f}, {model['HAC2_95CI_high']:+.5f}] | "
                f"{model['minimum_snap_safe_pixels_per_map']}–"
                f"{model['maximum_snap_safe_pixels_per_map']} |"
            )
    lines += [
        "",
        "HC3 is the primary heteroskedasticity-robust interval. Newey-West HAC(1/2) rows are finite-lag serial-correlation sensitivities; with only 18 maps they are not a fitted relaxation-noise model.",
        "",
        "![Normalized shape](figures/normalized_shape_chronology.png)",
        "",
        "The planned Spearman comparisons of absolute and normalized 50 nm force against far-field slope are:",
        "",
        "| Location | metric | rho | BH q |",
        "|---|---|---:|---:|",
    ]
    for row in associations:
        lines.append(f"| {row['location']} | {row['metric']} | {row['spearman_rho']:+.3f} | {row['BH_q_across_four_planned_associations']:.4g} |")
    lines += [
        "",
        "These correlations diagnose baseline coupling; they are not causal corrections.",
        "",
        "## Evidence boundary",
        "",
        "- The primary claim concerns measured finite-speed force/history and model-free normalized shape. It does not yet assign a unique mechanism.",
        "- Location, surface relaxation, solution/contamination evolution, contact history, baseline drift and optical sensitivity remain possible contributors.",
        "- The two locations are analyzed separately. Agreement in direction is replication of an association, not one continuous relaxation curve.",
        "- Pixel p-values are retained only as naive diagnostics. Physical-row sign flips and row-cluster bootstrap intervals are the conservative spatial summaries; the map/pair remains the experimental unit.",
        "- A separate PB analysis may report apparent/model-conditioned lambda_D and |psi| only after snap-safe fit-window checks.",
        "",
        "## Numerical and provenance checks",
        "",
        f"- Raw archives: `{validation['raw_archives']}`; ZIP CRC failures: `0`; primary maps: `{validation['primary_maps']}`; approach/retract pairs: `{validation['approach_curves']}`.",
        f"- Independent fixed-pixel reconstruction maximum absolute difference: `{validation['fixed_pixel_reconstruction_max_abs_difference_pN']:.3e} pN`.",
        f"- Force arrays: `{validation['force_array_shape']}`; finite target-bin fraction: `{validation['finite_target_bin_fraction']:.6f}`.",
        f"- All time/speed designs had full rank; maximum reported condition number: `{validation['maximum_time_design_condition_number']:.3f}`.",
        f"- All snap-safe sensitivity designs had full rank; the minimum retained per-map count was `{validation['minimum_snap_safe_pixels_per_map']}` pixels.",
        f"- Existing PB/base model self-check: `{validation['existing_model_self_check']}` (used only as a library check here).",
        f"- Random seed for row-cluster bootstrap: `{RANDOM_SEED}`; bootstrap samples follow the repository advanced-analysis setting `{advanced.BOOTSTRAP_SAMPLES}`.",
        "- `provenance.json` records input hashes, units, code identity, software versions and claim scope; `artifact_manifest.sha256` covers the model-free package and intentionally excludes the separately manifested `apparent_pb/` package.",
    ]
    (RESULTS / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def postprocess_saved_results() -> None:
    """Add derived sensitivities without repeating raw JPK reconstruction."""
    required = (
        "all_pixel_force_curves.npz", "map_inventory_QC.csv", "pixel_QC.csv",
        "normalized_force_pixels.csv", "palindrome_history_summary.csv",
        "same_speed_full_interval_changes.csv", "time_adjusted_models.csv",
        "baseline_coupling_associations.csv", "water_contact_sensitivity_curves.csv",
        "provenance.json",
    )
    missing = [name for name in required if not (RESULTS / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Postprocess inputs are missing: {missing}")

    snap_models, snap_pairs = snap_safe_sensitivity_from_saved()
    write_csv(RESULTS / "snap_safe_time_models.csv", snap_models)
    write_csv(RESULTS / "snap_safe_palindrome_summary.csv", snap_pairs)
    records = read_csv(RESULTS / "map_inventory_QC.csv")
    contact = read_csv(RESULTS / "water_contact_sensitivity_curves.csv")
    alternatives = {
        str(span): float(np.median([
            row[f"sensitivity_{int(span)}nm_nm_per_V"] for row in contact
            if row["retained_for_global_consensus"]
        ]))
        for span in pilot.MAP_CONTACT_CHECKS_NM
    }
    with (RESULTS / "provenance.json").open(encoding="utf-8") as stream:
        provenance = json.load(stream)
    validation = provenance["validation"]
    validation.update(
        {
            "all_snap_safe_designs_full_rank": all(
                row["design_rank"] == row["design_parameters"] for row in snap_models
            ),
            "minimum_snap_safe_pixels_per_map": min(
                row["minimum_snap_safe_pixels_per_map"] for row in snap_models
            ),
        }
    )
    if not validation["all_snap_safe_designs_full_rank"]:
        raise RuntimeError("At least one snap-safe time/speed design lost rank")
    provenance["postprocessed_utc"] = datetime.now(timezone.utc).isoformat()
    provenance["snap_safe_sensitivity"] = (
        "exclude each pixel unless the 5 nm bin lower edge is at least 5 nm "
        "beyond its detected approach snap; normalized branch applies the "
        "guard to its 20 nm denominator bin"
    )
    code_key = Path(__file__).resolve().relative_to(ROOT).as_posix()
    provenance["code_sha256"][code_key] = sha256_file(Path(__file__).resolve())
    (RESULTS / "provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    render_report(
        provenance["water_InvOLS_nm_per_V"], alternatives, records,
        read_csv(RESULTS / "palindrome_history_summary.csv"),
        read_csv(RESULTS / "same_speed_full_interval_changes.csv"),
        read_csv(RESULTS / "time_adjusted_models.csv"),
        read_csv(RESULTS / "baseline_coupling_associations.csv"),
        snap_models, snap_pairs, validation,
    )
    create_manifest()
    print(f"Refreshed saved postprocessing in {RESULTS}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--postprocess-only", action="store_true",
        help="reuse validated saved reconstruction and refresh derived sensitivities/report",
    )
    arguments = parser.parse_args()
    if arguments.postprocess_only:
        postprocess_saved_results()
        return
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)

    raw_inventory, primary_inventory = build_inventory()
    write_csv(RESULTS / "raw_archive_inventory.csv", raw_inventory)
    sources, records = prepare_sources(primary_inventory)
    contact_rows, sensitivity, alternatives = attach_calibration(sources, records)
    print(f"Current-water hard-contact InvOLS: {sensitivity * 1e9:.8f} nm/V", flush=True)
    curve_rows, map_rows, force_rows, matrices, reference_error = reconstruct(sources, records)

    pixel_slices, map_slices = selected_force_tables(sources, records, matrices)
    normalized_pixels, normalized_summary = normalized_shape(sources, records, matrices, curve_rows)
    pair_rows, pair_aggregate = palindrome_analysis(sources, records, matrices, rng)
    endpoint_rows = same_speed_endpoint_analysis(sources, records, matrices, rng)
    model_rows = time_models(records, map_rows, map_slices, normalized_summary)
    association_rows = qc_associations(map_rows, map_slices, normalized_summary)

    output_tables = {
        "water_contact_sensitivity_curves.csv": contact_rows,
        "map_inventory_QC.csv": map_rows,
        "pixel_QC.csv": curve_rows,
        "map_force_curves.csv": force_rows,
        "pixel_force_slices.csv": pixel_slices,
        "map_force_slices.csv": map_slices,
        "palindrome_pair_history.csv": pair_rows,
        "palindrome_history_summary.csv": pair_aggregate,
        "same_speed_full_interval_changes.csv": endpoint_rows,
        "time_adjusted_models.csv": model_rows,
        "normalized_force_pixels.csv": normalized_pixels,
        "normalized_force_summary.csv": normalized_summary,
        "baseline_coupling_associations.csv": association_rows,
    }
    for filename, rows in output_tables.items():
        write_csv(RESULTS / filename, rows)

    force_stack = np.stack([
        matrices[source.path.relative_to(ROOT).as_posix() + "|line"] for source in sources
    ])
    constant_stack = np.stack([
        matrices[source.path.relative_to(ROOT).as_posix() + "|constant"] for source in sources
    ])
    np.savez_compressed(
        RESULTS / "all_pixel_force_curves.npz",
        distance_nm=pilot.BIN_CENTERS_NM,
        global_primary_order=np.arange(1, 37),
        instrument_scan_number=np.asarray([record["instrument_scan_number"] for record in records]),
        location=np.asarray([record["location"] for record in records]),
        location_order=np.asarray([record["location_order"] for record in records]),
        point_index=np.arange(64),
        line_pN=force_stack,
        constant_pN=constant_stack,
    )

    snap_models, snap_pairs = snap_safe_sensitivity_from_saved()
    write_csv(RESULTS / "snap_safe_time_models.csv", snap_models)
    write_csv(RESULTS / "snap_safe_palindrome_summary.csv", snap_pairs)

    plot_chronology(records, map_slices)
    plot_palindrome(pair_rows)
    plot_qc(records, map_rows)
    plot_normalized(normalized_summary)

    finite_target = force_stack[:, :, TARGET_INDICES]
    validation = {
        "raw_archives": len(raw_inventory),
        "raw_total_bytes": int(sum(row["size_bytes"] for row in raw_inventory)),
        "raw_zip_crc_pass": all(row["zip_crc_ok"] for row in raw_inventory),
        "primary_maps": len(sources),
        "excluded_test_maps": sum(row["analysis_role"] == "test_excluded" for row in raw_inventory),
        "excluded_survey_reposition_maps": sum(row["analysis_role"] == "survey_reposition_excluded" for row in raw_inventory),
        "approach_curves": sum(len(source.curves) for source in sources),
        "approach_parser_skips": sum(source.skipped_curves for source in sources),
        "retract_curves": int(sum(row["retract_curves"] for row in map_rows)),
        "retract_parser_skips": int(sum(row["retract_skipped_curves"] for row in map_rows)),
        "force_array_shape": list(force_stack.shape),
        "finite_target_bin_fraction": float(np.mean(np.isfinite(finite_target))),
        "fixed_pixel_reconstruction_max_abs_difference_pN": reference_error,
        "global_water_InvOLS_nm_per_V": sensitivity * 1e9,
        "retained_contact_fits": sum(row["retained_for_global_consensus"] for row in contact_rows),
        "total_contact_fits": len(contact_rows),
        "all_time_designs_full_rank": all(
            row["design_rank"] == row["design_parameters"] for row in model_rows
        ),
        "maximum_time_design_condition_number": max(row["design_condition_number"] for row in model_rows),
        "all_snap_safe_designs_full_rank": all(
            row["design_rank"] == row["design_parameters"] for row in snap_models
        ),
        "minimum_snap_safe_pixels_per_map": min(
            row["minimum_snap_safe_pixels_per_map"] for row in snap_models
        ),
        "existing_model_self_check": base.self_test(),
    }
    if validation["approach_curves"] != 2284 or validation["retract_curves"] != 2284:
        raise RuntimeError(f"Unexpected primary curve count: {validation}")
    if validation["approach_parser_skips"] or validation["retract_parser_skips"]:
        raise RuntimeError("Primary parser skipped curves")
    if force_stack.shape != (36, 64, 100) or constant_stack.shape != force_stack.shape:
        raise RuntimeError("Unexpected force-array shape")
    if not validation["all_time_designs_full_rank"]:
        raise RuntimeError("At least one time/speed design lost rank")
    if not validation["all_snap_safe_designs_full_rank"]:
        raise RuntimeError("At least one snap-safe time/speed design lost rank")

    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": "10-09-26 two-location model-free time/history analysis",
        "sample_medium": "pure water",
        "sample_medium_authority": "user confirmed in current task",
        "temperature_C": TEMPERATURE_C,
        "temperature_authority": "user confirmed in current task",
        "cantilever": "D4",
        "spring_constant_N_per_m": SPRING_CONSTANT_N_PER_M,
        "spring_constant_repeatability_sd_N_per_m": SPRING_CONSTANT_REPEATABILITY_SD_N_PER_M,
        "spring_constant_authority": "user confirmed reuse of independent 27-08-26 D4 air calibration",
        "embedded_calibration_not_used": {
            "InvOLS_nm_per_V": EMBEDDED_INVOLS_NM_PER_V,
            "spring_constant_N_per_m": EMBEDDED_SPRING_CONSTANT_N_PER_M,
        },
        "probe_radius_m": PROBE_RADIUS_M,
        "probe_radius_authority": "user confirmed reuse in current task",
        "water_InvOLS_nm_per_V": sensitivity * 1e9,
        "water_InvOLS_method": "global median after per-map hard-contact validity and MAD retention over 36 primary maps; scan 36 retains its 44 saved curves",
        "hard_contact": {
            "primary_span_nm": pilot.MAP_CONTACT_SPAN_NM,
            "sensitivity_spans_nm": list(pilot.MAP_CONTACT_CHECKS_NM),
            "r2_min": pilot.MAP_CONTACT_R2_MIN,
            "valid_range_nm_per_V": [35.0, 90.0],
            "minimum_valid_curves_per_complete_map": 56,
            "partial_scan_36_minimum_valid_fraction": 0.875,
        },
        "force_definition": "F_pN = k_N_per_m * InvOLS_m_per_V * far-corrected_vDeflection_V * 1e12",
        "separation_definition": "D_nm = (measuredHeight_m + InvOLS_m_per_V * corrected_vDeflection_V - hard_contact_height_m) * 1e9",
        "far_baseline": "primary robust line on initial max(80,ceil(0.2*N)) samples; secondary far-window constant",
        "binning": "5 nm precontact raw-sample medians; no interpolation or extrapolation",
        "primary_target_distances_nm": TARGET_DISTANCES_NM.tolist(),
        "experimental_unit": "map/palindrome pair; 8x8 pixels are paired spatial observations, not independent experimental replicates",
        "location_scope": "time and same-pixel inference within A or B only; intervening geometry change excludes continuous same-pixel A-to-B inference",
        "normalization": "(F(D)-F(200))/(F(20)-F(200)); require denominator > max(3*per-pixel far noise,50 pN)",
        "snap_safe_sensitivity": "exclude each pixel unless the 5 nm bin lower edge is at least 5 nm beyond its detected approach snap; normalized branch applies the guard to its 20 nm denominator bin",
        "time_model": "map-median OLS with HC3 covariance; time scaled to full within-location span; measured gap speed included; Newey-West HAC(1/2) reported as ordered-map sensitivity",
        "random_seed": RANDOM_SEED,
        "bootstrap_samples": advanced.BOOTSTRAP_SAMPLES,
        "validation": validation,
        "raw_input_sha256": {row["source"]: row["sha256"] for row in raw_inventory},
        "code_sha256": {
            path.relative_to(ROOT).as_posix(): sha256_file(path)
            for path in (
                Path(__file__).resolve(), Path(base.__file__).resolve(), Path(pilot.__file__).resolve(),
                Path(chronology.__file__).resolve(), Path(reference.__file__).resolve(), Path(advanced.__file__).resolve(),
            )
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
    }
    (RESULTS / "provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    render_report(
        sensitivity * 1e9, alternatives, records, pair_aggregate, endpoint_rows,
        model_rows, association_rows, snap_models, snap_pairs, validation,
    )
    create_manifest()
    print(json.dumps(validation, indent=2), flush=True)
    print(f"Wrote {RESULTS}", flush=True)


if __name__ == "__main__":
    main()

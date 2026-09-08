#!/usr/bin/env python3
"""Reproduce the 08-09-26 water-map quality comparison.

Reads the saved reconstruction and surface-fit tables. Writes only the
surface_fit/quality_comparison directory; the primary fit artifacts are inputs.
"""

from __future__ import annotations

from datetime import datetime, timezone
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.stats import spearmanr

import analyze_08_09_26_fixed_pixel as chronology
import fit_08_09_26_water as fit
import fit_glycerol_surface_forces as base


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "analysis" / "water_08_09_26_results"
SURFACE = DATA / "surface_fit"
OUT = SURFACE / "quality_comparison"
MAD_SCALE = 1.4826


def read(path: Path) -> list[dict]:
    return chronology.read_derived_csv(path)


def mad(values) -> float:
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    if not len(a):
        return float("nan")
    return float(MAD_SCALE * np.median(np.abs(a - np.median(a))))


def median(values) -> float:
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    return float(np.median(a)) if len(a) else float("nan")


def percentile(values, q: float) -> float:
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    return float(np.percentile(a, q)) if len(a) else float("nan")


def relative_rmse(row: dict) -> float:
    signal = float(row["force_signal_pN"])
    return 100.0 * float(row["residual_rmse_pN"]) / signal if signal > 0 else float("nan")


def spearman_row(scope: str, speed, rows: list[dict], key: str) -> dict:
    x = np.array([r["acquisition_order"] for r in rows], dtype=float)
    y = np.array([r[key] for r in rows], dtype=float)
    result = spearmanr(x, y)
    return {"scope": scope, "nominal_speed_um_per_s": speed,
            "parameter": key, "n_maps": len(rows),
            "spearman_rho": float(result.statistic),
            "inference_boundary": "descriptive_serial_association; adjacent maps are not independent experiments"}


def common_25_fits(inputs: list[dict], primary: dict[int, dict]) -> tuple[list[dict], list[dict]]:
    results, predictions = [], []
    for order in range(19, 37):
        source = primary[order]
        selected = [r for r in inputs if r["acquisition_order"] == order
                    and r["baseline_method"] == "linear_drift_corrected"
                    and 25 <= r["distance_nm"] <= 250 and r["pixel_count"] >= 48
                    and np.isfinite(r["force_median_pN"])
                    and np.isfinite(r["fit_weight_scale_pN"])]
        if len(selected) != 46:
            raise AssertionError(f"Map {order}: expected 46 bins at 25-250 nm")
        task = {key: source[key] for key in (
            "acquisition_order", "block", "source", "map_start_time", "map_end_time",
            "map_midpoint_time", "nominal_speed_um_per_s", "safe_dmin_nm",
            "detected_snap_max_nm", "detected_snap_median_nm", "detected_snap_p95_nm",
            "snap_detected_pixels", "noise_floor_pN")}
        task.update({"variant": "common_25_250nm", "baseline_method": "linear_drift_corrected",
                     "dmin_nm": 25.0, "dmax_nm": 250.0,
                     "minimum_pixels_per_fitted_bin": min(r["pixel_count"] for r in selected),
                     "hydro_speed_um_s": 0.0, "no_slip_correction_applied_to_pixels": False,
                     "distance_nm": [r["distance_nm"] for r in selected],
                     "force_pN": [r["force_median_pN"] for r in selected],
                     "sigma_pN": [r["fit_weight_scale_pN"] for r in selected]})
        result = fit.fit_one(task)
        results.append(result)
        parameters = fit.result_parameters(result)
        d = np.asarray(task["distance_nm"], dtype=float)
        observed = np.asarray(task["force_pN"], dtype=float)
        predicted = fit.predict(d, parameters)
        sigma = np.asarray(task["sigma_pN"], dtype=float)
        for distance, obs, pred, weight in zip(d, observed, predicted, sigma):
            predictions.append({"acquisition_order": order,
                "nominal_speed_um_per_s": result["nominal_speed_um_per_s"],
                "distance_nm": float(distance), "observed_force_pN": float(obs),
                "predicted_force_pN": float(pred), "residual_observed_minus_fit_pN": float(obs-pred),
                "fit_weight_scale_pN": float(weight)})
    return results, predictions


def spatial_tables(distance: np.ndarray, line: np.ndarray, constant: np.ndarray,
                   primary: dict[int, dict]) -> tuple[list[dict], dict[tuple[int, str], dict], list[dict]]:
    details, summaries, references = [], {}, []
    common = (distance >= 80) & (distance <= 250)
    reference_window = (distance >= 200) & (distance <= 250)
    for order in range(1, 37):
        for method, cube in (("linear_drift_corrected", line),
                             ("far_constant_referenced", constant)):
            values = cube[order-1]
            reference = np.nanmedian(values[:, reference_window], axis=1)
            referenced = values - reference[:, None]
            for point, level in enumerate(reference):
                references.append({"acquisition_order": order, "point_index": point,
                    "baseline_method": method, "reference_window_min_nm": 200.0,
                    "reference_window_max_nm": 250.0,
                    "per_pixel_reference_median_pN": float(level),
                    "finite_reference_bins": int(np.sum(np.isfinite(values[point, reference_window])))})
            per_d_raw, per_d_ref = [], []
            for j in np.flatnonzero(common):
                raw_mad = mad(values[:, j])
                ref_mad = mad(referenced[:, j])
                per_d_raw.append(raw_mad)
                per_d_ref.append(ref_mad)
                details.append({"acquisition_order": order, "baseline_method": method,
                    "distance_nm": float(distance[j]),
                    "finite_pixel_count": int(np.sum(np.isfinite(values[:, j]))),
                    "force_pixel_median_pN": median(values[:, j]),
                    "spatial_mad_pN": raw_mad,
                    "reference_subtracted_spatial_mad_pN": ref_mad,
                    "reference_definition": "subtract each pixel median force over 200-250 nm before spatial MAD"})
            pmin = float(primary[order]["dmin_nm"])
            pwindow = (distance >= pmin) & (distance <= 250)
            primary_raw = [mad(values[:, j]) for j in np.flatnonzero(pwindow)]
            summaries[(order, method)] = {
                "common_80_250_spatial_mad_distance_median_pN": median(per_d_raw),
                "common_80_250_reference_subtracted_spatial_mad_distance_median_pN": median(per_d_ref),
                "primary_window_spatial_mad_distance_median_pN": median(primary_raw),
                "aggregation": "MAD across pixels at each distance, then median across distance bins"}
    return details, summaries, references


def build_map_rows(inventory: list[dict], pixels: list[dict], primary: dict[int, dict],
                   variants: dict[tuple[int, str], dict], spatial: dict[tuple[int, str], dict]) -> list[dict]:
    rows = []
    for meta in inventory:
        order = int(meta["acquisition_order"])
        curve = [r for r in pixels if r["acquisition_order"] == order]
        far_noise = [float(r["far_noise_pN"]) for r in curve]
        contact_r2 = [float(r["contact_fit_R2"]) for r in curve
                      if np.isfinite(float(r["contact_fit_R2"]))]
        inv_ols = [float(r["contact_invOLS_nm_per_V"]) for r in curve
                   if np.isfinite(float(r["contact_invOLS_nm_per_V"]))]
        snap = [float(r["approach_snap_distance_nm"]) for r in curve
                if r["approach_snap_detected"] and np.isfinite(float(r["approach_snap_distance_nm"]))]
        p = primary[order]
        c = variants[(order, "common_80_250nm")]
        line = spatial[(order, "linear_drift_corrected")]
        const = spatial[(order, "far_constant_referenced")]
        rows.append({
            "acquisition_order": order, "block": int(meta["block"]),
            "nominal_speed_um_per_s": float(meta["nominal_speed_um_per_s"]),
            "source": meta["source"], "source_sha256": meta["sha256"],
            "approach_curves": int(meta["approach_curves"]),
            "approach_skipped_curves": int(meta["skipped_curves"]),
            "retract_curves": int(meta["retract_curves"]),
            "retract_skipped_curves": int(meta["retract_skipped_curves"]),
            "instrument_failure_flag_count": int(meta["instrument_failure_flag_count"]),
            "valid_contact_curves": int(meta["valid_contact_curves"]),
            "contact_fit_r2_median": median(contact_r2), "contact_fit_r2_min": float(np.min(contact_r2)),
            "contact_invOLS_median_nm_per_V": float(meta["map_contact_invOLS_median_nm_per_V"]),
            "contact_invOLS_mad_nm_per_V": float(meta["map_contact_invOLS_mad_nm_per_V"]),
            "terminal_load_median_nN": float(meta["terminal_load_median_nN"]),
            "far_noise_pixel_median_pN": median(far_noise), "far_noise_pixel_p95_pN": percentile(far_noise, 95),
            "far_noise_pixel_max_pN": float(np.max(far_noise)),
            "far_noise_pixels_over_30pN": int(np.sum(np.asarray(far_noise) > 30)),
            "far_noise_pixels_over_50pN": int(np.sum(np.asarray(far_noise) > 50)),
            "far_slope_median_pN_per_100nm": float(meta["far_slope_median_pN_per_100nm"]),
            "far_slope_mad_pN_per_100nm": float(meta["far_slope_mad_pN_per_100nm"]),
            "snap_detected_pixels": len(snap), "snap_detected_fraction": len(snap) / 64,
            "snap_distance_median_nm": float(meta["approach_snap_distance_median_nm"]),
            "snap_distance_p95_nm": percentile(snap, 95),
            "snap_distance_max_nm": float(np.max(snap)) if snap else float("nan"),
            "primary_dmin_nm": float(p["dmin_nm"]), "primary_dmax_nm": float(p["dmax_nm"]),
            "primary_r2": float(p["r2"]), "primary_rmse_pN": float(p["residual_rmse_pN"]),
            "primary_force_signal_pN": float(p["force_signal_pN"]),
            "primary_relative_rmse_percent": relative_rmse(p),
            "primary_lambda_D_apparent_nm": float(p["lambda_D_apparent_nm"]),
            "primary_potential_magnitude_apparent_mV": float(p["surface_potential_magnitude_apparent_mV"]),
            "primary_fit_quality_pass": bool(p["fit_quality_pass"]), "primary_quality_flags": p["quality_flags"],
            "common_80_250_r2": float(c["r2"]), "common_80_250_rmse_pN": float(c["residual_rmse_pN"]),
            "common_80_250_force_signal_pN": float(c["force_signal_pN"]),
            "common_80_250_relative_rmse_percent": relative_rmse(c),
            "common_80_250_fit_quality_pass": bool(c["fit_quality_pass"]),
            "common_80_250_quality_flags": c["quality_flags"],
            "common_80_250_line_spatial_mad_pN": line["common_80_250_spatial_mad_distance_median_pN"],
            "common_80_250_line_reference_subtracted_spatial_mad_pN": line["common_80_250_reference_subtracted_spatial_mad_distance_median_pN"],
            "common_80_250_constant_spatial_mad_pN": const["common_80_250_spatial_mad_distance_median_pN"],
            "common_80_250_constant_reference_subtracted_spatial_mad_pN": const["common_80_250_reference_subtracted_spatial_mad_distance_median_pN"],
            "primary_window_line_spatial_mad_pN": line["primary_window_spatial_mad_distance_median_pN"],
            "primary_window_constant_spatial_mad_pN": const["primary_window_spatial_mad_distance_median_pN"]})
    return rows


def aggregate_row(scope: str, label, rows: list[dict], common25: dict[int, dict]) -> dict:
    def med(key):
        return median([r[key] for r in rows])
    lam = np.array([r["primary_lambda_D_apparent_nm"] for r in rows], dtype=float)
    psi = np.array([r["primary_potential_magnitude_apparent_mV"] for r in rows], dtype=float)
    common = [common25[r["acquisition_order"]] for r in rows if r["acquisition_order"] in common25]
    def common_med(key):
        return median([r[key] for r in common]) if common else None
    return {"scope": scope, "label": label, "first_map": min(r["acquisition_order"] for r in rows),
        "last_map": max(r["acquisition_order"] for r in rows), "n_maps": len(rows),
        "nominal_speed_um_per_s": rows[0]["nominal_speed_um_per_s"] if scope == "phase_speed" else "",
        "approach_curves_sum": sum(r["approach_curves"] for r in rows),
        "approach_skipped_curves_sum": sum(r["approach_skipped_curves"] for r in rows),
        "retract_curves_sum": sum(r["retract_curves"] for r in rows),
        "retract_skipped_curves_sum": sum(r["retract_skipped_curves"] for r in rows),
        "instrument_failure_flags_sum": sum(r["instrument_failure_flag_count"] for r in rows),
        "far_noise_map_median_pN": med("far_noise_pixel_median_pN"),
        "contact_fit_r2_map_median": med("contact_fit_r2_median"),
        "contact_invOLS_map_median_nm_per_V": med("contact_invOLS_median_nm_per_V"),
        "terminal_load_map_median_nN": med("terminal_load_median_nN"),
        "snap_distance_map_median_nm": med("snap_distance_median_nm"),
        "primary_dmin_min_nm": min(r["primary_dmin_nm"] for r in rows),
        "primary_dmin_max_nm": max(r["primary_dmin_nm"] for r in rows),
        "primary_r2_map_median": med("primary_r2"),
        "primary_rmse_map_median_pN": med("primary_rmse_pN"),
        "primary_force_signal_map_median_pN": med("primary_force_signal_pN"),
        "primary_relative_rmse_map_median_percent": med("primary_relative_rmse_percent"),
        "primary_quality_pass_count": sum(bool(r["primary_fit_quality_pass"]) for r in rows),
        "primary_lambda_sample_cv_percent": float(100*np.std(lam, ddof=1)/np.mean(lam)) if len(lam) > 1 else float("nan"),
        "primary_potential_sample_cv_percent": float(100*np.std(psi, ddof=1)/np.mean(psi)) if len(psi) > 1 else float("nan"),
        "common_80_250_r2_map_median": med("common_80_250_r2"),
        "common_80_250_rmse_map_median_pN": med("common_80_250_rmse_pN"),
        "common_80_250_force_signal_map_median_pN": med("common_80_250_force_signal_pN"),
        "common_80_250_relative_rmse_map_median_percent": med("common_80_250_relative_rmse_percent"),
        "common_80_250_quality_pass_count": sum(bool(r["common_80_250_fit_quality_pass"]) for r in rows),
        "common_80_250_line_spatial_mad_map_median_pN": med("common_80_250_line_spatial_mad_pN"),
        "common_80_250_line_reference_subtracted_spatial_mad_map_median_pN": med("common_80_250_line_reference_subtracted_spatial_mad_pN"),
        "common_80_250_constant_spatial_mad_map_median_pN": med("common_80_250_constant_spatial_mad_pN"),
        "common_80_250_constant_reference_subtracted_spatial_mad_map_median_pN": med("common_80_250_constant_reference_subtracted_spatial_mad_pN"),
        "common_25_250_fit_count": len(common),
        "common_25_250_lambda_map_median_nm": common_med("lambda_D_apparent_nm"),
        "common_25_250_potential_map_median_mV": common_med("surface_potential_magnitude_apparent_mV"),
        "common_25_250_r2_map_median": common_med("r2"),
        "common_25_250_rmse_map_median_pN": common_med("residual_rmse_pN"),
        "common_25_250_quality_pass_count": sum(bool(r["fit_quality_pass"]) for r in common) if common else None,
        "aggregation": "compute each metric per map, then take the median over maps; CV uses sample SD with ddof=1 divided by mean"}


def baseline_comparison(distance: np.ndarray, line: np.ndarray, constant: np.ndarray,
                        inventory: list[dict]) -> list[dict]:
    common = (distance >= 80) & (distance <= 250)
    j100 = int(np.flatnonzero(distance == 100)[0])
    rows = []
    for order in range(1, 37):
        delta = line[order-1, :, common] - constant[order-1, :, common]
        meta = inventory[order-1]
        rows.append({"acquisition_order": order, "block": int(meta["block"]),
            "nominal_speed_um_per_s": float(meta["nominal_speed_um_per_s"]),
            "line_map_median_force_at_100nm_pN": median(line[order-1, :, j100]),
            "constant_map_median_force_at_100nm_pN": median(constant[order-1, :, j100]),
            "line_minus_constant_at_100nm_pN": median(line[order-1, :, j100])-median(constant[order-1, :, j100]),
            "common_80_250_line_minus_constant_pixel_bin_median_pN": median(delta.ravel()),
            "common_80_250_line_minus_constant_pixel_bin_mad_pN": mad(delta.ravel()),
            "common_80_250_line_minus_constant_pixel_bin_p95_abs_pN": percentile(np.abs(delta.ravel()), 95),
            "independence_boundary": "pixel-distance cells are repeated, spatially and serially correlated measurements"})
    return rows


def selected_diagnostics(map_rows: list[dict], variants: dict[tuple[int, str], dict]) -> list[dict]:
    by_order = {r["acquisition_order"]: r for r in map_rows}
    line14 = variants[(14, "primary")]
    const14 = variants[(14, "constant_baseline")]
    r25 = by_order[25]
    entries = [
        (14, "baseline_sensitivity", "primary_line_lambda_D_apparent_nm", line14["lambda_D_apparent_nm"], "nm"),
        (14, "baseline_sensitivity", "constant_lambda_D_apparent_nm", const14["lambda_D_apparent_nm"], "nm"),
        (14, "baseline_sensitivity", "constant_vs_line_lambda_relative_change_percent", 100*(const14["lambda_D_apparent_nm"]/line14["lambda_D_apparent_nm"]-1), "percent"),
        (14, "baseline_sensitivity", "primary_line_potential_magnitude_apparent_mV", line14["surface_potential_magnitude_apparent_mV"], "mV"),
        (14, "baseline_sensitivity", "constant_potential_magnitude_apparent_mV", const14["surface_potential_magnitude_apparent_mV"], "mV"),
        (14, "baseline_sensitivity", "constant_vs_line_potential_relative_change_percent", 100*(const14["surface_potential_magnitude_apparent_mV"]/line14["surface_potential_magnitude_apparent_mV"]-1), "percent"),
        (14, "baseline_sensitivity", "primary_line_r2", line14["r2"], "dimensionless"),
        (14, "baseline_sensitivity", "constant_r2", const14["r2"], "dimensionless"),
        (25, "spatial_far_field_anomaly", "primary_window_line_spatial_mad_distance_median_pN", r25["primary_window_line_spatial_mad_pN"], "pN"),
        (25, "spatial_far_field_anomaly", "far_slope_mad_pN_per_100nm", r25["far_slope_mad_pN_per_100nm"], "pN_per_100nm"),
        (25, "spatial_far_field_anomaly", "far_noise_pixel_median_pN", r25["far_noise_pixel_median_pN"], "pN"),
        (25, "spatial_far_field_anomaly", "far_noise_pixel_p95_pN", r25["far_noise_pixel_p95_pN"], "pN"),
        (25, "spatial_far_field_anomaly", "far_noise_pixels_over_50pN", r25["far_noise_pixels_over_50pN"], "pixels")]
    return [{"acquisition_order": order, "category": category, "metric": metric,
             "value": float(value), "unit": unit} for order, category, metric, value, unit in entries]


def main() -> None:
    fit.configure_model()
    OUT.mkdir(parents=True, exist_ok=True)
    parent_manifest = SURFACE / "artifact_manifest.sha256"
    parent_manifest_before = base.sha256_file(parent_manifest)

    inventory = read(DATA / "map_inventory_QC.csv")
    pixels = read(DATA / "pixel_QC.csv")
    inputs = read(SURFACE / "map_median_fit_inputs.csv")
    all_fits = read(SURFACE / "fit_parameters_all_variants.csv")
    sensitivity = read(SURFACE / "parameter_sensitivity.csv")
    if len(inventory) != 36 or len(pixels) != 2304 or len(inputs) != 10800 or len(all_fits) != 216:
        raise AssertionError("Unexpected saved input row counts")
    primary = {int(r["acquisition_order"]): r for r in all_fits if r["variant"] == "primary"}
    variants = {(int(r["acquisition_order"]), r["variant"]): r for r in all_fits}
    if set(primary) != set(range(1, 37)):
        raise AssertionError("Primary fit map coverage changed")

    with np.load(DATA / "all_pixel_force_curves.npz") as archive:
        distance = archive["distance_nm"].copy()
        line = archive["line_pN"].copy()
        constant = archive["constant_pN"].copy()
        if not np.array_equal(archive["acquisition_order"], np.arange(1, 37)):
            raise AssertionError("NPZ map ordering changed")
        if not np.array_equal(archive["point_index"], np.arange(64)):
            raise AssertionError("NPZ pixel ordering changed")

    common25_rows, common25_predictions = common_25_fits(inputs, primary)
    common25 = {int(r["acquisition_order"]): r for r in common25_rows}
    spatial_detail, spatial_summary, pixel_references = spatial_tables(
        distance, line, constant, primary)
    map_rows = build_map_rows(inventory, pixels, primary, variants, spatial_summary)
    baseline_rows = baseline_comparison(distance, line, constant, inventory)
    baseline_by_order = {r["acquisition_order"]: r for r in baseline_rows}

    group_rows = []
    for label, lo, hi in (("maps_1_18", 1, 18), ("maps_19_36", 19, 36)):
        group_rows.append(aggregate_row("phase", label,
            [r for r in map_rows if lo <= r["acquisition_order"] <= hi], common25))
    block_rows = []
    for block in range(1, 7):
        subset = [r for r in map_rows if r["block"] == block]
        row = aggregate_row("six_map_block", f"block_{block}", subset, common25)
        orders = {r["acquisition_order"] for r in subset}
        row["far_noise_pooled_pixel_median_pN"] = median(
            [r["far_noise_pN"] for r in pixels if r["acquisition_order"] in orders])
        row["line_force_at_100nm_map_median_pN"] = median(
            [baseline_by_order[r["acquisition_order"]]["line_map_median_force_at_100nm_pN"] for r in subset])
        row["constant_force_at_100nm_map_median_pN"] = median(
            [baseline_by_order[r["acquisition_order"]]["constant_map_median_force_at_100nm_pN"] for r in subset])
        block_rows.append(row)
    speed_rows = []
    for label, lo, hi in (("maps_1_18", 1, 18), ("maps_19_36", 19, 36)):
        for speed in (1.0, 2.0, 4.0):
            subset = [r for r in map_rows if lo <= r["acquisition_order"] <= hi
                      and r["nominal_speed_um_per_s"] == speed]
            speed_rows.append(aggregate_row("phase_speed", f"{label}_speed_{speed:g}", subset, common25))

    trend_rows = []
    for key in ("lambda_D_apparent_nm", "surface_potential_magnitude_apparent_mV"):
        trend_rows.append(spearman_row("maps_19_36_all_speeds", "", common25_rows, key))
        for speed in (1.0, 2.0, 4.0):
            subset = [r for r in common25_rows if r["nominal_speed_um_per_s"] == speed]
            trend_rows.append(spearman_row("maps_19_36_within_speed", speed, subset, key))

    method_rows = []
    for label, lo, hi in (("maps_1_18", 1, 18), ("maps_19_36", 19, 36)):
        for variant in ("constant_baseline", "lower_plus_10nm", "upper_200nm", "common_80_250nm", "no_slip_hydrodynamic"):
            subset = [r for r in sensitivity if lo <= r["acquisition_order"] <= hi and r["variant"] == variant]
            lchange = [r["lambda_relative_change_percent"] for r in subset]
            pchange = [r["potential_relative_change_percent"] for r in subset]
            method_rows.append({"phase": label, "variant": variant, "n_maps": len(subset),
                "fit_quality_pass_count": sum(bool(r["fit_quality_pass"]) for r in subset),
                "lambda_relative_change_median_percent": median(lchange),
                "lambda_absolute_relative_change_median_percent": median(np.abs(lchange)),
                "potential_relative_change_median_percent": median(pchange),
                "potential_absolute_relative_change_median_percent": median(np.abs(pchange)),
                "reference": "each variant relative to the same map primary fit"})

    selected = selected_diagnostics(map_rows, variants)
    m25 = next(r for r in map_rows if r["acquisition_order"] == 25)
    other_post = [r for r in map_rows if 19 <= r["acquisition_order"] <= 36 and r["acquisition_order"] != 25]
    c25_lam = [r["lambda_D_apparent_nm"] for r in common25_rows]
    c25_psi = [r["surface_potential_magnitude_apparent_mV"] for r in common25_rows]
    summary = {
        "comparison": "first 18 versus last 18 of 36 sequential 8x8 water force maps",
        "headline": "Later maps have closer snap-in, lower far-field noise, and lower spatial MAD, making their recovered near-field branch better suited to the present PB fit. The shared 80-250 nm tail does not constrain PB parameters reliably in either group.",
        "group_summary": {r["label"]: r for r in group_rows},
        "common_25_250_maps_19_36": {
            "maps": 18, "bins_per_map": 46,
            "quality_pass_count": sum(bool(r["fit_quality_pass"]) for r in common25_rows),
            "lambda_D_median_nm": median(c25_lam), "lambda_D_min_nm": min(c25_lam), "lambda_D_max_nm": max(c25_lam),
            "potential_magnitude_median_mV": median(c25_psi), "potential_magnitude_min_mV": min(c25_psi), "potential_magnitude_max_mV": max(c25_psi)},
        "map_14_baseline_sensitivity": {r["metric"]: r["value"] for r in selected if r["acquisition_order"] == 14},
        "map_25_spatial_far_field_anomaly": {
            **{r["metric"]: r["value"] for r in selected if r["acquisition_order"] == 25},
            "other_maps_19_36_primary_window_line_spatial_mad_min_pN": min(r["primary_window_line_spatial_mad_pN"] for r in other_post),
            "other_maps_19_36_primary_window_line_spatial_mad_max_pN": max(r["primary_window_line_spatial_mad_pN"] for r in other_post),
            "other_maps_19_36_far_slope_mad_min_pN_per_100nm": min(r["far_slope_mad_pN_per_100nm"] for r in other_post),
            "other_maps_19_36_far_slope_mad_max_pN_per_100nm": max(r["far_slope_mad_pN_per_100nm"] for r in other_post)},
        "interpretation_boundaries": [
            "Maps, pixels, and adjacent distance bins are serially or spatially correlated and are not treated as independent experiments.",
            "MAD is 1.4826*median(abs(x-median(x))). Spatial summaries compute MAD across pixels at each distance, then median across distance bins, then median across maps for group summaries.",
            "Per-pixel reference subtraction removes each curve's median force over 200-250 nm before spatial MAD; it separates constant curve offsets from shape variation but does not identify the remainder as independent random noise.",
            "Primary-window and common-window fits answer different quality questions. Their force-signal values must not be interpreted as force changes at a fixed distance.",
            "Fit pass and numerical consistency do not establish equilibrium PB validity, a causal liquid change, or experimental confidence intervals."]}

    outputs = {
        "common_25_250_fit_parameters.csv": common25_rows,
        "common_25_250_fit_predictions.csv": common25_predictions,
        "quality_per_map.csv": map_rows,
        "quality_group_summary.csv": group_rows,
        "quality_six_map_block_summary.csv": block_rows,
        "quality_phase_speed_summary.csv": speed_rows,
        "common_25_250_trend_diagnostics.csv": trend_rows,
        "method_sensitivity_summary.csv": method_rows,
        "spatial_spread_by_distance.csv": spatial_detail,
        "per_pixel_200_250nm_reference.csv": pixel_references,
        "baseline_comparison_per_map.csv": baseline_rows,
        "selected_map_diagnostics.csv": selected}
    for name, rows in outputs.items():
        base.write_csv(OUT / name, rows)
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")

    input_paths = [DATA / "all_pixel_force_curves.npz", DATA / "map_inventory_QC.csv",
        DATA / "pixel_QC.csv", SURFACE / "map_median_fit_inputs.csv",
        SURFACE / "fit_parameters_per_map.csv", SURFACE / "fit_parameters_all_variants.csv",
        SURFACE / "parameter_sensitivity.csv", SURFACE / "provenance.json"]
    validation = {"input_rows": {"maps": len(inventory), "pixels": len(pixels),
        "fit_inputs": len(inputs), "all_primary_and_variant_fits": len(all_fits)},
        "output_rows": {name: len(rows) for name, rows in outputs.items()},
        "all_maps_capture_complete": all(r["approach_curves"] == 64 and r["approach_skipped_curves"] == 0
             and r["retract_curves"] == 64 and r["retract_skipped_curves"] == 0 for r in map_rows),
        "all_maps_instrument_failure_flag_count_zero": all(r["instrument_failure_flag_count"] == 0 for r in map_rows),
        "common_25_250_all_46_bins": all(r["n_points"] == 46 for r in common25_rows),
        "common_25_250_quality_pass_count": sum(bool(r["fit_quality_pass"]) for r in common25_rows),
        "parent_surface_fit_manifest_sha256_before": parent_manifest_before,
        "parent_manifest_hash_scope": "execution-time non-mutation guard only; the parent manifest may later be regenerated to index this comparison directory, so this is not a current-package assertion"}
    provenance = {"created_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "reproducible descriptive quality comparison; no primary fit files modified",
        "definitions": {"MAD": "1.4826*median(abs(x-median(x)))",
            "group_aggregation": "derive each metric per map, then median across maps",
            "sample_CV": "100*sample_standard_deviation(ddof=1)/mean",
            "relative_fit_residual_percent": "100*RMSE/force_signal",
            "force_signal": "max(force)-median(last max(3,floor(N/5)) fitted distance bins)",
            "spatial_MAD": "MAD across pixels separately at each distance, then median across distance bins",
            "reference_subtracted_spatial_MAD": "subtract each pixel median at 200-250 nm, then MAD across pixels at each distance, then median across distance bins",
            "fit_R2_RMSE": "unweighted force residuals; optimizer uses saved spatial-MAD weights and soft_l1 loss"},
        "fit_method": "reuse fit_08_09_26_water.fit_one, nonlinear PB Derjaguin model, line baseline, 25-250 nm, 46 bins, saved weights",
        "input_sha256": {p.relative_to(ROOT).as_posix(): base.sha256_file(p) for p in input_paths},
        "raw_capture_identity": [{"acquisition_order": int(r["acquisition_order"]),
            "source": r["source"], "sha256": r["sha256"]} for r in inventory],
        "code_sha256": {p.relative_to(ROOT).as_posix(): base.sha256_file(p)
            for p in (Path(__file__), Path(fit.__file__), Path(base.__file__), Path(chronology.__file__))},
        "model_constants": {"model": fit.MODEL, "probe_radius_m": fit.RADIUS_M,
            "temperature_C": base.TEMPERATURE_C, "epsilon_r": base.EPSILON_R[0], "Hamaker_J": base.HAMAKER_J},
        "independence_boundary": summary["interpretation_boundaries"], "validation": validation,
        "software": {"python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "scipy": scipy.__version__}}
    if base.sha256_file(parent_manifest) != parent_manifest_before:
        raise AssertionError("Parent surface-fit manifest changed during comparison")
    provenance["validation"]["parent_surface_fit_manifest_sha256_after"] = base.sha256_file(parent_manifest)
    (OUT / "provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    manifest = [f"{base.sha256_file(path)}  {path.relative_to(OUT).as_posix()}"
        for path in sorted(OUT.iterdir()) if path.is_file() and path.name != "artifact_manifest.sha256"]
    (OUT / "artifact_manifest.sha256").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUT), "files": sorted([p.name for p in OUT.iterdir()]),
        "common_25_250_quality_pass": validation["common_25_250_quality_pass_count"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()

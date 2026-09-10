#!/usr/bin/env python3
"""Secondary apparent nonlinear-PB fits for the 10-09-26 time series.

Fits are performed map by map on the median of reconstructed pixel curves.
They describe finite-speed, model-conditioned apparent lambda_D and potential
magnitude.  They do not convert the model-free time/history association into a
bulk electrolyte or electrokinetic measurement.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import platform
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import analyze_10_09_26_time_dependence as time_analysis  # noqa: E402
import fit_08_09_26_water as prior_pb  # noqa: E402
import fit_glycerol_surface_forces as base  # noqa: E402


DATA = ROOT / "analysis" / "time_dependence_10_09_26_results"
OUT = DATA / "apparent_pb"
FIGURES = OUT / "figures"
RADIUS_M = time_analysis.PROBE_RADIUS_M
K_N_PER_M = time_analysis.SPRING_CONSTANT_N_PER_M
VARIANTS = (
    "primary",
    "constant_baseline",
    "lower_plus_10nm",
    "upper_200nm",
    "common_150_250nm",
    "no_slip_hydrodynamic",
)


def read_csv(path: Path) -> list[dict]:
    return prior_pb.read_csv(path)


def json_default(value):
    """Convert NumPy scalar diagnostics without silently stringifying objects."""
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def configure_model() -> None:
    prior_pb.RADIUS_M = RADIUS_M
    prior_pb.configure_model()


def hydro_pn(distance_nm: np.ndarray, speed_um_s: float) -> np.ndarray:
    configure_model()
    eta_pa_s = base.cheng_viscosity_mPa_s(0.0, time_analysis.TEMPERATURE_C) * 1e-3
    return (
        6.0 * np.pi * eta_pa_s * RADIUS_M**2 * speed_um_s * 1e-6
        / (np.asarray(distance_nm, dtype=np.float64) * 1e-9) * 1e12
    )


def drainage_corrected_pixels(
    meta: dict, pixel_rows: list[dict], distance_nm: np.ndarray, force_matrix: np.ndarray
) -> tuple[np.ndarray, list[dict]]:
    """Apply the repository's far-line-projected ideal no-slip sensitivity."""
    source = base.load_source(ROOT / meta["source"], 0)
    if source.sha256 != meta["sha256"] or len(source.curves) != int(meta["approach_curves"]):
        raise RuntimeError("Raw source identity changed during hydrodynamic sensitivity")
    lookup = {int(row["point_index"]): row for row in pixel_rows}
    output = np.full_like(force_matrix, np.nan)
    diagnostics: list[dict] = []
    for curve in source.curves:
        point = int(curve.point_index)
        row = lookup[point]
        height = curve.measured_height_m
        far = base.fit_far_field_drift(height, curve.deflection_V)
        voltage = curve.deflection_V - base.baseline_voltage(height, curve.deflection_V, far)
        sensitivity = float(row["sensitivity_used_nm_per_V"]) * 1e-9
        contact_height_m = float(row["contact_height_um"]) * 1e-6
        raw_distance_nm = (height + sensitivity * voltage - contact_height_m) * 1e9
        speed = float(row["gap_speed_20_200nm_um_per_s"])
        if not np.isfinite(speed) or speed <= 0:
            raise RuntimeError("Measured gap speed is unavailable for no-slip sensitivity")
        far_distance = raw_distance_nm[: far.n_points]
        if np.any(far_distance <= 0):
            raise RuntimeError("No-slip far-field samples include nonpositive separation")
        scanner_coordinate_nm = (height[: far.n_points] - contact_height_m) * 1e9
        slope, intercept, _, _ = base.robust_line(
            scanner_coordinate_nm, hydro_pn(far_distance, speed)
        )
        # F[pN] / (k[N/m]*1000) is cantilever displacement in nm.
        binned_scanner_coordinate_nm = distance_nm - force_matrix[point] / (K_N_PER_M * 1000.0)
        projected = hydro_pn(distance_nm, speed) - (
            slope * binned_scanner_coordinate_nm + intercept
        )
        output[point] = force_matrix[point] - projected
        diagnostics.append(
            {
                "location": meta["location"],
                "location_order": meta["location_order"],
                "instrument_scan_number": meta["instrument_scan_number"],
                "point_index": point,
                "gap_speed_um_per_s": speed,
                "far_sample_count": far.n_points,
                "far_distance_min_nm": float(np.min(far_distance)),
                "far_distance_max_nm": float(np.max(far_distance)),
                "hydro_far_line_slope_pN_per_nm": slope,
                "hydro_far_line_intercept_pN": intercept,
            }
        )
    return output, diagnostics


def fit_task(task: dict) -> dict:
    configure_model()
    return prior_pb.fit_one(task)


def prepare() -> tuple[
    list[dict], list[dict], list[dict], list[dict], list[dict], dict
]:
    inventory = read_csv(DATA / "map_inventory_QC.csv")
    pixels = read_csv(DATA / "pixel_QC.csv")
    with np.load(DATA / "all_pixel_force_curves.npz", allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    distance = np.asarray(arrays["distance_nm"], dtype=np.float64)
    if arrays["line_pN"].shape != (36, 64, 100) or arrays["constant_pN"].shape != (36, 64, 100):
        raise RuntimeError("Unexpected upstream force-array shape")
    if not np.array_equal(arrays["global_primary_order"], np.arange(1, 37)):
        raise RuntimeError("Upstream primary map order changed")
    if any(not np.isclose(row["spring_constant_used_N_per_m"], K_N_PER_M, rtol=1e-12) for row in inventory):
        raise RuntimeError("Upstream D4 force calibration changed")

    inputs: list[dict] = []
    tasks: list[dict] = []
    exclusions: list[dict] = []
    drainage_rows: list[dict] = []
    for index, meta in enumerate(inventory):
        order = index + 1
        pixel_rows = [row for row in pixels if row["acquisition_order"] == order]
        if len(pixel_rows) < 48:
            exclusions.append(
                {
                    "location": meta["location"],
                    "location_order": meta["location_order"],
                    "instrument_scan_number": meta["instrument_scan_number"],
                    "available_pixels": len(pixel_rows),
                    "reason": "fewer_than_48_saved_pixels_no_PB_fit",
                }
            )
            continue
        snap = np.asarray(
            [
                row["approach_snap_distance_nm"] for row in pixel_rows
                if row["approach_snap_detected"] and np.isfinite(row["approach_snap_distance_nm"])
            ],
            dtype=np.float64,
        )
        max_snap = float(np.max(snap)) if snap.size else 0.0
        safe_dmin = max(20.0, 5.0 * math.ceil((max_snap + 7.5) / 5.0))
        drained, diagnostics = drainage_corrected_pixels(
            meta, pixel_rows, distance, arrays["line_pN"][index]
        )
        drainage_rows.extend(diagnostics)
        summaries: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, float]] = {}
        for baseline, values in (
            ("linear_drift_corrected", arrays["line_pN"][index]),
            ("far_constant_referenced", arrays["constant_pN"][index]),
            ("linear_then_no_slip_corrected", drained),
        ):
            median = np.nanmedian(values, axis=0)
            spread = 1.4826 * np.nanmedian(np.abs(values - median), axis=0)
            counts = np.sum(np.isfinite(values), axis=0)
            far = finite(median[(distance >= 400.0) & (distance <= 500.0)])
            noise = max(3.0, float(base.robust_mad(far)))
            sigma = np.maximum(spread, noise)
            summaries[baseline] = median, sigma, counts, noise
            for j, d in enumerate(distance):
                inputs.append(
                    {
                        "acquisition_order": order,
                        "location": meta["location"],
                        "location_order": meta["location_order"],
                        "instrument_scan_number": meta["instrument_scan_number"],
                        "map_midpoint_time": meta["map_midpoint_time"],
                        "relative_time_min": meta["relative_time_min"],
                        "nominal_speed_um_per_s": meta["nominal_speed_um_per_s"],
                        "baseline_method": baseline,
                        "distance_nm": float(d),
                        "force_median_pN": float(median[j]),
                        "pixel_spread_mad_pN": float(spread[j]),
                        "fit_weight_scale_pN": float(sigma[j]),
                        "pixel_count": int(counts[j]),
                        "safe_dmin_nm": safe_dmin,
                        "detected_snap_max_nm": max_snap,
                    }
                )

        for variant in VARIANTS:
            lower = safe_dmin
            upper = 250.0
            baseline = "linear_drift_corrected"
            if variant == "constant_baseline":
                baseline = "far_constant_referenced"
            elif variant == "lower_plus_10nm":
                lower += 10.0
            elif variant == "upper_200nm":
                upper = 200.0
            elif variant == "common_150_250nm":
                lower = 150.0
            elif variant == "no_slip_hydrodynamic":
                baseline = "linear_then_no_slip_corrected"
            force, sigma, counts, noise = summaries[baseline]
            selected = (
                (distance >= lower) & (distance <= upper) & (counts >= 48)
                & np.isfinite(force) & np.isfinite(sigma) & (sigma > 0)
            )
            if np.count_nonzero(selected) < 12:
                exclusions.append(
                    {
                        "location": meta["location"],
                        "location_order": meta["location_order"],
                        "instrument_scan_number": meta["instrument_scan_number"],
                        "variant": variant,
                        "available_bins": int(np.count_nonzero(selected)),
                        "reason": "fewer_than_12_snap_safe_fit_bins",
                    }
                )
                continue
            if variant != "common_150_250nm" and np.any(distance[selected] - 2.5 < max_snap + 5.0 - 1e-12):
                raise RuntimeError("PB fit window violates the five-nm snap guard")
            tasks.append(
                {
                    "acquisition_order": order,
                    "location": meta["location"],
                    "location_order": meta["location_order"],
                    "instrument_scan_number": meta["instrument_scan_number"],
                    "source": meta["source"],
                    "map_start_time": meta["map_start_time"],
                    "map_end_time": meta["map_end_time"],
                    "map_midpoint_time": meta["map_midpoint_time"],
                    "relative_time_min": meta["relative_time_min"],
                    "nominal_speed_um_per_s": meta["nominal_speed_um_per_s"],
                    "measured_gap_speed_um_per_s": meta["gap_speed_20_200nm_median_um_per_s"],
                    "variant": variant,
                    "baseline_method": baseline,
                    "dmin_nm": float(np.min(distance[selected])),
                    "dmax_nm": float(np.max(distance[selected])),
                    "safe_dmin_nm": safe_dmin,
                    "detected_snap_max_nm": max_snap,
                    "detected_snap_median_nm": quantile(snap, 0.5),
                    "snap_detected_pixels": int(snap.size),
                    "minimum_pixels_per_fitted_bin": int(np.min(counts[selected])),
                    "noise_floor_pN": noise,
                    "hydro_speed_um_s": 0.0,
                    "no_slip_correction_applied_to_pixels": variant == "no_slip_hydrodynamic",
                    "distance_nm": distance[selected].tolist(),
                    "force_pN": force[selected].tolist(),
                    "sigma_pN": sigma[selected].tolist(),
                }
            )
    validation = {
        "upstream_maps": len(inventory),
        "PB_eligible_maps": len({task["acquisition_order"] for task in tasks}),
        "fit_tasks": len(tasks),
        "excluded_map_scans": sorted({row["instrument_scan_number"] for row in exclusions if row.get("reason") == "fewer_than_48_saved_pixels_no_PB_fit"}),
        "guard_rule": "5 nm bin lower edge is at least 5 nm beyond maximum detected snap; true common variant uses 150-250 nm",
    }
    return inventory, inputs, tasks, exclusions, drainage_rows, validation


def finite(values) -> np.ndarray:
    selected = np.asarray(values, dtype=np.float64)
    return selected[np.isfinite(selected)]


def quantile(values, q: float) -> float:
    selected = finite(values)
    return float(np.quantile(selected, q)) if selected.size else float("nan")


def run_fits(tasks: list[dict], workers: int) -> list[dict]:
    if workers <= 1:
        return [fit_task(task) for task in tasks]
    results: list[dict] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(fit_task, task): task for task in tasks}
        for index, future in enumerate(as_completed(futures), 1):
            results.append(future.result())
            if index % 24 == 0 or index == len(tasks):
                print(f"Completed {index}/{len(tasks)} apparent-PB fits", flush=True)
    return sorted(results, key=lambda row: (row["acquisition_order"], VARIANTS.index(row["variant"])))


def reuse_fits(tasks: list[dict]) -> list[dict]:
    """Reuse completed optimizer outputs only after checking task identity."""
    path = OUT / "fit_parameters_all_variants.csv"
    if not path.is_file():
        raise FileNotFoundError(f"No completed fit table to reuse: {path}")
    results = read_csv(path)
    expected = {(task["acquisition_order"], task["variant"]): task for task in tasks}
    observed = {(row["acquisition_order"], row["variant"]): row for row in results}
    if len(results) != len(tasks) or set(observed) != set(expected):
        raise RuntimeError("Saved optimizer result keys do not match current fit tasks")
    for key, task in expected.items():
        row = observed[key]
        for field in ("baseline_method", "source"):
            if row[field] != task[field]:
                raise RuntimeError(f"Saved optimizer result changed at {key} {field}")
        for field in ("dmin_nm", "dmax_nm", "safe_dmin_nm", "detected_snap_max_nm"):
            if not np.isclose(row[field], task[field], rtol=0.0, atol=1e-12):
                raise RuntimeError(f"Saved optimizer result changed at {key} {field}")
    print(f"Reusing {len(results)} completed fits after task-identity checks", flush=True)
    return sorted(results, key=lambda row: (row["acquisition_order"], VARIANTS.index(row["variant"])))


def time_parameter_models(inventory: list[dict], results: list[dict]) -> list[dict]:
    primary = [row for row in results if row["variant"] == "primary"]
    rows: list[dict] = []
    for location in ("A", "B"):
        location_primary = sorted(
            [row for row in primary if row["location"] == location],
            key=lambda row: row["location_order"],
        )
        full_span = max(row["relative_time_min"] for row in inventory if row["location"] == location)
        for metric, unit in (
            ("lambda_D_apparent_nm", "nm"),
            ("surface_potential_magnitude_apparent_mV", "mV"),
        ):
            for subset_name, subset in (
                ("all_eligible_maps", location_primary),
                ("quality_pass_maps", [row for row in location_primary if row["fit_quality_pass"]]),
            ):
                if len(subset) < 8:
                    rows.append(
                        {
                            "location": location,
                            "metric": metric,
                            "unit": unit,
                            "subset": subset_name,
                            "n_maps": len(subset),
                            "status": "insufficient_maps_for_time_speed_model",
                        }
                    )
                    continue
                times = np.asarray([row["relative_time_min"] for row in subset], dtype=np.float64)
                time_fraction = times / full_span
                speeds = np.asarray([row["measured_gap_speed_um_per_s"] for row in subset], dtype=np.float64)
                values = np.asarray([row[metric] for row in subset], dtype=np.float64)
                design = np.column_stack([np.ones(len(subset)), time_fraction, speeds - np.mean(speeds)])
                fit = time_analysis.ols_hc3(design, values, np.asarray([0.0, 1.0, 0.0]))
                rows.append(
                    {
                        "location": location,
                        "metric": metric,
                        "unit": unit,
                        "subset": subset_name,
                        "n_maps": len(subset),
                        "full_location_span_min": full_span,
                        "time_adjusted_full_span_change": fit["estimate"],
                        "HC3_se": fit["hc3_se"],
                        "HC3_95CI_low": fit["hc3_ci_low"],
                        "HC3_95CI_high": fit["hc3_ci_high"],
                        "HC3_p": fit["hc3_p"],
                        "HAC1_se": fit["hac1_se"],
                        "HAC1_95CI_low": fit["hac1_ci_low"],
                        "HAC1_95CI_high": fit["hac1_ci_high"],
                        "HAC1_p": fit["hac1_p"],
                        "HAC2_se": fit["hac2_se"],
                        "HAC2_95CI_low": fit["hac2_ci_low"],
                        "HAC2_95CI_high": fit["hac2_ci_high"],
                        "HAC2_p": fit["hac2_p"],
                        "gap_speed_coefficient_per_um_s": float(fit["beta"][-1]),
                        "r2": fit["r2"],
                        "rmse": fit["rmse"],
                        "design_rank": fit["rank"],
                        "design_parameters": fit["parameters"],
                        "design_condition_number": fit["condition_number"],
                        "residual_lag1_correlation": fit["residual_lag1_correlation"],
                        "status": "descriptive_model_conditioned_chronology",
                    }
                )
    return rows


def sensitivity_summary(results: list[dict]) -> list[dict]:
    primary = {(row["acquisition_order"]): row for row in results if row["variant"] == "primary"}
    rows: list[dict] = []
    for location in ("A", "B"):
        for variant in VARIANTS[1:]:
            related = [row for row in results if row["location"] == location and row["variant"] == variant and row["acquisition_order"] in primary]
            lambda_change = np.asarray([100.0 * (row["lambda_D_apparent_nm"] / primary[row["acquisition_order"]]["lambda_D_apparent_nm"] - 1.0) for row in related])
            potential_change = np.asarray([100.0 * (row["surface_potential_magnitude_apparent_mV"] / primary[row["acquisition_order"]]["surface_potential_magnitude_apparent_mV"] - 1.0) for row in related])
            rows.append(
                {
                    "location": location,
                    "variant": variant,
                    "maps_with_variant": len(related),
                    "quality_pass_maps": sum(row["fit_quality_pass"] for row in related),
                    "lambda_relative_change_median_percent": quantile(lambda_change, 0.5),
                    "lambda_relative_change_q25_percent": quantile(lambda_change, 0.25),
                    "lambda_relative_change_q75_percent": quantile(lambda_change, 0.75),
                    "potential_relative_change_median_percent": quantile(potential_change, 0.5),
                    "potential_relative_change_q25_percent": quantile(potential_change, 0.25),
                    "potential_relative_change_q75_percent": quantile(potential_change, 0.75),
                }
            )
    return rows


def prediction_rows(results: list[dict]) -> list[dict]:
    rows: list[dict] = []
    configure_model()
    for result in results:
        distance = np.arange(result["dmin_nm"], result["dmax_nm"] + 0.1, 5.0)
        parameters = prior_pb.result_parameters(result)
        prediction = prior_pb.predict(distance, parameters, 0.0)
        for d, value in zip(distance, prediction, strict=True):
            rows.append(
                {
                    "acquisition_order": result["acquisition_order"],
                    "location": result["location"],
                    "location_order": result["location_order"],
                    "instrument_scan_number": result["instrument_scan_number"],
                    "variant": result["variant"],
                    "distance_nm": float(d),
                    "predicted_force_pN": float(value),
                }
            )
    return rows


def plot_parameters(results: list[dict]) -> None:
    primary = [row for row in results if row["variant"] == "primary"]
    figure, axes = plt.subplots(2, 2, figsize=(13.5, 8.5), sharex="row")
    for row_index, location in enumerate(("A", "B")):
        subset = sorted([row for row in primary if row["location"] == location], key=lambda row: row["location_order"])
        for column, (field, label) in enumerate(
            (("lambda_D_apparent_nm", "Apparent λD (nm)"), ("surface_potential_magnitude_apparent_mV", "Apparent |ψ| (mV)"))
        ):
            ax = axes[row_index, column]
            times = np.asarray([row["relative_time_min"] for row in subset])
            values = np.asarray([row[field] for row in subset])
            ax.plot(times, values, color="0.5", lw=1.0)
            for row, x, y in zip(subset, times, values, strict=True):
                ax.scatter(
                    [x], [y], s=45, facecolors=time_analysis.COLORS[row["nominal_speed_um_per_s"]] if row["fit_quality_pass"] else "none",
                    edgecolors=time_analysis.COLORS[row["nominal_speed_um_per_s"]], linewidths=1.5,
                )
            ax.set_title(f"Location {location}")
            ax.set_xlabel("Elapsed time within location (min)")
            ax.set_ylabel(label)
            ax.grid(alpha=0.2)
    figure.suptitle("Finite-speed, model-conditioned apparent PB parameters\nfilled: quality screen pass; open: flagged")
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(FIGURES / "apparent_parameters_vs_time.png", dpi=220)
    plt.close(figure)


def plot_fits(inputs: list[dict], results: list[dict]) -> None:
    primary = [row for row in results if row["variant"] == "primary"]
    figure, axes = plt.subplots(7, 5, figsize=(18, 22), squeeze=False)
    configure_model()
    for ax, result in zip(axes.flat, primary, strict=False):
        selected = [
            row for row in inputs
            if row["acquisition_order"] == result["acquisition_order"]
            and row["baseline_method"] == "linear_drift_corrected"
            and result["dmin_nm"] <= row["distance_nm"] <= result["dmax_nm"]
            and row["pixel_count"] >= 48
        ]
        distance = np.asarray([row["distance_nm"] for row in selected])
        force = np.asarray([row["force_median_pN"] for row in selected])
        prediction = prior_pb.predict(distance, prior_pb.result_parameters(result), 0.0)
        ax.plot(distance, force, "o", ms=2.5, color="#264653", label="data")
        ax.plot(distance, prediction, lw=1.3, color="#e76f51", label="fit")
        flag = "PASS" if result["fit_quality_pass"] else result["quality_flags"]
        ax.set_title(f"{result['location']}{result['location_order']:02d} scan {result['instrument_scan_number']}\nλ={result['lambda_D_apparent_nm']:.2f} nm, ψ={result['surface_potential_magnitude_apparent_mV']:.1f} mV, {flag}", fontsize=8)
        ax.grid(alpha=0.15)
        ax.tick_params(labelsize=7)
    for ax in axes.flat[len(primary):]:
        ax.axis("off")
    figure.supxlabel("Separation D (nm)")
    figure.supylabel("Map-median force (pN)")
    figure.suptitle("All eligible primary apparent-PB fits", y=0.998)
    figure.tight_layout(rect=(0.025, 0.02, 1, 0.992))
    figure.savefig(FIGURES / "all_primary_fits.png", dpi=180)
    plt.close(figure)


def render_report(
    inventory: list[dict], results: list[dict], exclusions: list[dict], models: list[dict],
    sensitivity: list[dict], validation: dict,
) -> None:
    primary = [row for row in results if row["variant"] == "primary"]
    lines = [
        "# 10-09-26 secondary apparent-PB analysis",
        "",
        "## Direct result and scope",
        "",
        "These fits are secondary to the model-free force/history result. They are finite-speed, map-median, equal-constant-potential nonlinear-PB + sphere-plane vdW fits with a fitted constant residual baseline. `lambda_D` and `|psi|` are apparent/model-conditioned; they are not validated bulk Debye length or electrokinetic zeta potential.",
        "",
        f"PB fitting retained `{validation['PB_eligible_maps']}/36` primary maps. Scan 36 has only 44 saved pixels and was excluded rather than lowering the established >=48-pixel rule. `{validation['primary_quality_pass_maps']}/{len(primary)}` primary fits passed the descriptive quality screen.",
        "",
        "![Apparent parameters](figures/apparent_parameters_vs_time.png)",
        "",
        "## Time-adjusted parameter diagnostics",
        "",
        "The model is parameter = intercept + full-location time fraction + measured gap speed. HC3 is the primary map-level interval; Newey-West HAC(1/2) is an ordered-map serial-correlation sensitivity. Quality-pass-only rows are a sensitivity branch and can be selection-biased.",
        "",
        "| Location | parameter | subset | n | full-span change | HC3 95% CI | HAC(2) 95% CI | HC3 p | HAC(2) p |",
        "|---|---|---|---:|---:|:---|:---|---:|---:|",
    ]
    for row in models:
        if row["status"].startswith("insufficient"):
            lines.append(f"| {row['location']} | {row['metric']} | {row['subset']} | {row['n_maps']} | — | — | — | — | — |")
        else:
            lines.append(
                f"| {row['location']} | {row['metric']} | {row['subset']} | {row['n_maps']} | {row['time_adjusted_full_span_change']:+.3f} {row['unit']} | [{row['HC3_95CI_low']:+.3f}, {row['HC3_95CI_high']:+.3f}] | [{row['HAC2_95CI_low']:+.3f}, {row['HAC2_95CI_high']:+.3f}] | {row['HC3_p']:.4g} | {row['HAC2_p']:.4g} |"
            )
    lines += [
        "",
        "## Per-map primary results",
        "",
        "The lower bound is set per map so the lower edge of every fitted bin lies at least 5 nm beyond the largest detected approach snap. A single long-distance event in A01 therefore moves its lower bound to 150 nm; its result must be read from the quality flags, not treated like a near-field fit.",
        "",
        "| Location-map | scan | speed / µm/s | time / min | window / nm | apparent λD / nm | apparent |ψ| / mV | R2 | status |",
        "|---|---:|---:|---:|:---|---:|---:|---:|---|",
    ]
    for row in primary:
        status = "PASS" if row["fit_quality_pass"] else row["quality_flags"]
        lines.append(
            f"| {row['location']}{row['location_order']:02d} | {row['instrument_scan_number']} | {row['nominal_speed_um_per_s']:g} | {row['relative_time_min']:.2f} | {row['dmin_nm']:.0f}–{row['dmax_nm']:.0f} | {row['lambda_D_apparent_nm']:.3f} | {row['surface_potential_magnitude_apparent_mV']:.3f} | {row['r2']:.4f} | {status} |"
        )
    lines += [
        "",
        "![All fits](figures/all_primary_fits.png)",
        "",
        "## Method sensitivity",
        "",
        "Variant changes are paired to each map's primary optimum; quartiles across maps are descriptive, not confidence intervals.",
        "",
        "| Location | variant | maps | pass | median lambda change | median potential change |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in sensitivity:
        lines.append(
            f"| {row['location']} | {row['variant']} | {row['maps_with_variant']} | {row['quality_pass_maps']} | {row['lambda_relative_change_median_percent']:+.2f}% | {row['potential_relative_change_median_percent']:+.2f}% |"
        )
    lines += [
        "",
        f"The shared 150–250 nm tail produced `{validation['common_window_quality_pass_fits']}/35` quality-pass fits; it does not independently identify the PB parameters.",
        "",
        "## Model and evidence boundary",
        "",
        f"- Fixed inputs: R={RADIUS_M*1e6:.6f} µm, T={time_analysis.TEMPERATURE_C:.1f} °C, water epsilon_r={base.EPSILON_R[0]:g}, A_H={base.HAMAKER_J:.3g} J; force reconstruction uses D4 k={K_N_PER_M:.9f} N/m and the upstream current-water InvOLS.",
        "- Reported parameter trends condition on fixed k, InvOLS, R, temperature, dielectric constant, Hamaker constant and boundary condition; these systematic uncertainties are not included in HC3/HAC intervals.",
        "- The surface model assumes symmetric 1:1 electrolyte, equal constant potential, Derjaguin geometry and a constant residual baseline. The force fit determines potential magnitude only.",
        "- `no_slip_hydrodynamic` subtracts a fixed ideal 6*pi*eta*R^2*U/D sensitivity after applying the same per-curve far-line projection. It does not identify slip or hydrodynamic amplitude from the data.",
        "- The 150–250 nm common-window branch tests whether a shared tail supports the parameters. High R2 alone does not override boundary, Jacobian, weak-signal or local-scale flags.",
        "- Adjacent distance bins and pixels are correlated. Local parameter scales are conditioning diagnostics relative to spatial MAD, not experimental confidence intervals.",
        "- Any parameter chronology must agree with the primary measured-force and normalized-shape evidence and remains unable to distinguish surface relaxation, solution evolution, contact history and instrumental effects uniquely.",
        "",
        "## Numerical checks",
        "",
        f"- Optimizer success: `{validation['optimizer_success_fits']}/{validation['fit_tasks_completed']}`; finite parameter rows: `{validation['finite_parameter_fits']}/{validation['fit_tasks_completed']}`.",
        f"- Maximum prediction re-evaluation difference: `{validation['prediction_closure_max_abs_pN']:.3e} pN`.",
        f"- Maximum relative two-start cost difference among primary fits: `{validation['primary_max_two_start_relative_cost_difference']:.3e}`; minimum snap-guard margin: `{validation['minimum_snap_guard_margin_nm']:.3f} nm`.",
        f"- Ideal no-slip check F(2U)/F(U): `{validation['hydrodynamic_double_speed_ratio']:.12f}`; F(20 nm)/F(50 nm): `{validation['hydrodynamic_20_to_50nm_ratio']:.12f}`.",
        f"- Existing nonlinear-PB/base self-check: `{validation['existing_model_self_check']}`.",
        "- `provenance.json` binds upstream artifacts, raw identities, code and software; `artifact_manifest.sha256` covers this directory.",
    ]
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def create_manifest() -> None:
    paths = [path for path in OUT.rglob("*") if path.is_file() and path.name != "artifact_manifest.sha256"]
    lines = [f"{time_analysis.sha256_file(path)}  {path.relative_to(OUT).as_posix()}" for path in sorted(paths)]
    (OUT / "artifact_manifest.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--reuse-fit-results", action="store_true",
        help="reuse the saved optimizer table after exact task-window identity checks",
    )
    arguments = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    configure_model()

    inventory, inputs, tasks, exclusions, drainage_rows, validation = prepare()
    print(f"Prepared {len(tasks)} PB tasks from {validation['PB_eligible_maps']} eligible maps", flush=True)
    results = reuse_fits(tasks) if arguments.reuse_fit_results else run_fits(tasks, arguments.workers)
    models = time_parameter_models(inventory, results)
    sensitivity = sensitivity_summary(results)
    predictions = prediction_rows(results)

    prior_pb.refresh_diagnostics(results, inputs)
    time_analysis.write_csv(OUT / "map_median_fit_inputs.csv", inputs)
    time_analysis.write_csv(OUT / "fit_parameters_all_variants.csv", results)
    time_analysis.write_csv(OUT / "fit_parameters_per_map.csv", [row for row in results if row["variant"] == "primary"])
    time_analysis.write_csv(OUT / "fit_exclusions.csv", exclusions)
    time_analysis.write_csv(OUT / "no_slip_correction_per_pixel.csv", drainage_rows)
    time_analysis.write_csv(OUT / "fitted_force_curves.csv", predictions)
    time_analysis.write_csv(OUT / "parameter_time_adjusted_models.csv", models)
    time_analysis.write_csv(OUT / "parameter_sensitivity_summary.csv", sensitivity)
    plot_parameters(results)
    plot_fits(inputs, results)

    primary = [row for row in results if row["variant"] == "primary"]
    configure_model()
    closure = 0.0
    for row in results:
        distance = np.arange(row["dmin_nm"], row["dmax_nm"] + 0.1, 5.0)
        direct = prior_pb.predict(distance, prior_pb.result_parameters(row), 0.0)
        saved = np.asarray([
            item["predicted_force_pN"] for item in predictions
            if item["acquisition_order"] == row["acquisition_order"] and item["variant"] == row["variant"]
        ])
        if direct.shape != saved.shape:
            raise RuntimeError("Prediction shape closure failed")
        closure = max(closure, float(np.max(np.abs(direct - saved))))
    hydro_u = hydro_pn(np.asarray([20.0, 50.0]), 1.0)
    primary_start_cost_spreads = []
    for row in primary:
        costs = np.asarray([float(value) for value in row["two_start_costs"].split(";")])
        primary_start_cost_spreads.append(
            float(np.ptp(costs) / np.max(costs)) if np.max(costs) > 0 else 0.0
        )
    validation.update(
        {
            "fit_tasks_completed": len(results),
            "optimizer_success_fits": sum(row["optimizer_success"] for row in results),
            "finite_parameter_fits": sum(
                np.all(np.isfinite([row["lambda_D_apparent_nm"], row["surface_potential_magnitude_apparent_mV"], row["baseline_pN"]]))
                for row in results
            ),
            "primary_quality_pass_maps": sum(row["fit_quality_pass"] for row in primary),
            "primary_flagged_maps": [
                {"scan": row["instrument_scan_number"], "flags": row["quality_flags"]}
                for row in primary if not row["fit_quality_pass"]
            ],
            "primary_max_two_start_relative_cost_difference": max(
                primary_start_cost_spreads
            ),
            "common_window_quality_pass_fits": sum(
                row["fit_quality_pass"]
                for row in results if row["variant"] == "common_150_250nm"
            ),
            "minimum_snap_guard_margin_nm": min(
                (row["dmin_nm"] - 2.5) - (row["detected_snap_max_nm"] + 5.0)
                for row in results
            ),
            "prediction_closure_max_abs_pN": closure,
            "hydrodynamic_double_speed_ratio": float(hydro_pn(np.asarray([20.0]), 2.0)[0] / hydro_pn(np.asarray([20.0]), 1.0)[0]),
            "hydrodynamic_20_to_50nm_ratio": float(hydro_u[0] / hydro_u[1]),
            "hydrodynamic_force_pN_at_1um_s": {"20nm": float(hydro_u[0]), "50nm": float(hydro_u[1])},
            "existing_model_self_check": base.self_test(),
        }
    )
    if validation["optimizer_success_fits"] != len(results) or validation["finite_parameter_fits"] != len(results):
        raise RuntimeError("At least one PB fit failed or returned nonfinite parameters")
    if closure > 1e-10:
        raise RuntimeError("Saved PB predictions do not close")
    if not np.isclose(validation["hydrodynamic_double_speed_ratio"], 2.0, rtol=1e-12):
        raise RuntimeError("Hydrodynamic speed scaling failed")
    if not np.isclose(validation["hydrodynamic_20_to_50nm_ratio"], 2.5, rtol=1e-12):
        raise RuntimeError("Hydrodynamic inverse-distance scaling failed")

    upstream_files = [
        DATA / name for name in (
            "all_pixel_force_curves.npz", "map_inventory_QC.csv", "pixel_QC.csv",
            "provenance.json", "artifact_manifest.sha256",
        )
    ]
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": "secondary finite-speed apparent nonlinear-PB fits for 10-09-26 two-location time series",
        "claim_scope": "apparent/model-conditioned lambda_D and potential magnitude; not bulk Debye length or zeta potential",
        "model": "equal constant-potential sphere-plane nonlinear 1:1 PB Derjaguin + sphere-plane vdW + fitted constant baseline",
        "potential_sign": "unidentified magnitude only",
        "probe_radius_m": RADIUS_M,
        "temperature_C": time_analysis.TEMPERATURE_C,
        "epsilon_r": base.EPSILON_R[0],
        "Hamaker_J": base.HAMAKER_J,
        "spring_constant_N_per_m": K_N_PER_M,
        "water_InvOLS_nm_per_V": inventory[0]["global_water_InvOLS_used_nm_per_V"],
        "aggregation": "median of available individual 5 nm pixel-bin medians; scan 36 excluded for fewer than 48 pixels",
        "weight": "spatial 1.4826*MAD floored by max(3 pN, robust MAD of map median at 400-500 nm); no sqrt(pixel count) reduction",
        "time_model": "map-level OLS with measured gap speed; HC3 primary and Newey-West HAC(1/2) ordered-map sensitivities",
        "fit_variants": list(VARIANTS),
        "fit_bounds": {
            "lambda_D_nm": [1.0, 1000.0],
            "potential_magnitude_mV": [0.1, 250.0],
            "baseline_pN": [-500.0, 500.0],
        },
        "snap_guard": validation["guard_rule"],
        "no_slip_sensitivity": "fixed ideal 6*pi*eta*R^2*U/D using per-pixel measured gap speed and far-line projection; not fitted or validated hydrodynamic amplitude",
        "validation": validation,
        "upstream_sha256": {path.relative_to(ROOT).as_posix(): time_analysis.sha256_file(path) for path in upstream_files},
        "code_sha256": {
            path.relative_to(ROOT).as_posix(): time_analysis.sha256_file(path)
            for path in (Path(__file__).resolve(), Path(prior_pb.__file__).resolve(), Path(time_analysis.__file__).resolve(), Path(base.__file__).resolve())
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
    }
    (OUT / "provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2, default=json_default) + "\n",
        encoding="utf-8",
    )
    render_report(inventory, results, exclusions, models, sensitivity, validation)
    create_manifest()
    print(json.dumps(validation, indent=2, default=json_default), flush=True)
    print(f"Wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()

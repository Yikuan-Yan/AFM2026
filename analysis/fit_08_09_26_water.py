#!/usr/bin/env python3
"""Fit 08-09-26 map medians with the repository's nonlinear PB evaluator.

These are finite-speed, equal-potential model parameters. Each fit uses one
map, not a fabricated three-speed equilibrium intercept. Original calibration
and force reconstruction artifacts are read-only inputs.
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
from matplotlib.lines import Line2D
import numpy as np
import scipy
from scipy.optimize import least_squares

import fit_glycerol_surface_forces as base
import analyze_08_09_26_fixed_pixel as chronology


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "analysis" / "water_08_09_26_results"
OUT = DATA / "surface_fit"
MODEL = "nonlinear_pb_derjaguin"
RADIUS_M = 4.54685e-6  # User explicitly confirmed this radius in this task.
LOWER = np.array([math.log(1.0), math.log(0.1), -500.0])
UPPER = np.array([math.log(1000.0), math.log(250.0), 500.0])
VARIANTS = ("primary", "constant_baseline", "lower_plus_10nm",
            "upper_200nm", "common_80_250nm", "no_slip_hydrodynamic")


def read_csv(path: Path) -> list[dict]:
    return chronology.read_derived_csv(path)


def configure_model() -> None:
    base.PROBE_RADIUS_M = RADIUS_M


def hydro_pN(distance_nm: np.ndarray, speed_um_s: float) -> np.ndarray:
    """Positive no-slip sphere-plane drainage force on approach, for sensitivity."""
    eta = base.cheng_viscosity_mPa_s(0.0, base.TEMPERATURE_C) * 1e-3
    return 6 * np.pi * eta * RADIUS_M**2 * speed_um_s * 1e-6 / (distance_nm * 1e-9) * 1e12


def predict(distance: np.ndarray, parameters: np.ndarray,
            speed_um_s: float = 0.0) -> np.ndarray:
    return (base.total_equilibrium_force_pN(distance, math.exp(parameters[0]),
            math.exp(parameters[1]), parameters[2], base.EPSILON_R[0], MODEL)
            + hydro_pN(distance, speed_um_s))


def drainage_corrected_pixels(meta: dict, pixel_rows: list[dict], distance: np.ndarray,
                              force_matrix: np.ndarray) -> tuple[np.ndarray, list[dict]]:
    """Approximate the same far-line operator on the theoretical drainage force.

    The 1/D term uses each curve's measured gap speed. Its far-field line is
    fitted on that curve's original far samples against scanner height. The
    line is evaluated using h-h0 = D-F/k on the saved distance bins. Subtract
    per pixel before the spatial median because medians are nonlinear.
    """
    source = base.load_source(ROOT / meta["source"], 0)
    if source.sha256 != meta["sha256"] or len(source.curves) != 64 or source.skipped_curves:
        raise AssertionError("Raw drainage input identity differs from reconstruction")
    corrected = np.full_like(force_matrix, np.nan)
    diagnostics = []
    lookup = {row["point_index"]: row for row in pixel_rows}
    for curve in source.curves:
        point = int(curve.point_index)
        row = lookup[point]
        height = curve.measured_height_m
        far = base.fit_far_field_drift(height, curve.deflection_V)
        voltage = curve.deflection_V - base.baseline_voltage(height, curve.deflection_V, far)
        sensitivity = float(row["sensitivity_used_nm_per_V"]) * 1e-9
        h0 = float(row["contact_height_um"]) * 1e-6
        raw_distance = (height + sensitivity * voltage - h0) * 1e9
        speed = float(row["gap_speed_20_200nm_um_per_s"])
        if not np.isfinite(speed) or speed <= 0:
            raise AssertionError("Missing measured gap speed for drainage sensitivity")
        far_distance = raw_distance[:far.n_points]
        if np.any(far_distance <= 0):
            raise AssertionError("Nonpositive far separation")
        x = (height[:far.n_points] - h0) * 1e9
        slope, intercept, _, _ = base.robust_line(x, hydro_pN(far_distance, speed))
        bin_height = distance - force_matrix[point] / (0.2736 * 1000)
        subtracted = hydro_pN(distance, speed) - (slope * bin_height + intercept)
        corrected[point] = force_matrix[point] - subtracted
        diagnostics.append({"acquisition_order": meta["acquisition_order"], "point_index": point,
            "gap_speed_um_per_s": speed, "far_distance_min_nm": float(np.min(far_distance)),
            "far_distance_max_nm": float(np.max(far_distance)), "far_sample_count": far.n_points,
            "hydro_far_line_slope_pN_per_nm": slope, "hydro_far_line_intercept_pN": intercept,
            "hydro_after_far_line_at_80nm_pN": float(subtracted[np.flatnonzero(distance == 80)[0]])})
    return corrected, diagnostics


def local_geometry(residual, parameters: np.ndarray, n_points: int) -> tuple:
    # Central differences of the unmodified weighted residual diagnose the
    # local model, independently of soft_l1's transformed optimizer Jacobian.
    steps = np.array([1e-4, 1e-4, 0.01])
    jacobian = np.column_stack([
        (residual(parameters + np.eye(3)[j] * steps[j]) -
         residual(parameters - np.eye(3)[j] * steps[j])) / (2 * steps[j])
        for j in range(3)])
    norms = np.linalg.norm(jacobian, axis=0)
    rank, condition, correlation = 0, float("inf"), float("nan")
    se = np.full(3, np.nan)
    if np.all(np.isfinite(norms)) and np.all(norms > 0):
        _, singular, vt = np.linalg.svd(jacobian / norms, full_matrices=False)
        rank = int(np.sum(singular > max(jacobian.shape) * np.finfo(float).eps * singular[0]))
        condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else float("inf")
        if rank == 3:
            # Invert in normalized coordinates so an extremely small absolute
            # derivative cannot be silently dropped by pinv and look certain.
            covariance = ((vt.T / singular**2) @ vt) / np.outer(norms, norms)
            # Spatial heterogeneity is the stated weight scale, not a Gaussian
            # noise variance to be reduced when a smooth median fits closely.
            covariance *= max(1.0, float(np.sum(residual(parameters)**2)) / (n_points - 3))
            se = np.sqrt(np.maximum(0, np.diag(covariance)))
            if se[0] > 0 and se[1] > 0:
                correlation = float(covariance[0, 1] / (se[0] * se[1]))
    return rank, condition, correlation, se


def fit_one(task: dict) -> dict:
    configure_model()
    distance = np.asarray(task["distance_nm"], dtype=float)
    force = np.asarray(task["force_pN"], dtype=float)
    sigma = np.asarray(task["sigma_pN"], dtype=float)
    speed = float(task["hydro_speed_um_s"])

    def residual(parameters: np.ndarray) -> np.ndarray:
        return (predict(distance, parameters, speed) - force) / sigma

    starts = []
    for lam, psi in ((20.0, 35.0), (100.0, 100.0)):
        baseline = float(np.clip(np.median(force[-5:] - hydro_pN(distance[-5:], speed)), -499, 499))
        result = least_squares(residual, [math.log(lam), math.log(psi), baseline],
                               bounds=(LOWER, UPPER), method="trf", loss="soft_l1",
                               f_scale=1.0, xtol=1e-8, ftol=1e-8, gtol=1e-8,
                               max_nfev=250)
        starts.append(result)
    best = min(starts, key=lambda result: result.cost)
    parameters = best.x
    prediction = predict(distance, parameters, speed)
    raw_residual = force - prediction
    sse = float(np.sum(raw_residual**2))
    sst = float(np.sum((force - np.mean(force))**2))
    r2 = 1 - sse / sst if sst > 0 else float("nan")
    lam, psi = np.exp(parameters[:2])
    rank, condition, correlation, se = local_geometry(residual, parameters, len(distance))
    boundary_fraction = np.minimum((parameters - LOWER) / (UPPER - LOWER),
                                   (UPPER - parameters) / (UPPER - LOWER))
    reasons = []
    if len(distance) < 12:
        reasons.append("fewer_than_12_bins")
    if not best.success:
        reasons.append("optimizer")
    if rank != 3 or condition > 1e6:
        reasons.append("jacobian")
    if np.any(boundary_fraction < 1e-4):
        reasons.append("parameter_boundary")
    if not np.isfinite(r2) or r2 < 0.70:
        reasons.append("r2_below_0.70")
    signal = float(np.max(force) - np.median(force[-max(3, len(force)//5):]))
    if signal < 5 * float(task["noise_floor_pN"]):
        reasons.append("weak_signal")
    # This is a conditioning warning, not a confidence statement about pixels.
    if not np.all(np.isfinite(se[:2])) or np.any(se[:2] > 0.5):
        reasons.append("broad_local_parameter_scale")
    result = {key: value for key, value in task.items()
              if key not in ("distance_nm", "force_pN", "sigma_pN")}
    result.update({
        "lambda_D_apparent_nm": float(lam),
        "surface_potential_magnitude_apparent_mV": float(psi),
        "baseline_pN": float(parameters[2]), "r2": r2,
        "residual_rmse_pN": math.sqrt(sse / len(distance)),
        "weighted_sse": float(np.sum(residual(parameters)**2)),
        "soft_l1_cost": float(best.cost), "optimizer_success": bool(best.success),
        "optimizer_nfev": int(best.nfev), "jacobian_rank": rank,
        "normalized_jacobian_condition": condition,
        "local_log_lambda_scale": float(se[0]), "local_log_potential_scale": float(se[1]),
        "local_lambda_potential_correlation": correlation,
        "local_scale_interpretation": "spatial_weight_curvature_no_downscaling_not_experimental_CI",
        "fit_quality_pass": not reasons, "quality_flags": ";".join(reasons),
        "n_points": len(distance), "force_signal_pN": signal,
        "two_start_costs": ";".join(f"{item.cost:.9g}" for item in starts),
        "data_speed_groups": 1, "zero_speed_extrapolated": False,
        "model": MODEL, "probe_radius_um": RADIUS_M * 1e6,
        "temperature_model_C": base.TEMPERATURE_C,
    })
    return result


def prepare() -> tuple[list[dict], list[dict], list[dict], dict]:
    inventory = read_csv(DATA / "map_inventory_QC.csv")
    pixels = read_csv(DATA / "pixel_QC.csv")
    curves_csv = read_csv(DATA / "map_force_curves.csv")
    with np.load(DATA / "all_pixel_force_curves.npz") as archive:
        data = {key: archive[key] for key in archive.files}
    distance = data["distance_nm"]
    if not (np.array_equal(data["acquisition_order"], np.arange(1, 37)) and
            np.array_equal(data["point_index"], np.arange(64)) and
            data["line_pN"].shape == (36, 64, len(distance))):
        raise AssertionError("Unexpected map/pixel axes")
    if any(float(row["spring_constant_used_N_per_m"]) != 0.2736 for row in inventory):
        raise AssertionError("Cantilever-2 calibration not present")
    inputs, tasks, drainage_diagnostics = [], [], []
    mismatch = 0.0
    lookup = {(int(r["acquisition_order"]), r["baseline_method"], float(r["distance_nm"])): r
              for r in curves_csv}
    for index, meta in enumerate(inventory):
        order = index + 1
        pixel_rows = [row for row in pixels if row["acquisition_order"] == order]
        snap_distances = [float(row["approach_snap_distance_nm"]) for row in pixel_rows
                          if row["approach_snap_detected"]
                          and np.isfinite(float(row["approach_snap_distance_nm"]))]
        max_snap = max(snap_distances, default=0.0)
        dmin = max(20.0, 5 * math.ceil((max_snap + 7.5) / 5))
        start, stop = (datetime.fromisoformat(meta[k]) for k in ("map_start_time", "map_end_time"))
        midpoint = (start + (stop - start) / 2).isoformat()
        summaries = {}
        drained, diagnostics = drainage_corrected_pixels(meta, pixel_rows, distance, data["line_pN"][index])
        drainage_diagnostics.extend(diagnostics)
        for baseline, values in (("linear_drift_corrected", data["line_pN"][index]),
                                 ("far_constant_referenced", data["constant_pN"][index]),
                                 ("linear_then_no_slip_corrected", drained)):
            median = np.nanmedian(values, axis=0)
            spread = 1.4826 * np.nanmedian(np.abs(values - median), axis=0)
            counts = np.sum(np.isfinite(values), axis=0)
            far = median[(distance >= 400) & (distance <= 500)]
            noise = max(3.0, float(base.robust_mad(far)))
            sigma = np.maximum(spread, noise)
            for j, d in enumerate(distance):
                if baseline != "linear_then_no_slip_corrected":
                    original = lookup[(order, baseline, float(d))]
                    if np.isfinite(median[j]) and np.isfinite(float(original["force_median_pN"])):
                        mismatch = max(mismatch, abs(median[j] - float(original["force_median_pN"])))
                    elif np.isfinite(median[j]) != np.isfinite(float(original["force_median_pN"])):
                        raise AssertionError("Reconstructed median finite mask mismatch")
                    if int(original["pixel_count"]) != int(counts[j]):
                        raise AssertionError("Pixel counts disagree with upstream output")
                inputs.append({"acquisition_order": order, "map_midpoint_time": midpoint,
                    "nominal_speed_um_per_s": meta["nominal_speed_um_per_s"],
                    "baseline_method": baseline, "distance_nm": float(d),
                    "force_median_pN": float(median[j]), "pixel_spread_mad_pN": float(spread[j]),
                    "fit_weight_scale_pN": float(sigma[j]), "pixel_count": int(counts[j]),
                    "safe_dmin_nm": dmin, "detected_snap_max_nm": max_snap,
                    "in_primary_window": bool(dmin <= d <= 250 and counts[j] >= 48)})
            summaries[baseline] = (median, sigma, counts, noise)
        for variant in VARIANTS:
            lower = dmin + 10 if variant == "lower_plus_10nm" else dmin
            if variant == "common_80_250nm":
                lower = max(80, dmin)
            upper = 200.0 if variant == "upper_200nm" else 250.0
            baseline = "far_constant_referenced" if variant == "constant_baseline" else "linear_drift_corrected"
            if variant == "no_slip_hydrodynamic":
                baseline = "linear_then_no_slip_corrected"
            force, sigma, counts, noise = summaries[baseline]
            selected = (distance >= lower) & (distance <= upper) & (counts >= 48) & np.isfinite(force) & np.isfinite(sigma)
            if np.sum(selected) < 12:
                raise AssertionError(f"Map {order} has insufficient safe bins")
            if np.any(distance[selected] - 2.5 <= max_snap + 5):
                # Equality is allowed at the planned 5 nm guard edge.
                if np.any(distance[selected] - 2.5 < max_snap + 5):
                    raise AssertionError("Unsafe snap-in overlap")
            tasks.append({"acquisition_order": order, "block": meta["block"],
                "source": meta["source"], "map_start_time": meta["map_start_time"],
                "map_end_time": meta["map_end_time"], "map_midpoint_time": midpoint,
                "nominal_speed_um_per_s": float(meta["nominal_speed_um_per_s"]),
                "variant": variant, "baseline_method": baseline,
                "dmin_nm": float(lower), "dmax_nm": upper,
                "safe_dmin_nm": dmin, "detected_snap_max_nm": max_snap,
                "detected_snap_median_nm": float(np.median(snap_distances)),
                "detected_snap_p95_nm": float(np.percentile(snap_distances, 95)),
                "snap_detected_pixels": len(snap_distances),
                "minimum_pixels_per_fitted_bin": int(np.min(counts[selected])),
                "noise_floor_pN": noise,
                "hydro_speed_um_s": 0.0,
                "no_slip_correction_applied_to_pixels": variant == "no_slip_hydrodynamic",
                "distance_nm": distance[selected].tolist(), "force_pN": force[selected].tolist(),
                "sigma_pN": sigma[selected].tolist()})
        if order % 6 == 0:
            print(f"Prepared maps 1–{order}; raw hashes and far-baseline drainage operator checked", flush=True)
    if mismatch > 1e-9:
        raise AssertionError(f"Upstream medians mismatch: {mismatch}")
    base.write_csv(OUT / "no_slip_correction_per_pixel.csv", drainage_diagnostics)
    return inventory, inputs, tasks, {"map_median_max_difference_from_upstream_pN": mismatch,
        "maps": 36, "approach_curves": 2304, "fit_tasks": len(tasks),
        "guard_rule": "bin lower edge at least 5 nm beyond largest detected snap per map"}


def result_parameters(row: dict) -> np.ndarray:
    return np.array([math.log(float(row["lambda_D_apparent_nm"])),
                     math.log(float(row["surface_potential_magnitude_apparent_mV"])),
                     float(row["baseline_pN"])])


def refresh_diagnostics(results: list[dict], inputs: list[dict]) -> None:
    """Re-evaluate geometry/flags at saved optima without refitting the data."""
    for index, row in enumerate(results):
        selected = [r for r in inputs if r["acquisition_order"] == row["acquisition_order"]
                    and r["baseline_method"] == row["baseline_method"]
                    and row["dmin_nm"] <= r["distance_nm"] <= row["dmax_nm"]
                    and r["pixel_count"] >= 48 and np.isfinite(r["force_median_pN"])]
        if len(selected) != row["n_points"]:
            raise AssertionError("Saved fit input count changed")
        d = np.array([r["distance_nm"] for r in selected])
        y = np.array([r["force_median_pN"] for r in selected])
        sigma = np.array([r["fit_weight_scale_pN"] for r in selected])
        parameters = result_parameters(row)
        residual = lambda p: (predict(d, p, row["hydro_speed_um_s"]) - y) / sigma
        if not np.isclose(np.sum(residual(parameters)**2), row["weighted_sse"], rtol=1e-8, atol=1e-10):
            raise AssertionError("Saved optimum no longer matches its weighted residual")
        rank, condition, correlation, se = local_geometry(residual, parameters, len(d))
        flags = [flag for flag in row["quality_flags"].split(";")
                 if flag and flag not in ("jacobian", "broad_local_parameter_scale")]
        if rank != 3 or condition > 1e6:
            flags.append("jacobian")
        if not np.all(np.isfinite(se[:2])) or np.any(se[:2] > 0.5):
            flags.append("broad_local_parameter_scale")
        row.update({"jacobian_rank": rank, "normalized_jacobian_condition": condition,
                    "local_lambda_potential_correlation": correlation,
                    "local_log_lambda_scale": float(se[0]), "local_log_potential_scale": float(se[1]),
                    "local_scale_interpretation": "spatial_weight_curvature_no_downscaling_not_experimental_CI",
                    "fit_quality_pass": not flags, "quality_flags": ";".join(flags)})
        if (index + 1) % 54 == 0:
            print(f"Refreshed diagnostics {index+1}/{len(results)} at unchanged optima", flush=True)


def finish(inventory: list[dict], inputs: list[dict], results: list[dict], validation: dict) -> None:
    configure_model()
    primary = [row for row in results if row["variant"] == "primary"]
    by_key = {(row["acquisition_order"], row["variant"]): row for row in results}
    sensitivity = []
    for row in primary:
        related = [by_key[(row["acquisition_order"], variant)] for variant in VARIANTS[1:]]
        for other in related:
            sensitivity.append({"acquisition_order": row["acquisition_order"],
                "variant": other["variant"], "fit_quality_pass": other["fit_quality_pass"],
                "lambda_relative_change_percent": 100 * (other["lambda_D_apparent_nm"] / row["lambda_D_apparent_nm"] - 1),
                "potential_relative_change_percent": 100 * (other["surface_potential_magnitude_apparent_mV"] / row["surface_potential_magnitude_apparent_mV"] - 1)})
    base.write_csv(OUT / "fit_parameters_all_variants.csv", results)
    base.write_csv(OUT / "fit_parameters_per_map.csv", primary)
    base.write_csv(OUT / "parameter_sensitivity.csv", sensitivity)
    prediction_rows = []
    for row in results:
        for d, y in zip(np.arange(row["dmin_nm"], row["dmax_nm"] + 1, 5),
                        predict(np.arange(row["dmin_nm"], row["dmax_nm"] + 1, 5), result_parameters(row), row["hydro_speed_um_s"])):
            prediction_rows.append({"acquisition_order": row["acquisition_order"],
                "variant": row["variant"], "baseline_method": row["baseline_method"],
                "no_slip_correction_applied_to_pixels": row["no_slip_correction_applied_to_pixels"],
                "distance_nm": float(d), "predicted_force_pN": float(y)})
    base.write_csv(OUT / "fitted_force_curves.csv", prediction_rows)
    figures(inventory, inputs, primary, by_key)
    validation["existing_model_self_checks"] = base.self_test()
    validation["primary_quality_pass_maps"] = [r["acquisition_order"] for r in primary if r["fit_quality_pass"]]
    validation["primary_flagged_maps"] = [{"map": r["acquisition_order"], "flags": r["quality_flags"]} for r in primary if not r["fit_quality_pass"]]
    input_files = [DATA / name for name in ("all_pixel_force_curves.npz", "map_inventory_QC.csv", "pixel_QC.csv", "map_force_curves.csv", "provenance.json")]
    provenance = {"created_utc": datetime.now(timezone.utc).isoformat(),
        "method": "single-map median force fit; finite approach speed; equal constant potential nonlinear PB + fixed van der Waals + constant fitted offset",
        "potential_sign": "unidentified; magnitude only; not an electrokinetic zeta measurement",
        "aggregation": "median of 64 individual 5 nm curve-bin medians; no averaging of fitted pixel parameters",
        "weight": "spatial 1.4826*MAD, floored by max(3 pN, MAD of map median at 400-500 nm); no sqrt(64) reduction",
        "local_parameter_scale": "inverse weighted model information from column-normalized SVD; multiplied by max(reduced weighted SSE,1); heterogeneity-relative diagnostic, not an experimental confidence interval",
        "parameter_bounds": {"lambda_nm": [1, 1000], "potential_magnitude_mV": [0.1, 250], "offset_pN": [-500, 500]},
        "radius_m": RADIUS_M, "radius_source": "user explicitly requested reuse of R=4.54685 um; not newly calibrated",
        "spring_constant_N_per_m": 0.2736, "spring_constant_supplied_uncertainty_N_per_m": 0.0063,
        "temperature_C": base.TEMPERATURE_C, "temperature_source": "repository model assumption; not a new map temperature measurement",
        "epsilon_r": base.EPSILON_R[0], "Hamaker_J": base.HAMAKER_J,
        "no_slip_sensitivity": "subtracts 6*pi*eta*R^2*measured per-pixel gap speed/D after applying an approximate matching far-line projection; per-pixel correction before median; chi fixed to 1; not experimentally verified",
        "validation": validation,
        "input_sha256": {p.relative_to(ROOT).as_posix(): base.sha256_file(p) for p in input_files},
        "code_sha256": {p.relative_to(ROOT).as_posix(): base.sha256_file(p) for p in (Path(__file__), Path(base.__file__), Path(chronology.__file__))},
        "software": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__}}
    (OUT / "provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
    write_report(primary, results, validation)
    manifest = [f"{base.sha256_file(p)}  {p.relative_to(OUT).as_posix()}" for p in sorted(OUT.rglob("*"))
                if p.is_file() and p.name != "artifact_manifest.sha256"]
    (OUT / "artifact_manifest.sha256").write_text("\n".join(manifest) + "\n", encoding="utf-8")


def figures(inventory: list[dict], inputs: list[dict], primary: list[dict], by_key: dict) -> None:
    out = OUT / "figures"
    out.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 170, "svg.fonttype": "none"})
    times = chronology.measurement_times(inventory, "map_midpoint_time")
    keys = ("lambda_D_apparent_nm", "surface_potential_magnitude_apparent_mV")
    labels = (r"Apparent Debye length $\lambda_D$ (nm)", r"Apparent surface potential $|\psi_s|$ (mV)")
    fig, axes = plt.subplots(2, 1, figsize=(13.2, 8.1), sharex=True, layout="constrained")
    for ax, key, label in zip(axes, keys, labels):
        values = np.array([row[key] for row in primary])
        good = np.array([bool(row["fit_quality_pass"]) for row in primary])
        unusable = np.array([any(flag in row["quality_flags"] for flag in
                           ("parameter_boundary", "r2_below", "optimizer", "jacobian", "weak_signal"))
                           for row in primary])
        weak = ~good & ~unusable
        ax.plot(times, np.where(good, values, np.nan), color="0.65", lw=0.9, zorder=1)
        for speed, color in chronology.COLORS.items():
            idx = np.array([r["nominal_speed_um_per_s"] == speed for r in primary])
            ax.scatter(times[idx & good], values[idx & good], color=color, edgecolor="0.25", lw=0.5, s=42, zorder=3)
            ax.scatter(times[idx & weak], values[idx & weak], facecolors="none", edgecolors=color, lw=1.8, s=52, zorder=4)
        ax.scatter(times[unusable], np.full(np.sum(unusable), 0.035), transform=ax.get_xaxis_transform(),
                   marker="x", color="#ad3f36", s=43, lw=1.4, zorder=4)
        for index in np.flatnonzero(unusable):
            ax.annotate(f"#{index+1}", (times[index], 0.035), xycoords=ax.get_xaxis_transform(),
                        xytext=(3, 7), textcoords="offset points", fontsize=8, color="#ad3f36")
        for number in (1, 6, 12, 18, 19, 25, 30, 36):
            if unusable[number-1]:
                continue
            ax.annotate(f"#{number}", (times[number-1], values[number-1]), xytext=(2, 9), textcoords="offset points", fontsize=8, color="0.3")
        chronology.axes_chronology(ax, inventory, "map_midpoint_time")
        ax.set_ylabel(label)
        ax.set_ylim(0, 1.22 * np.max(values[~unusable]))
        ax.grid(axis="y", alpha=0.2)
    legend = [Line2D([], [], marker="o", color="none", markerfacecolor=c, markeredgecolor="0.3", label=f"{u:g} µm/s") for u, c in chronology.COLORS.items()]
    legend.append(Line2D([], [], marker="o", color="none", markeredgecolor="0.4", label="Weakly constrained optimum"))
    legend.append(Line2D([], [], marker="x", color="none", markeredgecolor="#ad3f36", label="Unusable fit (bottom rail)"))
    axes[0].legend(handles=legend, loc="upper right", frameon=False, ncol=2, fontsize=9)
    axes[0].set_title("08 Sep 2026 water | nonlinear PB fits to each map's median force", loc="left", fontweight="bold", pad=18)
    axes[1].set_xlabel("Acquisition time (instrument UTC+02:00; map midpoint)")
    axes[0].set_xlabel("")
    fig.get_layout_engine().set(rect=(0, 0.075, 1, 0.925))
    fig.text(0.07, 0.017, "64 curves/map · k = 0.2736 N/m · R = 4.54685 µm · safe lower limit varies by map, upper limit 250 nm\nFinite-speed parameters. Open circles show uncertain optima; bottom crosses are status markers, not measured parameter values.", fontsize=9, color="0.3")
    for suffix in ("png", "svg"):
        fig.savefig(out / f"debye_length_surface_potential_vs_time.{suffix}")
    plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.3), sharex=True, layout="constrained")
    comparisons = (("constant_baseline", "Constant-baseline sensitivity"),
                   ("lower_plus_10nm", "Raise lower fit limit by 10 nm"),
                   ("no_slip_hydrodynamic", "No-slip drainage sensitivity"))
    for col, (variant, title) in enumerate(comparisons):
        other = [by_key[(r["acquisition_order"], variant)] for r in primary]
        for row, (key, label) in enumerate(zip(keys, labels)):
            ax = axes[row, col]
            ax.plot(times, [r[key] for r in primary], "o-", color="#263b53", ms=3, lw=1, label="Primary")
            ax.plot(times, [r[key] for r in other], "s--", color="#d07835", ms=3, lw=1, label=variant.replace("_", " "))
            primary_flagged = np.array([not r["fit_quality_pass"] for r in primary])
            ax.scatter(times[primary_flagged], np.array([r[key] for r in primary])[primary_flagged],
                       marker="x", color="#a52a2a", s=28, zorder=4)
            flagged = np.array([not r["fit_quality_pass"] for r in other])
            ax.scatter(times[flagged], np.array([r[key] for r in other])[flagged], marker="x", color="#a52a2a", s=28, zorder=4)
            chronology.axes_chronology(ax, inventory, "map_midpoint_time")
            ax.set_ylim(bottom=0)
            ax.grid(axis="y", alpha=0.2)
            if col == 0:
                ax.set_ylabel(label)
            if row == 0:
                ax.set_title(title, fontsize=11)
            else:
                ax.set_xlabel("Instrument time (UTC+02:00)")
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle("Sensitivity of optimizer outputs | x = failed residual / heterogeneity-scale check; flagged values are not reliable estimates", fontsize=12)
    fig.savefig(out / "parameter_sensitivity_vs_time.png")
    plt.close(fig)
    fig, axes = plt.subplots(6, 6, figsize=(17, 16), sharex=True, layout="constrained")
    for index, (ax, fit) in enumerate(zip(axes.flat, primary)):
        curve = [r for r in inputs if r["acquisition_order"] == index + 1 and r["baseline_method"] == "linear_drift_corrected"]
        d = np.array([r["distance_nm"] for r in curve])
        y = np.array([r["force_median_pN"] for r in curve])
        show = (d >= fit["dmin_nm"]) & (d <= 250)
        color = chronology.COLORS[fit["nominal_speed_um_per_s"]]
        ax.scatter(d[show], y[show], color=color, s=9, alpha=0.9)
        grid = np.linspace(fit["dmin_nm"], 250, 250)
        ax.plot(grid, predict(grid, result_parameters(fit)), color="#24364a", lw=1.1)
        ax.axhline(0, color="0.6", lw=0.5)
        ax.set_title(f"#{index+1} | {fit['nominal_speed_um_per_s']:g} µm/s | {fit['dmin_nm']:g}–250 nm", fontsize=9, loc="left")
        ax.text(0.97, 0.91, f"λ={fit['lambda_D_apparent_nm']:.1f} nm\n|ψ|={fit['surface_potential_magnitude_apparent_mV']:.1f} mV\nR²={fit['r2']:.3f}", transform=ax.transAxes, ha="right", va="top", fontsize=8,
                color="#a52a2a" if not fit["fit_quality_pass"] else "0.2")
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.15)
        if index // 6 == 5:
            ax.set_xlabel("D (nm)", fontsize=9)
        if index % 6 == 0:
            ax.set_ylabel("Force (pN)", fontsize=9)
    fig.suptitle("All 36 map-median force curves and primary PB fits | red parameter labels = weak or unusable estimates", fontsize=13)
    fig.savefig(out / "all_map_PB_fits.png")
    plt.close(fig)
    fig, axes = plt.subplots(2, 1, figsize=(13.2, 7), sharex=True, layout="constrained")
    axes[0].plot(times, [r["dmin_nm"] for r in primary], "o-", color="#263b53", label="Safe fit lower limit")
    axes[0].plot(times, [r["detected_snap_max_nm"] for r in primary], ".--", color="#b35e3f", label="Largest detected snap-in distance")
    axes[0].set_ylabel("Distance (nm)")
    axes[0].legend(frameon=False)
    axes[1].plot(times, [r["r2"] for r in primary], "o-", color="#263b53")
    axes[1].axhline(0.7, ls=":", color="#b35e3f", label="R² quality threshold")
    axes[1].set_ylabel("Fit R²")
    axes[1].legend(frameon=False)
    for ax in axes:
        chronology.axes_chronology(ax, inventory, "map_midpoint_time")
        ax.grid(axis="y", alpha=0.2)
    axes[1].set_xlabel("Acquisition time (instrument UTC+02:00)")
    axes[0].set_title("Fit-window and residual diagnostics", loc="left", fontweight="bold")
    fig.savefig(out / "fit_window_quality_vs_time.png")
    plt.close(fig)


def write_report(primary: list[dict], results: list[dict], validation: dict) -> None:
    lines = ["# 08-09-26 水中 Debye length 与 surface potential 拟合", "",
        "本次计算复用 `analysis/fit_glycerol_surface_forces.py` 的非线性 PB 球—平板力模型。每个点对应一张 map 的中位力曲线拟合；不是 64 个逐曲线拟合参数的中位数。横轴是实际 map 起止时间的中点（仪器 UTC+02:00）。", "",
        "## 结果", "",
        f"36 张 map 均完成拟合；{sum(bool(r['fit_quality_pass']) for r in primary)}/36 通过残差和相对于像素异质性尺度的参数约束检查，即第 19–36 张。筛查通过不等同于验证平衡电双层模型。",
        "",
        "| Map 范围 | 表观 λ_D 中位数 [最小, 最大] / nm | 表观 ψ_s 幅值中位数 [最小, 最大] / mV |", "|---|---:|---:|"]
    for low, high in ((1, 12), (13, 18), (19, 36)):
        selected = [r for r in primary if low <= r["acquisition_order"] <= high and r["fit_quality_pass"]]
        if selected:
            l = [r["lambda_D_apparent_nm"] for r in selected]
            v = [r["surface_potential_magnitude_apparent_mV"] for r in selected]
            lines.append(f"| {low}–{high} | {np.median(l):.2f} [{min(l):.2f}, {max(l):.2f}] | {np.median(v):.2f} [{min(v):.2f}, {max(v):.2f}] |")
    lines.extend(["", "第 1–18 张的最优拟合参数保留在 CSV，但都未通过上述约束检查。主图空心点表示参数约束较弱的最优值；第 7、11、15、18 张因参数触边、残差或弱信号问题，仅在底部标红叉。红叉是状态标记，不是零值，也不是原始数据缺失。第 7、11、15 张达到约 250 mV 的优化上界，不能解释成真实电势达到了 250 mV。", "",
        "第 1–18 张起始距离为 40–80 nm，第 19–36 张为 20–25 nm，两段可用距离范围不同。统一为 80–250 nm 后，36 张均未通过像素异质性尺度下的参数约束检查。因此不能依据主图断言第 18→19 处 Debye length 或表面电势发生了可定量确认的跳变。上述表中范围是 map 间范围，不是置信区间。", "",
        "![参数随时间](figures/debye_length_surface_potential_vs_time.png)", "",
        "## 输入、模型与拟合窗口", "",
        "- 输入：已按 cantilever2 的 `k=0.2736 ± 0.0063 N/m` 及水中全局 InvOLS 重建的 36×64 条 approach 曲线。逐曲线取 5 nm 距离箱中位数，再在每张 map 内取像素中位数。保留缺测；拟合箱要求至少 48 个像素。与上游中位数 CSV 的逐箱数值及像素数已核对。",
        "- 几何与介质：`R=4.54685 µm` 为用户本次明确授权沿用，未重新测量；`T=25.6 °C` 沿用仓库模型，非本批元数据实测温度；水 `ε_r=78.5`，`A_H=2.4×10⁻²¹ J`。",
        "- 模型：对称 1:1 电解质、球与平面具有相同恒定表面电势；`F(D)=F_PB(D; λ_D, |ψ_s|)−A_H R/(6D²)+B`。B 为各 map 的常数残余基线。力模型只识别电势幅值；未赋予负号，也未将它视为电动学 ζ 电势测量。",
        "- 主拟合窗口上限 250 nm。下限为 `max(20, 5 ceil((max_detected_snap_nm+7.5)/5)) nm`：每个所用距离箱下边缘比该 map 最大已检测跳入距离至少高 5 nm。下限依 map 改变，避免混入机械跳跃；它并不证明其他自由支上的力完全满足 PB。",
        "- 优化沿用原脚本的 log(λ_D)、log(|ψ_s|)、B 参数化、两组初值、soft_l1 损失及边界：λ_D=1–1000 nm，|ψ_s|=0.1–250 mV，B=−500–500 pN。本适配脚本按实际优化的 soft_l1 cost 选择初值结果。未伪造 `speed_groups=3`，未调用原脚本要求三速度截距输入的 `fit_equilibrium()`。",
        "- 权重尺度采用像素力的 `1.4826×MAD`，底限为 `max(3 pN, 400–500 nm 中位力曲线的 robust MAD)`。由于已保存距离到 500 nm，底限参考窗不同于原脚本的 500–850 nm。没有将 64 个空间相关像素当作 64 次独立实验，也没有除以 √64。",
        "", "## 解释边界与敏感性", "",
        "这些是每张 map 在 1/2/4 µm/s 有限接近速度下、所选力模型和基线条件下的表观参数。未把随时间变化的各 map 合并为平衡力，也未做三速度零速外推。", "",
        "保存了五种敏感性结果：相同距离坐标的常数远场基线；下限增加 10 nm；上限降至 200 nm；所有 map 共用 80–250 nm；扣除 `6πηR²U/D` 的 no-slip drainage（U 取逐曲线实测 gap speed，η 沿用水的 Cheng correlation）；理论阻力先在各曲线原始远场采样点按 scanner height 拟合并扣除线性部分，再从逐像素力中扣除，最后取中位数。该远场投影近似匹配实验基线算子；没有假定 median(F−Fhyd)=median(F)−median(Fhyd)。最后一种只是固定物理模型的敏感性，不证明无滑移或阻力修正正确。", "",
        "`fit_quality_pass` 检查优化结束、参数边界、局部 Jacobian、R²≥0.70、信号≥5×底噪及局部 log 参数尺度≤0.5。局部尺度从按像素 MAD 加权的 Jacobian 得到，不因中位曲线残差小就将尺度缩小；协方差诊断因子为 `max(加权残差平方和/(箱数−3), 1)`。该阈值是一项保守的描述性筛查。CSV 的局部尺度不是实验置信区间；未将距离箱当作独立实验，也未包含 InvOLS、k、距离零点、球半径及表面边界条件的全部系统误差。", "",
        "no-slip 比较的 U 是已有 20–200 nm gap-speed 回归的代理量；早期该区间包含跳跃运动，因此它只用于近似敏感性，不用于验证黏性阻力机制。", "",
        "### 第 19–36 张的敏感性摘要", "",
        "下表中位数使用这 18 张的最优值；检查通过列单独给出通过数量。变体间差异不是置信区间。", "",
        "| 方案 | λ_D 中位数 / nm | 电势幅值中位数 / mV | 检查通过 |", "|---|---:|---:|---:|"])
    names = {"primary": "主拟合", "constant_baseline": "常数远场基线", "lower_plus_10nm": "起始距离增加 10 nm",
             "upper_200nm": "终止距离降至 200 nm", "no_slip_hydrodynamic": "近似 no-slip 修正", "common_80_250nm": "共用 80–250 nm"}
    for variant in ("primary", "constant_baseline", "lower_plus_10nm", "upper_200nm", "no_slip_hydrodynamic", "common_80_250nm"):
        selected = [r for r in results if r["variant"] == variant and r["acquisition_order"] >= 19]
        passed = sum(bool(r["fit_quality_pass"]) for r in selected)
        if passed:
            lam = np.median([r["lambda_D_apparent_nm"] for r in selected])
            psi = np.median([r["surface_potential_magnitude_apparent_mV"] for r in selected])
            lines.append(f"| {names[variant]} | {lam:.2f} | {psi:.2f} | {passed}/18 |")
        else:
            lines.append(f"| {names[variant]} | 不可可靠确定 | 不可可靠确定 | 0/18 |")
    lines.extend(["", "起始距离增加 10 nm 的变体中，第 25 张未通过异质性尺度检查；主拟合并未出现单像素图中那样的离群参数。", "",
        "![敏感性](figures/parameter_sensitivity_vs_time.png)", "", "![全部拟合](figures/all_map_PB_fits.png)", "",
        "## 文件", "", "- [DATA_QUALITY_COMPARISON.md](DATA_QUALITY_COMPARISON.md)：前 18 张与后 18 张的内部及组间质量比较；后组统一 25–250 nm 的补充诊断保存在 `quality_comparison/`，与本报告的原主拟合分开。", "- `fit_parameters_per_map.csv`：主结果，36 行。", "- `fit_parameters_all_variants.csv`：全部 216 个拟合及诊断。", "- `map_median_fit_inputs.csv`：逐距离箱输入、权重、像素数与安全窗。", "- `fitted_force_curves.csv`：各拟合曲线。", "- `parameter_sensitivity.csv`：相对主结果变化。", "- `provenance.json`、`artifact_manifest.sha256`：来源、代码、环境及校验。", "",
        "重现：在仓库根目录运行 `python -X utf8 analysis/fit_08_09_26_water.py --workers 4`。"])
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    reuse = parser.add_mutually_exclusive_group()
    reuse.add_argument("--refresh-diagnostics", action="store_true",
                       help="Recompute local geometry and plots at saved optima, with input hashes checked")
    reuse.add_argument("--export-only", action="store_true",
                       help="Export saved fits and diagnostics without recomputing them")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    configure_model()
    if args.refresh_diagnostics or args.export_only:
        saved = json.loads((OUT / "provenance.json").read_text(encoding="utf-8"))
        for name, digest in saved["input_sha256"].items():
            if base.sha256_file(ROOT / name) != digest:
                raise AssertionError(f"Input changed: {name}")
        if (saved["radius_m"] != RADIUS_M or saved["temperature_C"] != base.TEMPERATURE_C
                or saved["epsilon_r"] != base.EPSILON_R[0] or saved["Hamaker_J"] != base.HAMAKER_J
                or saved["code_sha256"][Path(base.__file__).relative_to(ROOT).as_posix()] != base.sha256_file(Path(base.__file__))):
            raise AssertionError("Model constants changed; full refit required")
        inventory = read_csv(DATA / "map_inventory_QC.csv")
        inputs = read_csv(OUT / "map_median_fit_inputs.csv")
        results = read_csv(OUT / "fit_parameters_all_variants.csv")
        if args.refresh_diagnostics:
            refresh_diagnostics(results, inputs)
        finish(inventory, inputs, results, saved["validation"])
        print("Exports refreshed; optimized parameters unchanged", flush=True)
        return
    inventory, inputs, tasks, validation = prepare()
    base.write_csv(OUT / "map_median_fit_inputs.csv", inputs)
    print(f"Validated {len(inventory)} map medians; running {len(tasks)} fits", flush=True)
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(fit_one, task) for task in tasks]
        for future in as_completed(futures):
            results.append(future.result())
            if len(results) % 18 == 0:
                print(f"Completed {len(results)}/{len(tasks)} fits", flush=True)
    results.sort(key=lambda row: (row["acquisition_order"], VARIANTS.index(row["variant"])))
    finish(inventory, inputs, results, validation)
    print(json.dumps({"output": str(OUT), "primary_quality_pass": validation["primary_quality_pass_maps"],
                      "flagged": validation["primary_flagged_maps"]}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()

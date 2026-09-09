#!/usr/bin/env python3
"""Fit every map2--5 force curve with an equal-sphere HHF/DLVO model.

This is a deliberately model-mismatched diagnostic: the experiment is a
sphere--plane geometry, while the user requested the reference repository's
equal-potential colloid--colloid formula.  Every one of the 64 signed F(D)
curves in each same-location map is fitted independently over 20--200 nm.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.constants import Boltzmann, elementary_charge, epsilon_0
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import analyze_31_08_26_maps as raw  # noqa: E402
import fit_glycerol_surface_forces as base  # noqa: E402


DATA_ROOT = ROOT / "31-08-26"
RESULTS = ROOT / "analysis" / "maps_31_08_26_colloid_colloid_results"
FIGURES = RESULTS / "figures"
CALIBRATION_SUMMARY = (
    ROOT / "analysis" / "palindrome_27_08_26_full_results" / "calibration_summary.csv"
)

FIT_MIN_NM = 20.0
FIT_MAX_NM = 200.0
EPSILON_R_WATER = 78.5
RADIUS_1_M = raw.CONTEXT_PROBE_RADIUS_M
RADIUS_2_M = raw.CONTEXT_PROBE_RADIUS_M
EFFECTIVE_RADIUS_M = 1.0 / (1.0 / RADIUS_1_M + 1.0 / RADIUS_2_M)
HAMAKER_J = base.HAMAKER_J
TEMPERATURE_K = base.TEMPERATURE_K
WATER_VISCOSITY_PA_S = (
    base.cheng_viscosity_mPa_s(0.0, raw.CONTEXT_TEMPERATURE_C) * 1e-3
)

LAMBDA_BOUNDS_NM = (0.5, 300.0)
ZETA_BOUNDS_MV = (0.1, 250.0)
BASELINE_BOUNDS_PN = (-500.0, 500.0)
MIN_NOISE_PN = 2.0
PRIMARY_VARIANT = "far_linear_fixed_two_sphere_hydrodynamics"
FIT_VARIANTS = {
    PRIMARY_VARIANT: ("|line", True),
    "far_linear_no_hydrodynamics": ("|line", False),
    "far_constant_fixed_two_sphere_hydrodynamics": ("|constant", True),
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def equal_sphere_force_components_pN(
    distance_nm: np.ndarray,
    lambda_D_nm: float,
    zeta_mV: float,
    gap_speed_um_per_s: float,
    include_hydrodynamics: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return equal-potential HHF EDL, vdW, and fixed hydrodynamic forces."""
    distance_m = np.asarray(distance_nm, dtype=np.float64) * 1e-9
    if not np.all(np.isfinite(distance_m)) or np.any(distance_m <= 0.0):
        raise ValueError("separation must be finite and positive")
    lambda_m = float(lambda_D_nm) * 1e-9
    zeta_V = abs(float(zeta_mV)) * 1e-3
    if not np.isfinite(lambda_m) or lambda_m <= 0.0:
        raise ValueError("Debye length must be finite and positive")
    if not np.isfinite(zeta_V):
        raise ValueError("surface potential must be finite")
    kappa = 1.0 / lambda_m
    x = kappa * distance_m
    reciprocal = np.where(x > 40.0, np.exp(-x), 1.0 / (np.exp(x) + 1.0))
    edl_N = (
        4.0
        * np.pi
        * EFFECTIVE_RADIUS_M
        * EPSILON_R_WATER
        * epsilon_0
        * kappa
        * zeta_V**2
        * reciprocal
    )
    vdw_N = -HAMAKER_J * EFFECTIVE_RADIUS_M / (6.0 * distance_m**2)
    if include_hydrodynamics:
        speed_m_per_s = float(gap_speed_um_per_s) * 1e-6
        if not np.isfinite(speed_m_per_s) or speed_m_per_s <= 0.0:
            raise ValueError("gap-closing speed must be finite and positive")
        hydrodynamic_N = (
            6.0
            * np.pi
            * WATER_VISCOSITY_PA_S
            * EFFECTIVE_RADIUS_M**2
            * speed_m_per_s
            / distance_m
        )
    else:
        hydrodynamic_N = np.zeros_like(distance_m)
    components = tuple(value * 1e12 for value in (edl_N, vdw_N, hydrodynamic_N))
    if not all(np.all(np.isfinite(value)) for value in components):
        raise FloatingPointError("colloid-colloid model produced non-finite force")
    return components  # type: ignore[return-value]


def fit_one_curve(
    distance_nm: np.ndarray,
    observed_force_pN: np.ndarray,
    gap_speed_um_per_s: float,
    far_noise_pN: float,
    include_hydrodynamics: bool,
) -> tuple[dict, np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    distance_nm = np.asarray(distance_nm, dtype=np.float64)
    observed_force_pN = np.asarray(observed_force_pN, dtype=np.float64)
    if distance_nm.shape != observed_force_pN.shape or distance_nm.size < 30:
        raise ValueError("a fit requires at least 30 aligned distance-force bins")
    if not np.all(np.isfinite(distance_nm)) or not np.all(
        np.isfinite(observed_force_pN)
    ):
        raise ValueError("fit arrays must be finite")
    noise_pN = max(abs(float(far_noise_pN)), MIN_NOISE_PN)
    baseline_start = float(np.median(observed_force_pN[-8:]))
    baseline_start = float(np.clip(baseline_start, *BASELINE_BOUNDS_PN))

    lower = np.asarray(
        [
            math.log(LAMBDA_BOUNDS_NM[0]),
            math.log(ZETA_BOUNDS_MV[0]),
            BASELINE_BOUNDS_PN[0],
        ],
        dtype=np.float64,
    )
    upper = np.asarray(
        [
            math.log(LAMBDA_BOUNDS_NM[1]),
            math.log(ZETA_BOUNDS_MV[1]),
            BASELINE_BOUNDS_PN[1],
        ],
        dtype=np.float64,
    )

    def residual(parameters: np.ndarray) -> np.ndarray:
        lambda_D_nm = math.exp(float(parameters[0]))
        zeta_mV = math.exp(float(parameters[1]))
        edl, vdw, hydrodynamic = equal_sphere_force_components_pN(
            distance_nm,
            lambda_D_nm,
            zeta_mV,
            gap_speed_um_per_s,
            include_hydrodynamics,
        )
        prediction = edl + vdw + hydrodynamic + float(parameters[2])
        return (prediction - observed_force_pN) / noise_pN

    candidates = []
    for lambda_start_nm in (3.0, 10.0, 30.0, 100.0):
        for zeta_start_mV in (10.0, 40.0, 120.0):
            result = least_squares(
                residual,
                x0=np.asarray(
                    [
                        math.log(lambda_start_nm),
                        math.log(zeta_start_mV),
                        baseline_start,
                    ],
                    dtype=np.float64,
                ),
                bounds=(lower, upper),
                loss="soft_l1",
                f_scale=1.0,
                xtol=1e-11,
                ftol=1e-11,
                gtol=1e-11,
                max_nfev=3000,
            )
            candidates.append(result)
    result = min(candidates, key=lambda item: float(item.cost))
    lambda_D_nm = math.exp(float(result.x[0]))
    zeta_mV = math.exp(float(result.x[1]))
    baseline_pN = float(result.x[2])
    components = equal_sphere_force_components_pN(
        distance_nm,
        lambda_D_nm,
        zeta_mV,
        gap_speed_um_per_s,
        include_hydrodynamics,
    )
    prediction = sum(components) + baseline_pN
    raw_residual = observed_force_pN - prediction
    rss = float(np.sum(raw_residual**2))
    total = float(np.sum((observed_force_pN - np.mean(observed_force_pN)) ** 2))
    r_squared = float("nan") if total <= 0.0 else 1.0 - rss / total
    lambda_boundary = bool(
        np.isclose(lambda_D_nm, LAMBDA_BOUNDS_NM[0], rtol=2e-3)
        or np.isclose(lambda_D_nm, LAMBDA_BOUNDS_NM[1], rtol=2e-3)
    )
    zeta_boundary = bool(
        np.isclose(zeta_mV, ZETA_BOUNDS_MV[0], rtol=2e-3)
        or np.isclose(zeta_mV, ZETA_BOUNDS_MV[1], rtol=2e-3)
    )
    baseline_boundary = bool(
        np.isclose(baseline_pN, BASELINE_BOUNDS_PN[0], rtol=0.0, atol=0.5)
        or np.isclose(baseline_pN, BASELINE_BOUNDS_PN[1], rtol=0.0, atol=0.5)
    )
    dimensionless_potential = (
        elementary_charge * zeta_mV * 1e-3 / (Boltzmann * TEMPERATURE_K)
    )
    fit_status = (
        "solver_failure"
        if not result.success
        else "parameter_boundary"
        if lambda_boundary or zeta_boundary or baseline_boundary
        else "interior_solution"
    )
    summary = {
        "fit_success": bool(result.success),
        "solver_status": int(result.status),
        "solver_message": str(result.message),
        "fit_status": fit_status,
        "n_distance_bins": int(distance_nm.size),
        "lambda_D_nm": lambda_D_nm,
        "zeta_magnitude_mV": zeta_mV,
        "zeta_signed_silica_mV": -zeta_mV,
        "baseline_pN": baseline_pN,
        "rmse_pN": float(math.sqrt(np.mean(raw_residual**2))),
        "r2": r_squared,
        "far_noise_scale_pN": noise_pN,
        "ordinary_reduced_chi2_diagnostic": float(
            np.sum((raw_residual / noise_pN) ** 2) / (distance_nm.size - 3)
        ),
        "soft_l1_cost": float(result.cost),
        "optimality": float(result.optimality),
        "nfev": int(result.nfev),
        "lambda_boundary": lambda_boundary,
        "zeta_boundary": zeta_boundary,
        "baseline_boundary": baseline_boundary,
        "dimensionless_potential_epsi_over_kT": dimensionless_potential,
        "reference_25mV_ceiling_satisfied": bool(zeta_mV <= 25.0),
    }
    return summary, prediction, components


def summarize_fits(fit_rows: list[dict]) -> list[dict]:
    summaries: list[dict] = []
    for variant in FIT_VARIANTS:
        for map_order in range(2, 6):
            selected = [
                row
                for row in fit_rows
                if row["fit_variant"] == variant
                and int(row["map_order"]) == map_order
                and bool(row["fit_success"])
            ]
            if not selected:
                raise RuntimeError(f"no successful fits for map{map_order} {variant}")
            interior = [row for row in selected if row["fit_status"] == "interior_solution"]
            record: dict[str, object] = {
                "fit_variant": variant,
                "map_order": map_order,
                "relative_time_min": float(selected[0]["relative_time_min"]),
                "curves": 64,
                "successful_fits": len(selected),
                "interior_fits": len(interior),
                "boundary_fraction": float(
                    np.mean([row["fit_status"] == "parameter_boundary" for row in selected])
                ),
                "reference_25mV_ceiling_satisfied_fraction": float(
                    np.mean([bool(row["reference_25mV_ceiling_satisfied"]) for row in selected])
                ),
            }
            for field in ("lambda_D_nm", "zeta_magnitude_mV", "baseline_pN", "rmse_pN", "r2"):
                values = np.asarray([float(row[field]) for row in selected], dtype=np.float64)
                record[f"{field}_median"] = float(np.median(values))
                record[f"{field}_q25"] = float(np.quantile(values, 0.25))
                record[f"{field}_q75"] = float(np.quantile(values, 0.75))
            summaries.append(record)
    return summaries


def plot_per_curve_fits(
    fit_rows: list[dict], prediction_rows: list[dict]
) -> None:
    primary_fits = [row for row in fit_rows if row["fit_variant"] == PRIMARY_VARIANT]
    primary_predictions = [
        row for row in prediction_rows if row["fit_variant"] == PRIMARY_VARIANT
    ]
    figure, axes = plt.subplots(2, 2, figsize=(14.0, 10.0), sharex=True)
    for axis, map_order in zip(axes.ravel(), range(2, 6), strict=True):
        map_predictions = [
            row for row in primary_predictions if int(row["map_order"]) == map_order
        ]
        observed_matrix = []
        fitted_matrix = []
        for point_index in range(64):
            rows = sorted(
                (
                    row
                    for row in map_predictions
                    if int(row["point_index"]) == point_index
                ),
                key=lambda row: float(row["distance_nm"]),
            )
            if len(rows) < 30:
                continue
            distance = np.asarray([row["distance_nm"] for row in rows], dtype=float)
            observed = np.asarray([row["observed_force_pN"] for row in rows], dtype=float)
            fitted = np.asarray([row["fitted_total_force_pN"] for row in rows], dtype=float)
            observed_matrix.append(observed)
            fitted_matrix.append(fitted)
            axis.plot(distance, observed, color="#457b9d", alpha=0.075, lw=0.75)
            axis.plot(distance, fitted, color="#e76f51", alpha=0.075, lw=0.75)
        observed_array = np.asarray(observed_matrix)
        fitted_array = np.asarray(fitted_matrix)
        axis.plot(distance, np.median(observed_array, axis=0), color="#1d3557", lw=2.5, label="median observed")
        axis.plot(distance, np.median(fitted_array, axis=0), color="#d1495b", lw=2.5, label="median fitted")
        fit_map = [row for row in primary_fits if int(row["map_order"]) == map_order]
        interior_fraction = np.mean([row["fit_status"] == "interior_solution" for row in fit_map])
        axis.set_title(f"Map {map_order}: 64 independent fits; interior {interior_fraction:.1%}")
        axis.axhline(0.0, color="0.35", lw=0.8)
        axis.grid(alpha=0.22)
        axis.set_xlabel("Separation D (nm)")
        axis.set_ylabel("Force (pN)")
    axes[0, 0].legend(frameon=False)
    figure.suptitle("Equal-sphere HHF/DLVO fits to all original map curves (20–200 nm)")
    figure.tight_layout()
    figure.savefig(FIGURES / "all_256_colloid_colloid_curve_fits.png", dpi=230)
    plt.close(figure)


def plot_parameter_distributions(summary_rows: list[dict]) -> None:
    primary = [row for row in summary_rows if row["fit_variant"] == PRIMARY_VARIANT]
    times = np.asarray([row["relative_time_min"] for row in primary], dtype=float)
    figure, axes = plt.subplots(2, 1, figsize=(9.5, 8.8), sharex=True)
    for axis, field, label, color in (
        (axes[0], "lambda_D_nm", "Apparent Debye length (nm)", "#2a9d8f"),
        (axes[1], "zeta_magnitude_mV", r"Apparent $|\psi|$ (mV)", "#d1495b"),
    ):
        median = np.asarray([row[f"{field}_median"] for row in primary], dtype=float)
        q25 = np.asarray([row[f"{field}_q25"] for row in primary], dtype=float)
        q75 = np.asarray([row[f"{field}_q75"] for row in primary], dtype=float)
        axis.errorbar(
            times,
            median,
            yerr=np.vstack([median - q25, q75 - median]),
            fmt="o-",
            color=color,
            capsize=4,
            lw=2.2,
        )
        axis.set_ylabel(label)
        axis.grid(alpha=0.23)
    axes[1].axhline(25.0, color="0.35", ls="--", lw=1.2, label="reference 25 mV ceiling")
    axes[1].legend(frameon=False)
    axes[1].set_xlabel("Time from start of Map 2 (min)")
    axes[1].set_xticks(times)
    axes[1].set_xticklabels([f"{value:.2f}" for value in times])
    figure.suptitle("Per-curve apparent equal-sphere fit parameters: median and IQR")
    figure.tight_layout()
    figure.savefig(FIGURES / "apparent_colloid_colloid_parameters_vs_time.png", dpi=230)
    plt.close(figure)


def render_report(
    spring_N_per_m: float,
    global_sensitivity_nm_per_V: float,
    summary_rows: list[dict],
) -> str:
    primary = [row for row in summary_rows if row["fit_variant"] == PRIMARY_VARIANT]
    lines = [
        "# 31-08-26: 256 independent colloid-colloid fits",
        "",
        "## Scope",
        "",
        "Map2–map5 each contribute all 64 original signed F(D) curves. Each curve is reconstructed from raw vDeflection/measuredHeight, binned every 5 nm, and fitted independently over 20–200 nm. No map-level median or four-point |F| time curve is used as a fit input.",
        "",
        f"Force reconstruction uses global liquid InvOLS `{global_sensitivity_nm_per_V:.6f} nm/V` and independently calibrated `k={spring_N_per_m:.9f} N/m`.",
        "",
        "## Requested colloid-colloid model",
        "",
        "This is a deliberate geometry-mismatch diagnostic because the experiment is silica sphere–plane, while this fit assumes two identical silica spheres with `R1=R2=4.546849 µm` and `R_eff=R/2=2.273424 µm`.",
        "",
        "The equal-potential linearized HHF and two-sphere Derjaguin terms are",
        "",
        r"`F_EDL = 4π R_eff ε κ ψ² / [exp(κD)+1]`,",
        "",
        r"`F_vdW = -A_H R_eff/(6D²)`,",
        "",
        r"`F_hyd = 6πη R_eff² v_gap/D`,",
        "",
        "with a fitted additive baseline. Primary fits include the fixed pure-water no-slip hydrodynamic term; no-hydrodynamics and far-constant branches are saved as sensitivity variants. Free parameters are apparent `lambda_D`, equal-surface `|psi|`, and baseline.",
        "",
        "## Primary per-map distributions",
        "",
        "Values are medians [IQR] over 64 independent curve fits.",
        "",
        "| map | time (min) | apparent lambda_D (nm) | apparent |psi| (mV) | baseline (pN) | RMSE (pN) | R2 | interior fits | <=25 mV |",
        "|---:|---:|:---|:---|:---|:---|:---|---:|---:|",
    ]
    for row in primary:
        lines.append(
            f"| {row['map_order']} | {row['relative_time_min']:.2f} | "
            f"{row['lambda_D_nm_median']:.2f} [{row['lambda_D_nm_q25']:.2f}, {row['lambda_D_nm_q75']:.2f}] | "
            f"{row['zeta_magnitude_mV_median']:.1f} [{row['zeta_magnitude_mV_q25']:.1f}, {row['zeta_magnitude_mV_q75']:.1f}] | "
            f"{row['baseline_pN_median']:+.1f} [{row['baseline_pN_q25']:+.1f}, {row['baseline_pN_q75']:+.1f}] | "
            f"{row['rmse_pN_median']:.1f} [{row['rmse_pN_q25']:.1f}, {row['rmse_pN_q75']:.1f}] | "
            f"{row['r2_median']:.3f} [{row['r2_q25']:.3f}, {row['r2_q75']:.3f}] | "
            f"{int(row['interior_fits'])}/64 | {row['reference_25mV_ceiling_satisfied_fraction']:.1%} |"
        )
    lines += [
        "",
        "## Interpretation boundary",
        "",
        "- These are per-curve apparent parameters under the user-requested wrong geometry, not equilibrium sphere–plane zeta potentials or bulk Debye lengths.",
        "- Equal-sphere conservative-force amplitudes are half the actual sphere–plane Derjaguin amplitudes at the same physical radius. The fitted potential therefore compensates for the geometry choice.",
        "- The HHF expression is linearized PB. Fits above the 25 mV reference ceiling are retained numerically but flagged as outside the reference implementation's operational potential range.",
        "- Adjacent 5 nm bins from one curve are correlated, and hydrodynamics/history/baseline remain model-conditioned. Across-curve IQR is more relevant than local optimizer precision.",
        "",
        "## Outputs",
        "",
        "- `per_curve_colloid_colloid_fits.csv`: all fit parameters and diagnostics.",
        "- `per_curve_colloid_colloid_predictions.csv`: observed and fitted values at every included distance bin.",
        "- `map_colloid_colloid_fit_summary.csv`: map-level medians and IQRs for all variants.",
        "- `figures/all_256_colloid_colloid_curve_fits.png` and `figures/apparent_colloid_colloid_parameters_vs_time.png`.",
        "- `provenance.json` and `artifact_manifest.sha256`.",
        "",
    ]
    return "\n".join(lines)


def create_manifest(paths: list[Path], destination: Path) -> None:
    unique = sorted(set(path.resolve() for path in paths), key=lambda path: str(path))
    destination.write_text(
        "\n".join(
            f"{sha256_file(path)}  {path.relative_to(ROOT)}" for path in unique
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    paths = sorted(DATA_ROOT.glob("*.jpk-force-map"))
    if len(paths) != 5:
        raise RuntimeError(f"expected five raw maps, got {len(paths)}")
    sources = [base.load_source(path.resolve(), 0) for path in paths]
    sources.sort(key=lambda source: source.timestamp)
    inferred, spring_N_per_m, spring_sd_N_per_m, fingerprints = raw.infer_cantilever(sources)
    contact_rows, global_sensitivity_m_per_V = raw.contact_analysis(sources)
    pixel_rows, map_rows, _, _, matrices = raw.analyze_maps(sources, spring_N_per_m)
    if inferred != "D5":
        raise RuntimeError(f"unexpected inferred cantilever {inferred}")
    reference_coordinates = raw.physical_coordinates_um(sources[1].path, 64)
    for source in sources[2:]:
        if not np.array_equal(
            reference_coordinates, raw.physical_coordinates_um(source.path, 64)
        ):
            raise RuntimeError("map2--5 coordinates are not exactly identical")

    distance_mask = (
        (raw.BIN_CENTERS_NM >= FIT_MIN_NM)
        & (raw.BIN_CENTERS_NM <= FIT_MAX_NM)
    )
    distance_nm = raw.BIN_CENTERS_NM[distance_mask]
    if distance_nm.size != 37:
        raise RuntimeError(f"expected 37 fit bins, got {distance_nm.size}")
    pixel_lookup = {
        (int(row["map_order"]), int(row["point_index"])): row for row in pixel_rows
    }
    relative_times = [
        (source.timestamp - sources[1].timestamp).total_seconds() / 60.0
        for source in sources[1:]
    ]

    fit_rows: list[dict] = []
    prediction_rows: list[dict] = []
    for map_order, source, relative_time_min in zip(
        range(2, 6), sources[1:], relative_times, strict=True
    ):
        source_key = raw.relative_source(source.path)
        for point_index in range(64):
            pixel = pixel_lookup[(map_order, point_index)]
            gap_speed = float(pixel["gap_speed_20_200nm_um_per_s"])
            far_noise = float(pixel["far_noise_pN"])
            row_index = int(pixel["row"])
            column_index = int(pixel["column"])
            for variant, (matrix_suffix, include_hydrodynamics) in FIT_VARIANTS.items():
                observed = np.asarray(
                    matrices[source_key + matrix_suffix][point_index, distance_mask],
                    dtype=np.float64,
                )
                valid = np.isfinite(observed)
                if np.count_nonzero(valid) < 30:
                    raise RuntimeError(
                        f"map{map_order} point{point_index}: only {np.count_nonzero(valid)} bins"
                    )
                selected_distance = distance_nm[valid]
                selected_observed = observed[valid]
                fit, prediction, components = fit_one_curve(
                    selected_distance,
                    selected_observed,
                    gap_speed,
                    far_noise,
                    include_hydrodynamics,
                )
                common = {
                    "source": source_key,
                    "map_order": map_order,
                    "relative_time_min": float(relative_time_min),
                    "point_index": point_index,
                    "row": row_index,
                    "column": column_index,
                    "fit_variant": variant,
                    "baseline_method": "linear_drift_corrected"
                    if matrix_suffix == "|line"
                    else "far_constant_referenced",
                    "include_fixed_two_sphere_hydrodynamics": include_hydrodynamics,
                    "gap_speed_um_per_s": gap_speed,
                    "fit_min_nm": FIT_MIN_NM,
                    "fit_max_nm": FIT_MAX_NM,
                }
                fit_rows.append({**common, **fit})
                edl, vdw, hydrodynamic = components
                for values in zip(
                    selected_distance,
                    selected_observed,
                    edl,
                    vdw,
                    hydrodynamic,
                    prediction,
                    strict=True,
                ):
                    d_value, observed_value, edl_value, vdw_value, hyd_value, predicted_value = values
                    prediction_rows.append(
                        {
                            **common,
                            "distance_nm": float(d_value),
                            "observed_force_pN": float(observed_value),
                            "fitted_edl_force_pN": float(edl_value),
                            "fixed_vdw_force_pN": float(vdw_value),
                            "fixed_hydrodynamic_force_pN": float(hyd_value),
                            "fitted_baseline_pN": float(fit["baseline_pN"]),
                            "fitted_total_force_pN": float(predicted_value),
                            "residual_observed_minus_fitted_pN": float(
                                observed_value - predicted_value
                            ),
                        }
                    )

    if len(fit_rows) != 4 * 64 * len(FIT_VARIANTS):
        raise RuntimeError(f"unexpected fit count {len(fit_rows)}")
    summary_rows = summarize_fits(fit_rows)
    write_csv(RESULTS / "per_curve_colloid_colloid_fits.csv", fit_rows)
    write_csv(
        RESULTS / "per_curve_colloid_colloid_predictions.csv", prediction_rows
    )
    write_csv(RESULTS / "map_colloid_colloid_fit_summary.csv", summary_rows)
    plot_per_curve_fits(fit_rows, prediction_rows)
    plot_parameter_distributions(summary_rows)
    (RESULTS / "REPORT.md").write_text(
        render_report(
            spring_N_per_m,
            global_sensitivity_m_per_V * 1e9,
            summary_rows,
        ),
        encoding="utf-8",
    )
    provenance = {
        "analysis": "all 64 curves per map fitted independently with requested equal-sphere colloid-colloid HHF/DLVO model",
        "claim_status": "apparent_model_mismatched_geometry_diagnostic",
        "input_hashes": {
            raw.relative_source(source.path): source.sha256 for source in sources
        },
        "raw_reconstruction_script": raw.relative_source(Path(raw.__file__).resolve()),
        "cantilever": {
            "inferred": inferred,
            "spring_constant_N_per_m": spring_N_per_m,
            "repeatability_sd_N_per_m": spring_sd_N_per_m,
            "fingerprints": fingerprints,
        },
        "global_liquid_InvOLS_nm_per_V": global_sensitivity_m_per_V * 1e9,
        "contact_fit_rows": len(contact_rows),
        "model": {
            "name": "equal_potential_equal_sphere_linear_HHF_plus_Derjaguin_vdW_and_optional_fixed_lubrication",
            "experimental_geometry": "sphere-plane",
            "fitted_geometry": "equal sphere-sphere at explicit user request",
            "R1_m": RADIUS_1_M,
            "R2_m": RADIUS_2_M,
            "R_eff_m": EFFECTIVE_RADIUS_M,
            "temperature_K": TEMPERATURE_K,
            "epsilon_r": EPSILON_R_WATER,
            "Hamaker_J": HAMAKER_J,
            "water_viscosity_Pa_s": WATER_VISCOSITY_PA_S,
            "primary_variant": PRIMARY_VARIANT,
            "fit_variants": FIT_VARIANTS,
            "fit_window_nm": [FIT_MIN_NM, FIT_MAX_NM],
            "distance_step_nm": 5.0,
            "free_parameters": ["lambda_D_nm", "zeta_magnitude_mV", "baseline_pN"],
            "loss": "soft_l1 with per-curve far-noise scale, minimum 2 pN",
            "parameter_bounds": {
                "lambda_D_nm": LAMBDA_BOUNDS_NM,
                "zeta_magnitude_mV": ZETA_BOUNDS_MV,
                "baseline_pN": BASELINE_BOUNDS_PN,
            },
        },
        "fit_counts": {
            "maps": 4,
            "curves_per_map": 64,
            "variants": len(FIT_VARIANTS),
            "total_fits": len(fit_rows),
            "primary_fits": sum(row["fit_variant"] == PRIMARY_VARIANT for row in fit_rows),
            "prediction_rows": len(prediction_rows),
        },
        "randomness": "none",
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "dependency_hashes": {
            raw.relative_source(Path(raw.__file__).resolve()): sha256_file(
                Path(raw.__file__).resolve()
            ),
            raw.relative_source(Path(base.__file__).resolve()): sha256_file(
                Path(base.__file__).resolve()
            ),
            raw.relative_source(CALIBRATION_SUMMARY): sha256_file(
                CALIBRATION_SUMMARY
            ),
        },
    }
    (RESULTS / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    artifacts = [
        path
        for path in RESULTS.rglob("*")
        if path.is_file() and path.name != "artifact_manifest.sha256"
    ]
    artifacts += [Path(__file__).resolve(), Path(raw.__file__).resolve(), CALIBRATION_SUMMARY]
    create_manifest(artifacts, RESULTS / "artifact_manifest.sha256")
    print(f"Wrote {RESULTS}")
    print(
        f"primary fits=256, all fits={len(fit_rows)}, "
        f"prediction rows={len(prediction_rows)}, R_eff={EFFECTIVE_RADIUS_M*1e6:.6f} um"
    )


if __name__ == "__main__":
    main()

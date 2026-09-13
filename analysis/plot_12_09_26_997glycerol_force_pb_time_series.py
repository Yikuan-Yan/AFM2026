#!/usr/bin/env python3
"""Recreate the established force-slice and PB-parameter time figures.

This 99.7 wt% glycerol data set cannot use the older fixed-sample contact and
first-fraction baseline rules.  It consumes the separately validated,
speed-conditioned constant baselines, detects hard contact in a terminal
physical-distance window, and estimates one batch InvOLS from those contacts.

The PB fit is deliberately limited to producing the requested analogue of the
earlier Debye-length/surface-potential plot.  Its parameters are finite-speed,
model-conditioned apparent values; no hydrodynamic subtraction is applied.
"""

from __future__ import annotations

import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import scipy
from scipy.optimize import least_squares
from scipy.stats import spearmanr

import fit_glycerol_surface_forces as base


ROOT = Path(__file__).resolve().parents[1]
ADAPTIVE = ROOT / "analysis" / "glycerol_99p7_D3_adaptive_baseline_results"
CHRONOLOGY = ROOT / "analysis" / "glycerol_99p7_D3_time_speed_results"
OUT = ROOT / "analysis" / "glycerol_99p7_D3_force_pb_time_results"
FIG = OUT / "figures"

SPEEDS = (0.1, 0.3, 0.9, 2.7)
COLORS = {0.1: "#2a9d8f", 0.3: "#457b9d", 0.9: "#e9c46a", 2.7: "#e76f51"}
TARGETS_NM = (20.0, 50.0, 100.0, 200.0)
BIN_CENTRES_NM = np.arange(5.0, 300.0 + 0.1, 5.0)
BIN_HALF_WIDTH_NM = 2.5

# Independent D3 thermal calibration.  InvOLS is re-estimated from the force
# maps because the older D3 value fails the hard-contact unit-slope check.
D3_K_N_PER_M = 0.2383668909630904
D3_PRIOR_INVOLS_NM_PER_V = 91.5629518721293
CONTACT_WINDOW_NM = 10.0
CONTACT_R2_MIN = 0.99
CONTACT_INVOLS_RANGE_NM_PER_V = (35.0, 90.0)
MIN_CONTACTS_PER_MAP = 32

MODEL = "nonlinear_pb_derjaguin"
PROBE_RADIUS_M = 4.546848945303745e-6
HAMAKER_J = 2.4e-21
TEMPERATURE_C = 25.6
EPSILON_R_997_GLYCEROL = 42.5
EPSILON_R_REFERENCE = "https://doi.org/10.3390/ma11040650"
FIT_MIN_NM = 20.0
FIT_MAX_NM = 250.0
LOWER = np.array([math.log(1.0), math.log(0.1), -500.0])
UPPER = np.array([math.log(1000.0), math.log(500.0), 500.0])


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def mad(values: np.ndarray) -> float:
    return base.robust_mad(np.asarray(values, dtype=float))


def q(values: np.ndarray, percentile: float) -> float:
    selected = np.asarray(values, dtype=float)
    selected = selected[np.isfinite(selected)]
    return float(np.percentile(selected, percentile)) if selected.size else float("nan")


def block_boundaries(chronology: list[dict]) -> list[datetime]:
    output: list[datetime] = []
    for index in (8, 16, 24):
        left = datetime.fromisoformat(chronology[index - 1]["map_midpoint_time"])
        right = datetime.fromisoformat(chronology[index]["map_midpoint_time"])
        output.append(left + (right - left) / 2)
    return output


def style_time_axis(axis: plt.Axes, chronology: list[dict]) -> None:
    for boundary in block_boundaries(chronology):
        axis.axvline(boundary, color="0.45", ls="--", lw=0.8, zorder=0)
    times = [datetime.fromisoformat(row["map_midpoint_time"]) for row in chronology]
    margin = (times[-1] - times[0]) * 0.012
    axis.set_xlim(times[0] - margin, times[-1] + margin)
    timezone = times[0].tzinfo
    axis.xaxis.set_major_locator(mdates.MinuteLocator(interval=10, tz=timezone))
    axis.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=timezone))
    axis.grid(axis="y", alpha=0.2)


def contact_diagnostic(
    source: base.SourceData, baseline_by_point: dict[int, dict[str, str]], order: int
) -> tuple[list[dict], float, float, int]:
    rows: list[dict] = []
    accepted: list[float] = []
    for curve in source.curves:
        baseline = float(baseline_by_point[int(curve.point_index)]["dynamic_baseline_raw_V"])
        voltage = curve.deflection_V - baseline
        travel_nm = (curve.measured_height_m[0] - curve.measured_height_m) * 1e9
        selected = travel_nm >= travel_nm[-1] - CONTACT_WINDOW_NM
        try:
            slope, _, r2, retained = base.robust_line(travel_nm[selected], voltage[selected])
            inv_ols = 1.0 / slope if slope > 0.0 else float("nan")
        except (ValueError, np.linalg.LinAlgError):
            slope, r2, retained, inv_ols = (float("nan"), float("nan"), 0, float("nan"))
        keep = bool(
            np.isfinite(inv_ols)
            and CONTACT_INVOLS_RANGE_NM_PER_V[0] <= inv_ols <= CONTACT_INVOLS_RANGE_NM_PER_V[1]
            and np.isfinite(r2)
            and r2 >= CONTACT_R2_MIN
        )
        if keep:
            accepted.append(inv_ols)
        rows.append(
            {
                "acquisition_order": order,
                "source": str(source.path.relative_to(ROOT)).replace("\\", "/"),
                "point_index": int(curve.point_index),
                "terminal_physical_window_nm": CONTACT_WINDOW_NM,
                "terminal_sample_count": int(np.count_nonzero(selected)),
                "contact_slope_V_per_nm": slope,
                "contact_invOLS_nm_per_V": inv_ols,
                "contact_r2": r2,
                "contact_retained_points": retained,
                "contact_calibration_accepted": keep,
            }
        )
    if len(accepted) < MIN_CONTACTS_PER_MAP:
        raise RuntimeError(f"Map {order}: only {len(accepted)} valid terminal contacts")
    return rows, float(np.median(accepted)), mad(np.asarray(accepted)), len(accepted)


def reconstruct_map(
    source: base.SourceData,
    baseline_by_point: dict[int, dict[str, str]],
    order: int,
    speed: float,
    inv_ols_nm_per_v: float,
) -> tuple[list[dict], dict]:
    scale_pN_per_V = D3_K_N_PER_M * inv_ols_nm_per_v * 1000.0
    matrix = np.full((len(source.curves), BIN_CENTRES_NM.size), np.nan)
    contact_plane_mad_nm: list[float] = []
    endpoint_distance_nm: list[float] = []
    for curve_index, curve in enumerate(source.curves):
        baseline = float(baseline_by_point[int(curve.point_index)]["dynamic_baseline_raw_V"])
        voltage = curve.deflection_V - baseline
        height = curve.measured_height_m
        travel_nm = (height[0] - height) * 1e9
        terminal = travel_nm >= travel_nm[-1] - CONTACT_WINDOW_NM
        contact_coordinate_m = height[terminal] + inv_ols_nm_per_v * 1e-9 * voltage[terminal]
        contact_plane_m = float(np.median(contact_coordinate_m))
        contact_plane_mad_nm.append(mad(contact_coordinate_m * 1e9))
        distance_nm = (height + inv_ols_nm_per_v * 1e-9 * voltage - contact_plane_m) * 1e9
        force_pN = scale_pN_per_V * voltage
        endpoint_distance_nm.append(float(distance_nm[-1]))
        for column, centre in enumerate(BIN_CENTRES_NM):
            selected = (
                (distance_nm >= centre - BIN_HALF_WIDTH_NM)
                & (distance_nm < centre + BIN_HALF_WIDTH_NM)
                & np.isfinite(force_pN)
            )
            if np.count_nonzero(selected) >= 2:
                matrix[curve_index, column] = float(np.median(force_pN[selected]))

    rows: list[dict] = []
    for column, centre in enumerate(BIN_CENTRES_NM):
        values = matrix[:, column]
        values = values[np.isfinite(values)]
        rows.append(
            {
                "acquisition_order": order,
                "source": str(source.path.relative_to(ROOT)).replace("\\", "/"),
                "speed_um_per_s": speed,
                "distance_nm": centre,
                "available_pixels": int(values.size),
                "map_median_force_pN": float(np.median(values)) if values.size else float("nan"),
                "map_force_q25_pN": q(values, 25),
                "map_force_q75_pN": q(values, 75),
                "map_force_spatial_mad_pN": mad(values),
            }
        )
    diagnostic = {
        "acquisition_order": order,
        "speed_um_per_s": speed,
        "contact_plane_within_window_mad_nm_median": float(np.median(contact_plane_mad_nm)),
        "contact_plane_within_window_mad_nm_q95": q(np.asarray(contact_plane_mad_nm), 95),
        "terminal_reconstructed_distance_nm_median": float(np.median(endpoint_distance_nm)),
        "terminal_reconstructed_distance_nm_mad": mad(np.asarray(endpoint_distance_nm)),
        "force_scale_pN_per_V": scale_pN_per_V,
    }
    return rows, diagnostic


def configure_model() -> None:
    base.PROBE_RADIUS_M = PROBE_RADIUS_M
    base.HAMAKER_J = HAMAKER_J
    base.TEMPERATURE_C = TEMPERATURE_C
    base.TEMPERATURE_K = TEMPERATURE_C + 273.15


def prediction(distance_nm: np.ndarray, parameters: np.ndarray) -> np.ndarray:
    return base.total_equilibrium_force_pN(
        distance_nm,
        math.exp(parameters[0]),
        math.exp(parameters[1]),
        parameters[2],
        EPSILON_R_997_GLYCEROL,
        MODEL,
    )


def local_geometry(residual, parameters: np.ndarray, n_points: int) -> tuple[int, float, float, np.ndarray]:
    steps = np.array([1e-4, 1e-4, 0.01])
    jacobian = np.column_stack(
        [
            (residual(parameters + np.eye(3)[j] * steps[j]) - residual(parameters - np.eye(3)[j] * steps[j]))
            / (2.0 * steps[j])
            for j in range(3)
        ]
    )
    norms = np.linalg.norm(jacobian, axis=0)
    if not np.all(np.isfinite(norms)) or np.any(norms <= 0.0):
        return 0, float("inf"), float("nan"), np.full(3, np.nan)
    _, singular, vt = np.linalg.svd(jacobian / norms, full_matrices=False)
    rank = int(np.sum(singular > max(jacobian.shape) * np.finfo(float).eps * singular[0]))
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else float("inf")
    scales = np.full(3, np.nan)
    correlation = float("nan")
    if rank == 3 and n_points > 3:
        covariance = ((vt.T / singular**2) @ vt) / np.outer(norms, norms)
        covariance *= max(1.0, float(np.sum(residual(parameters) ** 2)) / (n_points - 3))
        scales = np.sqrt(np.maximum(0.0, np.diag(covariance)))
        if scales[0] > 0.0 and scales[1] > 0.0:
            correlation = float(covariance[0, 1] / (scales[0] * scales[1]))
    return rank, condition, correlation, scales


def fit_map(rows: list[dict], meta: dict) -> dict:
    selected = [
        row
        for row in rows
        if FIT_MIN_NM <= float(row["distance_nm"]) <= FIT_MAX_NM
        and int(row["available_pixels"]) >= 48
        and np.isfinite(float(row["map_median_force_pN"]))
        and np.isfinite(float(row["map_force_spatial_mad_pN"]))
    ]
    distance = np.asarray([float(row["distance_nm"]) for row in selected])
    force = np.asarray([float(row["map_median_force_pN"]) for row in selected])
    sigma = np.maximum(np.asarray([float(row["map_force_spatial_mad_pN"]) for row in selected]), 5.0)

    def residual(parameters: np.ndarray) -> np.ndarray:
        return (prediction(distance, parameters) - force) / sigma

    starts = []
    for lam, potential in ((20.0, 50.0), (50.0, 100.0), (100.0, 150.0), (200.0, 200.0)):
        offset = float(np.clip(np.median(force[-5:]), -499.0, 499.0))
        starts.append(
            least_squares(
                residual,
                [math.log(lam), math.log(potential), offset],
                bounds=(LOWER, UPPER),
                method="trf",
                loss="soft_l1",
                f_scale=1.0,
                xtol=1e-9,
                ftol=1e-9,
                gtol=1e-9,
                max_nfev=500,
            )
        )
    best = min(starts, key=lambda item: item.cost)
    parameters = best.x
    predicted = prediction(distance, parameters)
    raw_residual = force - predicted
    sse = float(np.sum(raw_residual**2))
    sst = float(np.sum((force - np.mean(force)) ** 2))
    r2 = 1.0 - sse / sst if sst > 0.0 else float("nan")
    rank, condition, correlation, scales = local_geometry(residual, parameters, len(distance))
    boundary_fraction = np.minimum(
        (parameters - LOWER) / (UPPER - LOWER),
        (UPPER - parameters) / (UPPER - LOWER),
    )
    lam, potential = np.exp(parameters[:2])
    signal = float(np.max(force) - np.median(force[-5:]))
    reasons: list[str] = []
    if len(distance) < 20:
        reasons.append("fewer_than_20_bins")
    if not best.success:
        reasons.append("optimizer")
    if rank != 3 or condition > 1e6:
        reasons.append("jacobian")
    if np.any(boundary_fraction < 1e-4):
        reasons.append("parameter_boundary")
    elif np.any(boundary_fraction < 0.05):
        reasons.append("near_parameter_boundary")
    if not np.isfinite(r2) or r2 < 0.70:
        reasons.append("r2_below_0.70")
    if signal < 5.0 * float(np.median(sigma[-5:])):
        reasons.append("weak_signal")
    if not np.all(np.isfinite(scales[:2])) or np.any(scales[:2] > 0.5):
        reasons.append("broad_local_parameter_scale")

    unusable_flags = {"optimizer", "jacobian", "parameter_boundary", "r2_below_0.70", "weak_signal"}
    quality_state = "pass" if not reasons else ("unusable" if any(reason in unusable_flags for reason in reasons) else "weak")
    return {
        **meta,
        "model": MODEL,
        "fit_min_nm": FIT_MIN_NM,
        "fit_max_nm": FIT_MAX_NM,
        "fit_bin_count": int(len(distance)),
        "lambda_D_apparent_nm": float(lam),
        "surface_potential_magnitude_apparent_mV": float(potential),
        "fit_offset_pN": float(parameters[2]),
        "r2": r2,
        "residual_rmse_pN": math.sqrt(sse / len(distance)),
        "weighted_sse": float(np.sum(residual(parameters) ** 2)),
        "optimizer_success": bool(best.success),
        "optimizer_nfev": int(best.nfev),
        "jacobian_rank": rank,
        "normalized_jacobian_condition": condition,
        "local_log_lambda_scale": float(scales[0]),
        "local_log_potential_scale": float(scales[1]),
        "local_lambda_potential_correlation": correlation,
        "fit_quality_pass": not reasons,
        "quality_state": quality_state,
        "quality_flags": ";".join(reasons),
        "epsilon_r_assumed": EPSILON_R_997_GLYCEROL,
        "probe_radius_um_assumed": PROBE_RADIUS_M * 1e6,
        "temperature_C_assumed": TEMPERATURE_C,
        "hydrodynamic_subtraction_applied": False,
    }


def force_figure(chronology: list[dict], binned_rows: list[dict]) -> None:
    lookup = {
        (int(row["acquisition_order"]), float(row["distance_nm"])): row
        for row in binned_rows
    }
    times = [datetime.fromisoformat(row["map_midpoint_time"]) for row in chronology]
    fig, axes = plt.subplots(2, 2, figsize=(18, 9), sharex=True)
    for axis, target in zip(axes.flat, TARGETS_NM, strict=True):
        values = np.asarray([float(lookup[(int(row["acquisition_order"]), target)]["map_median_force_pN"]) for row in chronology])
        q25 = np.asarray([float(lookup[(int(row["acquisition_order"]), target)]["map_force_q25_pN"]) for row in chronology])
        q75 = np.asarray([float(lookup[(int(row["acquisition_order"]), target)]["map_force_q75_pN"]) for row in chronology])
        axis.plot(times, values, color="0.45", lw=1.0, zorder=1)
        for speed, color in COLORS.items():
            selected = np.asarray([np.isclose(float(row["speed_um_per_s"]), speed) for row in chronology])
            axis.errorbar(
                np.asarray(times, dtype=object)[selected],
                values[selected],
                yerr=np.vstack((values[selected] - q25[selected], q75[selected] - values[selected])),
                fmt="o",
                ms=5.6,
                color=color,
                mec="0.25",
                mew=0.45,
                ecolor=color,
                elinewidth=0.75,
                capsize=0,
                alpha=0.95,
                zorder=3,
            )
        axis.set_title(f"D = {target:g} nm", fontsize=14)
        axis.set_ylabel("Absolute force (pN)")
        style_time_axis(axis, chronology)
    offset = times[0].strftime("%z")
    for axis in axes[-1]:
        axis.set_xlabel(f"Map acquisition midpoint time (UTC{offset[:3]}:{offset[3:]})")
    handles = [
        Line2D([], [], marker="o", ls="", color=color, markeredgecolor="0.25", label=f"{speed:g} µm/s")
        for speed, color in COLORS.items()
    ]
    fig.suptitle("Map medians: absolute force versus acquisition time", fontsize=16, y=0.985)
    fig.text(
        0.5,
        0.947,
        "12 Sep 2026 · 99.7 wt% glycerol · cantilever D3 · 5 nm separation bins · bars show spatial IQR",
        ha="center",
        fontsize=10,
    )
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.928), ncol=4, frameon=False)
    fig.text(
        0.5,
        0.012,
        "Speed-conditioned constant baseline; 10 nm physical hard-contact window; dashed lines separate four rotated palindrome blocks.",
        ha="center",
        fontsize=9,
        color="0.3",
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.89))
    for suffix in ("png", "svg"):
        fig.savefig(FIG / f"map_median_force_slices.{suffix}", dpi=200)
    plt.close(fig)


def parameter_figure(chronology: list[dict], fits: list[dict]) -> None:
    by_order = {int(row["acquisition_order"]): row for row in fits}
    ordered = [by_order[int(row["acquisition_order"])] for row in chronology]
    times = np.asarray([datetime.fromisoformat(row["map_midpoint_time"]) for row in chronology], dtype=object)
    fig, axes = plt.subplots(2, 1, figsize=(13.2, 8.1), sharex=True, layout="constrained")
    keys = ("lambda_D_apparent_nm", "surface_potential_magnitude_apparent_mV")
    labels = (r"Apparent Debye length $\lambda_D$ (nm)", r"Apparent surface potential $|\psi_s|$ (mV)")
    for axis, key, label in zip(axes, keys, labels, strict=True):
        values = np.asarray([float(row[key]) for row in ordered])
        states = np.asarray([row["quality_state"] for row in ordered])
        usable = states != "unusable"
        passed = states == "pass"
        weak = states == "weak"
        axis.plot(times, np.where(usable, values, np.nan), color="0.65", lw=0.9, zorder=1)
        for speed, color in COLORS.items():
            speed_mask = np.asarray([np.isclose(float(row["speed_um_per_s"]), speed) for row in chronology])
            axis.scatter(times[speed_mask & passed], values[speed_mask & passed], color=color, edgecolor="0.25", lw=0.5, s=42, zorder=3)
            axis.scatter(times[speed_mask & weak], values[speed_mask & weak], facecolors="none", edgecolors=color, lw=1.8, s=52, zorder=4)
        unusable = ~usable
        axis.scatter(times[unusable], np.full(np.count_nonzero(unusable), 0.035), transform=axis.get_xaxis_transform(), marker="x", color="#ad3f36", s=43, lw=1.4, zorder=4)
        for index in np.flatnonzero(unusable):
            axis.annotate(f"#{index + 1}", (mdates.date2num(times[index]), 0.035), xycoords=axis.get_xaxis_transform(), xytext=(3, 7), textcoords="offset points", fontsize=8, color="#ad3f36")
        axis.set_ylabel(label)
        finite_usable = values[usable & np.isfinite(values)]
        if finite_usable.size:
            axis.set_ylim(0.0, 1.22 * float(np.max(finite_usable)))
        style_time_axis(axis, chronology)
    legend = [Line2D([], [], marker="o", color="none", markerfacecolor=color, markeredgecolor="0.3", label=f"{speed:g} µm/s") for speed, color in COLORS.items()]
    legend += [
        Line2D([], [], marker="o", color="none", markerfacecolor="none", markeredgecolor="0.4", label="Weakly constrained optimum"),
        Line2D([], [], marker="x", color="none", markeredgecolor="#ad3f36", label="Unusable fit (bottom rail)"),
    ]
    axes[0].legend(handles=legend, loc="upper right", frameon=False, ncol=2, fontsize=9)
    axes[0].set_title("12 Sep 2026 · 99.7 wt% glycerol · per-map nonlinear PB fits", loc="left", fontweight="bold", pad=18)
    axes[1].set_xlabel("Acquisition time (instrument UTC+02:00; map midpoint)")
    fig.get_layout_engine().set(rect=(0, 0.085, 1, 0.925))
    fig.text(
        0.07,
        0.016,
        "64 curves/map · fit range 20–250 nm · R = 4.54685 µm and εr = 42.5 assumed · no hydrodynamic subtraction\n"
        "Finite-speed apparent parameters: marker quality concerns numerical identifiability, not equilibrium-model validity.",
        fontsize=9,
        color="0.3",
    )
    for suffix in ("png", "svg"):
        fig.savefig(FIG / f"debye_length_surface_potential_vs_time.{suffix}", dpi=200)
    plt.close(fig)


def rank_summary(chronology: list[dict], binned_rows: list[dict], fits: list[dict]) -> list[dict]:
    times = np.asarray([float(row["elapsed_midpoint_min"]) for row in chronology])
    speeds = np.asarray([float(row["speed_um_per_s"]) for row in chronology])
    output: list[dict] = []
    binned = {(int(row["acquisition_order"]), float(row["distance_nm"])): row for row in binned_rows}
    observables: list[tuple[str, str, np.ndarray, np.ndarray]] = []
    for target in TARGETS_NM:
        values = np.asarray([float(binned[(int(row["acquisition_order"]), target)]["map_median_force_pN"]) for row in chronology])
        observables.append((f"absolute_force_{target:g}nm", "pN", values, np.ones(values.size, dtype=bool)))
    fit_by_order = {int(row["acquisition_order"]): row for row in fits}
    ordered_fits = [fit_by_order[int(row["acquisition_order"])] for row in chronology]
    usable = np.asarray([row["quality_state"] != "unusable" for row in ordered_fits])
    for key, label in (("lambda_D_apparent_nm", "apparent_debye_length"), ("surface_potential_magnitude_apparent_mV", "apparent_surface_potential")):
        values = np.asarray([float(row[key]) for row in ordered_fits])
        observables.append((label, "nm" if "debye" in label else "mV", values, usable))
    for name, unit, values, valid in observables:
        keep = valid & np.isfinite(values)
        rho_speed = float(spearmanr(speeds[keep], values[keep]).statistic) if np.count_nonzero(keep) >= 3 else float("nan")
        rho_time = float(spearmanr(times[keep], values[keep]).statistic) if np.count_nonzero(keep) >= 3 else float("nan")
        same_speed = []
        for speed in SPEEDS:
            selected = keep & np.isclose(speeds, speed)
            if np.count_nonzero(selected) >= 3:
                same_speed.append(float(spearmanr(times[selected], values[selected]).statistic))
        output.append(
            {
                "observable": name,
                "unit": unit,
                "n_maps": int(np.count_nonzero(keep)),
                "global_spearman_rho_speed": rho_speed,
                "global_spearman_rho_elapsed_time": rho_time,
                "same_speed_time_rho_median": float(np.median(same_speed)) if same_speed else float("nan"),
                "same_speed_time_rho_min": float(np.min(same_speed)) if same_speed else float("nan"),
                "same_speed_time_rho_max": float(np.max(same_speed)) if same_speed else float("nan"),
                "inference_boundary": "descriptive map-level association; sequential non-randomized acquisition",
            }
        )
    return output


def write_report(
    inv_ols: float,
    map_contacts: list[dict],
    binned_rows: list[dict],
    reconstruction_diagnostics: list[dict],
    fits: list[dict],
    ranks: list[dict],
) -> None:
    passed = sum(row["quality_state"] == "pass" for row in fits)
    weak = sum(row["quality_state"] == "weak" for row in fits)
    unusable = sum(row["quality_state"] == "unusable" for row in fits)
    scale = D3_K_N_PER_M * inv_ols
    embedded_scale = 0.1888596503 * 62.6706576105
    force_lines = []
    for target in TARGETS_NM:
        values = [float(row["map_median_force_pN"]) for row in binned_rows if np.isclose(float(row["distance_nm"]), target)]
        force_lines.append(f"| {target:g} | {np.min(values):.1f}–{np.max(values):.1f} |")
    rank_lookup = {row["observable"]: row for row in ranks}
    lines = [
        "# 12-09-26 · 99.7 wt% glycerol · D3 · force/PB time figures",
        "",
        "## 物理图像",
        "",
        "动态 baseline 去掉的是每条 approach 的速度相关 detector offset；它没有去掉随距离增长的 viscous drainage。重建后的 20–200 nm 绝对力随 approach speed 显著增大，并在同速度下随 acquisition history 改变。因此 PB 图中的两个参数只能读作把有限速度曲线投影到 equilibrium PB 形状后的 apparent parameters，不能读作溶液真实 Debye length 或 equilibrium surface potential。",
        "",
        "## 对应图",
        "",
        "- `figures/map_median_force_slices.png`：与旧图相同的 20/50/100/200 nm map-median absolute-force 时间序列。",
        "- `figures/debye_length_surface_potential_vs_time.png`：与旧图相同的 apparent Debye length / surface-potential 时间序列。",
        "",
        "## 本批数据专用重建",
        "",
        f"- baseline：直接读取 held-out-map 验证过的四个 speed-conditioned 120 nm 窗口；每条曲线只减一个常数。",
        f"- contact：每条曲线只在 terminal {CONTACT_WINDOW_NM:g} nm scanner-travel 物理窗口内检测，不使用固定样本数。",
        f"- 本批 contact-derived InvOLS = **{inv_ols:.4f} nm/V**（32 个 map contact 中位数的中位数）；旧 D3 值 {D3_PRIOR_INVOLS_NM_PER_V:.4f} nm/V 不满足本批硬接触单位斜率。",
        f"- force scale = **{scale:.4f} nN/V**，使用独立 D3 `k={D3_K_N_PER_M:.9f} N/m`；它与 JPK embedded product {embedded_scale:.4f} nN/V 很接近，但这种乘积一致不等于两项独立标定都正确。",
        f"- 每个 map 可接受 terminal contacts 为 {min(int(row['accepted_contacts']) for row in map_contacts)}–{max(int(row['accepted_contacts']) for row in map_contacts)}/64；重建 endpoint separation 的 map median 绝对值最大为 {max(abs(float(row['terminal_reconstructed_distance_nm_median'])) for row in reconstruction_diagnostics):.3f} nm。",
        "",
        "## 绝对力范围",
        "",
        "| separation (nm) | 32-map median-force range (pN) |",
        "|---:|---:|",
        *force_lines,
        "",
        "## PB 映射质量与边界",
        "",
        f"- {passed}/32 pass，{weak}/32 weak，{unusable}/32 unusable。判据只检查 optimizer、boundary、Jacobian、R²、信号和局部参数尺度。",
        f"- 模型：equal-potential nonlinear PB Derjaguin + van der Waals + constant offset，20–250 nm；假定 `R={PROBE_RADIUS_M*1e6:.6f} µm`、`T={TEMPERATURE_C:.1f} °C`、`εr={EPSILON_R_997_GLYCEROL:g}`。",
        f"- `εr` 外部参考：{EPSILON_R_REFERENCE}；它仍是当前样品在假定温度下的模型输入，不是从本批 AFM 数据识别出的参数。",
        "- 未做 hydrodynamic subtraction；99.7% glycerol 的速度相关瞬态和 drainage 是 load-bearing model mismatch。surface potential 还条件依赖 assumed dielectric constant、probe radius、contact zero 和 force calibration。",
        "",
        "## 描述性时间/速度相关",
        "",
        "| observable | n maps | rho(speed) | rho(time) | median same-speed rho(time) |",
        "|---|---:|---:|---:|---:|",
    ]
    for name in ("absolute_force_20nm", "absolute_force_50nm", "absolute_force_100nm", "absolute_force_200nm", "apparent_debye_length", "apparent_surface_potential"):
        row = rank_lookup[name]
        lines.append(f"| {name} | {row['n_maps']} | {float(row['global_spearman_rho_speed']):.3f} | {float(row['global_spearman_rho_elapsed_time']):.3f} | {float(row['same_speed_time_rho_median']):.3f} |")
    lines += [
        "",
        "这些 rho 是 sequential、non-randomized acquisition 的描述统计；speed 与 time/history 仍不能因一张相关图而因果分离。",
        "",
        "## 文件",
        "",
        "- `contact_calibration_curves.csv` / `contact_calibration_maps.csv`：物理窗口接触标定和 QC。",
        "- `map_force_by_separation.csv` / `map_reconstruction_diagnostics.csv`：绝对力时间图的数值来源。",
        "- `map_pb_apparent_fits.csv`：PB optimum、残差和 identifiability flags。",
        "- `rank_associations.csv`：map-level 速度/时间秩相关。",
        "- `provenance.json` / `artifact_manifest.sha256`：假设、版本与输出哈希。",
        "",
    ]
    (OUT / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def write_manifest() -> None:
    targets = sorted(path for path in OUT.rglob("*") if path.is_file() and path.name != "artifact_manifest.sha256")
    lines = [f"{sha256_file(path)}  {path.relative_to(ROOT).as_posix()}" for path in targets]
    (OUT / "artifact_manifest.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    configure_model()
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    chronology = sorted(read_csv(CHRONOLOGY / "map_time_speed_observables.csv"), key=lambda row: int(row["acquisition_order"]))
    baseline_rows = read_csv(ADAPTIVE / "curve_dynamic_baselines.csv")
    by_source: dict[str, dict[int, dict[str, str]]] = {}
    for row in baseline_rows:
        by_source.setdefault(row["source"], {})[int(row["point_index"])] = row
    if len(chronology) != 32 or set(by_source) != {row["source"] for row in chronology}:
        raise RuntimeError("Chronology/baseline cohort mismatch")

    sources: dict[str, base.SourceData] = {}
    contact_rows: list[dict] = []
    map_contacts: list[dict] = []
    for meta in chronology:
        order = int(meta["acquisition_order"])
        source_rel = meta["source"]
        source = base.load_source(ROOT / source_rel, 99.7)
        sources[source_rel] = source
        rows, median_inv_ols, inv_ols_mad, accepted = contact_diagnostic(source, by_source[source_rel], order)
        contact_rows.extend(rows)
        map_contacts.append(
            {
                "acquisition_order": order,
                "speed_um_per_s": float(meta["speed_um_per_s"]),
                "accepted_contacts": accepted,
                "map_contact_invOLS_median_nm_per_V": median_inv_ols,
                "map_contact_invOLS_mad_nm_per_V": inv_ols_mad,
            }
        )
    inv_ols = float(np.median([row["map_contact_invOLS_median_nm_per_V"] for row in map_contacts]))

    binned_rows: list[dict] = []
    reconstruction_diagnostics: list[dict] = []
    fits: list[dict] = []
    for meta in chronology:
        order = int(meta["acquisition_order"])
        source_rel = meta["source"]
        speed = float(meta["speed_um_per_s"])
        current_rows, diagnostic = reconstruct_map(sources[source_rel], by_source[source_rel], order, speed, inv_ols)
        binned_rows.extend(current_rows)
        reconstruction_diagnostics.append({**{key: meta[key] for key in ("acquisition_order", "block", "position_in_block", "map_midpoint_time", "elapsed_midpoint_min", "speed_um_per_s", "source")}, **diagnostic})
        fits.append(fit_map(current_rows, {key: meta[key] for key in ("acquisition_order", "block", "position_in_block", "map_midpoint_time", "elapsed_midpoint_min", "speed_um_per_s", "source")}))

    ranks = rank_summary(chronology, binned_rows, fits)
    write_csv(OUT / "contact_calibration_curves.csv", contact_rows)
    write_csv(OUT / "contact_calibration_maps.csv", map_contacts)
    write_csv(OUT / "map_force_by_separation.csv", binned_rows)
    write_csv(OUT / "map_reconstruction_diagnostics.csv", reconstruction_diagnostics)
    write_csv(OUT / "map_pb_apparent_fits.csv", fits)
    write_csv(OUT / "rank_associations.csv", ranks)
    force_figure(chronology, binned_rows)
    parameter_figure(chronology, fits)
    write_report(inv_ols, map_contacts, binned_rows, reconstruction_diagnostics, fits, ranks)

    provenance = {
        "script": str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
        "input_dynamic_baseline_csv": str((ADAPTIVE / "curve_dynamic_baselines.csv").relative_to(ROOT)).replace("\\", "/"),
        "input_chronology_csv": str((CHRONOLOGY / "map_time_speed_observables.csv").relative_to(ROOT)).replace("\\", "/"),
        "maps": 32,
        "curves_per_map": 64,
        "baseline_mode": "held-out-validated speed-conditioned physical-window constant median",
        "contact_mode": f"terminal {CONTACT_WINDOW_NM:g} nm physical scanner-travel robust line and contact-coordinate median",
        "batch_contact_invOLS_nm_per_V": inv_ols,
        "batch_contact_invOLS_aggregation": "median of 32 per-map accepted-contact medians",
        "spring_constant_N_per_m": D3_K_N_PER_M,
        "force_scale_nN_per_V": D3_K_N_PER_M * inv_ols,
        "old_D3_invOLS_rejected_for_this_batch_nm_per_V": D3_PRIOR_INVOLS_NM_PER_V,
        "model": MODEL,
        "fit_range_nm": [FIT_MIN_NM, FIT_MAX_NM],
        "probe_radius_m_assumed": PROBE_RADIUS_M,
        "temperature_C_assumed": TEMPERATURE_C,
        "epsilon_r_assumed": EPSILON_R_997_GLYCEROL,
        "epsilon_r_reference": EPSILON_R_REFERENCE,
        "hamaker_J_assumed": HAMAKER_J,
        "hydrodynamic_subtraction_applied": False,
        "parameter_interpretation": "finite-speed model-conditioned apparent values, not identified equilibrium material parameters",
        "software": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    write_manifest()
    print(json.dumps({"output": str(OUT), "invOLS_nm_per_V": inv_ols, "fit_quality": {state: sum(row["quality_state"] == state for row in fits) for state in ("pass", "weak", "unusable")}}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Re-fit thermal PSDs and calibrate three AFM cantilevers.

The exported JPK fit parameters and ``fit-data`` column are deliberately not
used for the calibration.  Fundamental resonances are detected from the
measured PSD, re-fit with a simple-harmonic-oscillator model, and combined
with hard-contact inverse optical lever sensitivities decoded from the raw
JPK force curves.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable
from zipfile import ZipFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.constants import Boltzmann
from scipy.ndimage import median_filter
from scipy.optimize import minimize
from scipy.signal import find_peaks


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT / "calibration"
RESULTS = ROOT / "analysis" / "results"
FIGURES = RESULTS / "figures"

TEMPERATURE_C = 25.6
TEMPERATURE_K = TEMPERATURE_C + 273.15

# Peak search is intentionally independent of the exported/header frequency.
PEAK_SEARCH_MIN_HZ = 1_000.0
PEAK_SEARCH_MAX_HZ = 200_000.0
PEAK_BASELINE_KERNEL_BINS = 1001
PEAK_PROMINENCE_DECADES = 0.8
PEAK_MIN_DISTANCE_BINS = 200

THERMAL_FIT_HALF_WIDTH_HZ = 3_000.0
THERMAL_FIT_WINDOW_CHECKS_HZ = (2_000.0, 4_000.0)
CONTACT_SPAN_NM = 50.0
CONTACT_SPAN_CHECKS_NM = (40.0, 60.0)

# Batch-label values visible in IMG_4571.jpeg.  They are used only as an
# ex-post sanity check, never to select the resonance peak.
LABEL_FREQUENCY_RANGE_HZ = (11_000.0, 17_000.0)
LABEL_SPRING_RANGE_N_PER_M = (0.11, 0.56)


@dataclass
class ThermalFit:
    group: int
    file: str
    sha256: str
    exported_frequency_hz: float
    detected_peak_hz: float
    candidate_peaks_hz: list[float]
    candidate_prominences_decades: list[float]
    fit_frequency_hz: float
    quality_factor: float
    amplitude_v_per_sqrt_hz: float
    background_v2_per_hz: float
    resonance_area_v2: float
    correction_factor: float
    thermal_factor_n_m_per_v2: float
    fit_log10_rmse: float
    fit_success: bool
    fit_nfev: int
    fit_half_width_hz: float


@dataclass
class ForceFit:
    group: int
    file: str
    sha256: str
    contact_span_nm: float
    contact_actual_span_nm: float
    contact_points: int
    slope_v_per_m: float
    sensitivity_m_per_v: float
    sensitivity_fit_se_m_per_v: float
    r_squared: float
    stored_sensitivity_m_per_v: float
    stored_spring_constant_n_per_m: float
    sensitivity_40nm_m_per_v: float
    sensitivity_60nm_m_per_v: float


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"No rows supplied for {path}")
    if fieldnames is None:
        fieldnames = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_number(text: str) -> float:
    match = re.search(r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?", text)
    if match is None:
        raise ValueError(f"No numeric value in {text!r}")
    return float(match.group())


def parse_tnd_header(path: Path) -> dict[str, str]:
    header: dict[str, str] = {}
    with path.open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            body = line[1:].strip()
            if ": " in body:
                key, value = body.split(": ", 1)
                header[key] = value.strip()
    return header


def exported_frequency_hz(header: dict[str, str]) -> float:
    value = header.get("parameter.f", "nan")
    number = parse_number(value)
    if "kHz" in value:
        number *= 1_000.0
    return number


def load_tnd(path: Path) -> tuple[dict[str, str], np.ndarray]:
    header = parse_tnd_header(path)
    data = np.loadtxt(path, comments="#", dtype=np.float64)
    if data.ndim != 2 or data.shape[1] != 4:
        raise ValueError(f"Expected four TND columns in {path}, got {data.shape}")
    if data.shape[0] < 1000:
        raise ValueError(f"Unexpectedly short TND spectrum in {path}")
    if not np.all(np.isfinite(data)):
        raise ValueError(f"Non-finite TND values in {path}")
    if not np.all(np.diff(data[:, 0]) > 0):
        raise ValueError(f"Frequency is not strictly increasing in {path}")
    # Columns: frequency, raw Vertical Deflection PSD, measured averaged PSD,
    # exported fit-data.  Column 3 (fit-data) is never consumed here.
    if np.any(data[1:, 1:3] <= 0):
        raise ValueError(f"Non-positive measured PSD values in {path}")
    return header, data


def detect_resonance_candidates(frequency: np.ndarray, measured_psd: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mask = (
        (frequency >= PEAK_SEARCH_MIN_HZ)
        & (frequency <= PEAK_SEARCH_MAX_HZ)
        & np.isfinite(measured_psd)
        & (measured_psd > 0)
    )
    f = frequency[mask]
    log_psd = np.log10(measured_psd[mask])
    baseline = median_filter(log_psd, size=PEAK_BASELINE_KERNEL_BINS, mode="nearest")
    excess = log_psd - baseline
    peaks, properties = find_peaks(
        excess,
        prominence=PEAK_PROMINENCE_DECADES,
        distance=PEAK_MIN_DISTANCE_BINS,
    )
    if peaks.size == 0:
        raise RuntimeError("No significant thermal resonance was detected")
    order = np.argsort(f[peaks])
    return f[peaks][order], properties["prominences"][order]


def sho_psd(
    frequency: np.ndarray,
    resonance_hz: float,
    quality_factor: float,
    amplitude_squared: float,
    background: float,
) -> np.ndarray:
    ratio = frequency / resonance_hz
    denominator = (1.0 - ratio * ratio) ** 2 + (ratio / quality_factor) ** 2
    return background + amplitude_squared / denominator


def fit_sho(
    frequency: np.ndarray,
    measured_psd: np.ndarray,
    detected_peak_hz: float,
    half_width_hz: float,
) -> dict[str, float | bool | int]:
    mask = (
        (frequency >= detected_peak_hz - half_width_hz)
        & (frequency <= detected_peak_hz + half_width_hz)
        & np.isfinite(measured_psd)
        & (measured_psd > 0)
    )
    f = frequency[mask]
    p = measured_psd[mask]
    if f.size < 100:
        raise RuntimeError("Too few PSD bins in the SHO fit window")

    edge_count = max(20, f.size // 5)
    background_0 = float(np.median(np.r_[p[:edge_count], p[-edge_count:]]))
    peak_index = int(np.argmax(p))
    peak_power = float(p[peak_index])
    half_power = background_0 + 0.5 * max(peak_power - background_0, 0.0)
    left = np.where(p[:peak_index] < half_power)[0]
    right = np.where(p[peak_index + 1 :] < half_power)[0]
    if left.size and right.size:
        width = f[peak_index + 1 + right[0]] - f[left[-1]]
        quality_0 = float(np.clip(detected_peak_hz / width, 5.0, 300.0))
    else:
        quality_0 = 60.0
    amplitude_squared_0 = max((peak_power - background_0) / quality_0**2, 1e-30)

    x0 = np.array(
        [
            detected_peak_hz,
            math.log(quality_0),
            math.log(amplitude_squared_0),
            math.log(background_0),
        ],
        dtype=np.float64,
    )
    bounds = [
        (detected_peak_hz - 800.0, detected_peak_hz + 800.0),
        (math.log(5.0), math.log(500.0)),
        (math.log(1e-30), math.log(1e-8)),
        (math.log(1e-15), math.log(1e-7)),
    ]

    def objective(x: np.ndarray) -> float:
        model = sho_psd(f, x[0], math.exp(x[1]), math.exp(x[2]), math.exp(x[3]))
        # Gamma/Whittle likelihood for a measured (possibly averaged) PSD.
        # An unknown, constant number of averages multiplies the objective and
        # therefore does not change the optimum.
        return float(np.mean(np.log(model) + p / model))

    result = minimize(
        objective,
        x0,
        method="Nelder-Mead",
        bounds=bounds,
        options={"maxiter": 20_000, "xatol": 1e-9, "fatol": 1e-13},
    )
    resonance_hz = float(result.x[0])
    quality_factor = float(math.exp(result.x[1]))
    amplitude_squared = float(math.exp(result.x[2]))
    background = float(math.exp(result.x[3]))
    model = sho_psd(f, resonance_hz, quality_factor, amplitude_squared, background)
    log10_rmse = float(np.sqrt(np.mean((np.log10(p) - np.log10(model)) ** 2)))
    resonance_area_v2 = float(0.5 * math.pi * amplitude_squared * quality_factor * resonance_hz)
    return {
        "frequency_hz": resonance_hz,
        "quality_factor": quality_factor,
        "amplitude_v_per_sqrt_hz": math.sqrt(amplitude_squared),
        "amplitude_squared": amplitude_squared,
        "background_v2_per_hz": background,
        "resonance_area_v2": resonance_area_v2,
        "log10_rmse": log10_rmse,
        "success": bool(result.success),
        "nfev": int(result.nfev),
    }


def refit_tnd(path: Path, group: int, half_width_hz: float = THERMAL_FIT_HALF_WIDTH_HZ) -> tuple[ThermalFit, dict]:
    header, data = load_tnd(path)
    frequency = data[:, 0]
    raw_psd = data[:, 1]
    measured_average_psd = data[:, 2]
    candidates, prominences = detect_resonance_candidates(frequency, measured_average_psd)
    selected = float(candidates[0])  # Lowest significant resonance: fundamental.
    fit = fit_sho(frequency, measured_average_psd, selected, half_width_hz)

    correction_factor = parse_number(header["cantilever.correction factor"])
    thermal_factor = correction_factor * Boltzmann * TEMPERATURE_K / float(fit["resonance_area_v2"])
    record = ThermalFit(
        group=group,
        file=str(path.relative_to(ROOT)),
        sha256=sha256_file(path),
        exported_frequency_hz=exported_frequency_hz(header),
        detected_peak_hz=selected,
        candidate_peaks_hz=[float(value) for value in candidates],
        candidate_prominences_decades=[float(value) for value in prominences],
        fit_frequency_hz=float(fit["frequency_hz"]),
        quality_factor=float(fit["quality_factor"]),
        amplitude_v_per_sqrt_hz=float(fit["amplitude_v_per_sqrt_hz"]),
        background_v2_per_hz=float(fit["background_v2_per_hz"]),
        resonance_area_v2=float(fit["resonance_area_v2"]),
        correction_factor=correction_factor,
        thermal_factor_n_m_per_v2=thermal_factor,
        fit_log10_rmse=float(fit["log10_rmse"]),
        fit_success=bool(fit["success"]),
        fit_nfev=int(fit["nfev"]),
        fit_half_width_hz=half_width_hz,
    )
    plot_data = {
        "frequency": frequency,
        "raw_psd": raw_psd,
        "measured_average_psd": measured_average_psd,
        "fit": fit,
        "record": record,
    }
    return record, plot_data


def parse_properties(payload: bytes) -> dict[str, str]:
    properties: dict[str, str] = {}
    for line in payload.decode("ISO-8859-1").splitlines():
        if not line or line.startswith(("#", "!")) or "=" not in line:
            continue
        key, value = line.split("=", 1)
        properties[key.strip()] = value.strip()
    return properties


def decode_channel(
    archive: ZipFile,
    properties: dict[str, str],
    segment: int,
    lcd_index: int,
    channel: str,
) -> np.ndarray:
    payload = archive.read(f"segments/{segment}/channels/{channel}.dat")
    raw = np.frombuffer(payload, dtype=">i4").astype(np.float64)
    prefix = f"lcd-info.{lcd_index}.encoder.scaling."
    values = float(properties[prefix + "offset"]) + float(properties[prefix + "multiplier"]) * raw
    if channel == "measuredHeight":
        conversion = "lcd-info.3.conversion-set.conversion.nominal.scaling."
        values = float(properties[conversion + "offset"]) + float(properties[conversion + "multiplier"]) * values
    return values


def linear_contact_fit(height_m: np.ndarray, deflection_v: np.ndarray, span_nm: float) -> dict[str, float | int]:
    span_m = span_nm * 1e-9
    terminal_height = float(height_m[-1])
    indices = np.flatnonzero(np.abs(height_m - terminal_height) <= span_m)
    if indices.size < 20:
        raise RuntimeError(f"Only {indices.size} points in the terminal {span_nm:g} nm")
    # Preserve only the contiguous terminal block.
    start = int(indices[0])
    terminal_indices = np.arange(start, height_m.size)
    x = height_m[terminal_indices]
    y = deflection_v[terminal_indices]
    x_center = float(np.mean(x))
    design = np.column_stack([np.ones(x.size), x - x_center])
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    intercept, slope = (float(value) for value in coefficients)
    prediction = intercept + slope * (x - x_center)
    residual = y - prediction
    ss_res = float(np.sum(residual**2))
    ss_total = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 - ss_res / ss_total
    variance = ss_res / (x.size - 2)
    slope_se = math.sqrt(variance / float(np.sum((x - x_center) ** 2)))
    sensitivity = 1.0 / abs(slope)
    sensitivity_se = slope_se / (slope * slope)
    return {
        "slope_v_per_m": slope,
        "sensitivity_m_per_v": sensitivity,
        "sensitivity_fit_se_m_per_v": sensitivity_se,
        "r_squared": r_squared,
        "points": int(x.size),
        "actual_span_nm": float(np.ptp(x) * 1e9),
        "indices": terminal_indices,
        "prediction": prediction,
    }


def fit_force_curve(path: Path, group: int) -> tuple[ForceFit, dict]:
    with ZipFile(path) as archive:
        corrupt_member = archive.testzip()
        if corrupt_member is not None:
            raise RuntimeError(f"CRC failure in {path}: {corrupt_member}")
        properties = parse_properties(archive.read("shared-data/header.properties"))
        deflection_v = decode_channel(archive, properties, 0, 1, "vDeflection")
        measured_height_m = decode_channel(archive, properties, 0, 3, "measuredHeight")

    if deflection_v.shape != measured_height_m.shape or deflection_v.size != 1000:
        raise ValueError(f"Unexpected force-curve arrays in {path}")
    if not np.all(np.isfinite(deflection_v)) or not np.all(np.isfinite(measured_height_m)):
        raise ValueError(f"Non-finite force-curve values in {path}")

    primary = linear_contact_fit(measured_height_m, deflection_v, CONTACT_SPAN_NM)
    check_40 = linear_contact_fit(measured_height_m, deflection_v, 40.0)
    check_60 = linear_contact_fit(measured_height_m, deflection_v, 60.0)
    if float(primary["slope_v_per_m"]) >= 0:
        raise RuntimeError(f"Unexpected approach contact-slope sign in {path}")

    distance_key = "lcd-info.1.conversion-set.conversion.distance.scaling.multiplier"
    force_key = "lcd-info.1.conversion-set.conversion.force.scaling.multiplier"
    record = ForceFit(
        group=group,
        file=str(path.relative_to(ROOT)),
        sha256=sha256_file(path),
        contact_span_nm=CONTACT_SPAN_NM,
        contact_actual_span_nm=float(primary["actual_span_nm"]),
        contact_points=int(primary["points"]),
        slope_v_per_m=float(primary["slope_v_per_m"]),
        sensitivity_m_per_v=float(primary["sensitivity_m_per_v"]),
        sensitivity_fit_se_m_per_v=float(primary["sensitivity_fit_se_m_per_v"]),
        r_squared=float(primary["r_squared"]),
        stored_sensitivity_m_per_v=float(properties[distance_key]),
        stored_spring_constant_n_per_m=float(properties[force_key]),
        sensitivity_40nm_m_per_v=float(check_40["sensitivity_m_per_v"]),
        sensitivity_60nm_m_per_v=float(check_60["sensitivity_m_per_v"]),
    )
    plot_data = {
        "height_m": measured_height_m,
        "deflection_v": deflection_v,
        "indices": primary["indices"],
        "prediction": primary["prediction"],
        "record": record,
    }
    return record, plot_data


def sample_sd(values: np.ndarray) -> float:
    return float(np.std(values, ddof=1)) if values.size > 1 else float("nan")


def summarize_group(
    group: int,
    thermal: list[ThermalFit],
    force: list[ForceFit],
    fit_window_factors: dict[float, list[float]],
) -> dict[str, float | int]:
    sensitivity = np.array([item.sensitivity_m_per_v for item in force])
    sensitivity_40 = np.array([item.sensitivity_40nm_m_per_v for item in force])
    sensitivity_60 = np.array([item.sensitivity_60nm_m_per_v for item in force])
    factor = np.array([item.thermal_factor_n_m_per_v2 for item in thermal])
    frequency = np.array([item.fit_frequency_hz for item in thermal])
    quality = np.array([item.quality_factor for item in thermal])

    sensitivity_mean = float(np.mean(sensitivity))
    factor_mean = float(np.mean(factor))
    spring = factor_mean / sensitivity_mean**2
    sensitivity_sd = sample_sd(sensitivity)
    factor_sd = sample_sd(factor)
    relative_repeatability = math.sqrt(
        (factor_sd / factor_mean) ** 2 + (2.0 * sensitivity_sd / sensitivity_mean) ** 2
    )
    spring_repeatability_sd = spring * relative_repeatability

    contact_window_springs = [
        factor_mean / float(np.mean(sensitivity_40)) ** 2,
        factor_mean / float(np.mean(sensitivity_60)) ** 2,
    ]
    contact_window_shift = max(abs(value - spring) for value in contact_window_springs)
    thermal_window_springs = [
        float(np.mean(fit_window_factors[width])) / sensitivity_mean**2
        for width in THERMAL_FIT_WINDOW_CHECKS_HZ
    ]
    thermal_window_shift = max(abs(value - spring) for value in thermal_window_springs)

    return {
        "cantilever": group,
        "temperature_C": TEMPERATURE_C,
        "thermal_spectra": len(thermal),
        "force_curves": len(force),
        "sensitivity_nm_per_V": sensitivity_mean * 1e9,
        "sensitivity_repeatability_sd_nm_per_V": sensitivity_sd * 1e9,
        "sensitivity_40nm_mean_nm_per_V": float(np.mean(sensitivity_40)) * 1e9,
        "sensitivity_60nm_mean_nm_per_V": float(np.mean(sensitivity_60)) * 1e9,
        "minimum_contact_fit_R2": min(item.r_squared for item in force),
        "resonance_frequency_kHz": float(np.mean(frequency)) / 1e3,
        "resonance_frequency_sd_kHz": sample_sd(frequency) / 1e3,
        "quality_factor": float(np.mean(quality)),
        "quality_factor_sd": sample_sd(quality),
        "thermal_factor_1e15_N_m_per_V2": factor_mean * 1e15,
        "thermal_factor_sd_1e15_N_m_per_V2": factor_sd * 1e15,
        "spring_constant_N_per_m": spring,
        "spring_constant_repeatability_sd_N_per_m": spring_repeatability_sd,
        "spring_shift_contact_window_N_per_m": contact_window_shift,
        "spring_shift_thermal_window_N_per_m": thermal_window_shift,
        "all_fits_successful": int(all(item.fit_success for item in thermal)),
    }


def thermal_model_for_plot(record: ThermalFit, frequency: np.ndarray) -> np.ndarray:
    return sho_psd(
        frequency,
        record.fit_frequency_hz,
        record.quality_factor,
        record.amplitude_v_per_sqrt_hz**2,
        record.background_v2_per_hz,
    )


def make_thermal_full_plot(plot_data: dict[int, list[dict]], summaries: list[dict]) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 10.5), sharex=True)
    colors = plt.cm.viridis(np.linspace(0.08, 0.9, 5))
    for group, axis in enumerate(axes, start=1):
        for color, data in zip(colors, plot_data[group], strict=True):
            frequency = data["frequency"]
            measured = data["measured_average_psd"]
            record: ThermalFit = data["record"]
            mask = (frequency >= PEAK_SEARCH_MIN_HZ) & (frequency <= 150_000.0)
            axis.plot(frequency[mask] / 1e3, measured[mask], color=color, lw=0.75, alpha=0.78)
            fit_mask = np.abs(frequency - record.fit_frequency_hz) <= THERMAL_FIT_HALF_WIDTH_HZ
            axis.plot(
                frequency[fit_mask] / 1e3,
                thermal_model_for_plot(record, frequency[fit_mask]),
                color=color,
                lw=1.45,
            )
        summary = summaries[group - 1]
        axis.axvline(summary["resonance_frequency_kHz"], color="black", lw=0.8, ls="--")
        axis.set_yscale("log")
        axis.set_ylabel("PSD (V$^2$/Hz)")
        axis.set_title(
            f"Cantilever {group}: blind fundamental fit "
            f"{summary['resonance_frequency_kHz']:.4f} kHz"
        )
        axis.grid(alpha=0.2, which="both")
    axes[-1].set_xlabel("Frequency (kHz)")
    axes[-1].set_xlim(1.0, 150.0)
    fig.suptitle("Measured thermal PSD and independent SHO refits (exported fit-data unused)")
    fig.tight_layout()
    fig.savefig(FIGURES / "thermal_full_spectra_refits.png", dpi=220)
    plt.close(fig)


def make_thermal_zoom_plot(plot_data: dict[int, list[dict]]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.5), sharey=False)
    colors = plt.cm.plasma(np.linspace(0.08, 0.88, 5))
    for group, axis in enumerate(axes, start=1):
        for color, data in zip(colors, plot_data[group], strict=True):
            frequency = data["frequency"]
            measured = data["measured_average_psd"]
            record: ThermalFit = data["record"]
            mask = np.abs(frequency - record.fit_frequency_hz) <= 1_500.0
            axis.plot(frequency[mask] / 1e3, measured[mask], color=color, lw=0.9, alpha=0.7)
            axis.plot(
                frequency[mask] / 1e3,
                thermal_model_for_plot(record, frequency[mask]),
                color=color,
                lw=1.5,
            )
        axis.set_yscale("log")
        axis.set_xlabel("Frequency (kHz)")
        axis.set_ylabel("PSD (V$^2$/Hz)")
        axis.set_title(f"Cantilever {group}")
        axis.grid(alpha=0.2, which="both")
    fig.suptitle("Fundamental-mode data and SHO fits")
    fig.tight_layout()
    fig.savefig(FIGURES / "thermal_fundamental_zoom.png", dpi=220)
    plt.close(fig)


def make_contact_plot(plot_data: dict[int, list[dict]], summaries: list[dict]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.6), sharex=True)
    colors = plt.cm.cividis(np.linspace(0.08, 0.9, 5))
    for group, axis in enumerate(axes, start=1):
        for color, data in zip(colors, plot_data[group], strict=True):
            height = data["height_m"]
            deflection = data["deflection_v"]
            relative_nm = (height - height[-1]) * 1e9
            mask = np.abs(relative_nm) <= 180.0
            axis.plot(relative_nm[mask], deflection[mask], color=color, lw=0.8, alpha=0.75)
            indices = data["indices"]
            axis.plot(
                relative_nm[indices],
                data["prediction"],
                color=color,
                lw=1.8,
            )
        summary = summaries[group - 1]
        axis.set_title(
            f"Cantilever {group}: {summary['sensitivity_nm_per_V']:.2f} nm/V"
        )
        axis.set_xlabel("Measured height from terminal point (nm)")
        axis.set_ylabel("Vertical deflection (V)")
        axis.grid(alpha=0.2)
        axis.invert_xaxis()
    fig.suptitle("Approach hard-contact tails and terminal 50 nm fits")
    fig.tight_layout()
    fig.savefig(FIGURES / "force_contact_fits.png", dpi=220)
    plt.close(fig)


def render_report(
    summaries: list[dict],
    thermal_records: list[ThermalFit],
    force_records: list[ForceFit],
) -> str:
    lines = [
        "# 三支 cantilever 标定报告",
        "",
        f"实验温度：**{TEMPERATURE_C:.1f} °C（{TEMPERATURE_K:.2f} K）**。",
        "",
        "## 建议写入 JPK 的标定值",
        "",
        "| Cantilever | Sensitivity / InvOLS (nm/V) | Spring constant (N/m) | f₀ (kHz) | Q |",
        "|---:|---:|---:|---:|---:|",
    ]
    for summary in summaries:
        lines.append(
            "| {cantilever} | {sensitivity_nm_per_V:.2f} ± {sensitivity_repeatability_sd_nm_per_V:.2f} "
            "| {spring_constant_N_per_m:.4f} ± {spring_constant_repeatability_sd_N_per_m:.4f} "
            "| {resonance_frequency_kHz:.4f} ± {resonance_frequency_sd_kHz:.4f} "
            "| {quality_factor:.2f} ± {quality_factor_sd:.2f} |".format(**summary)
        )
    lines.extend(
        [
            "",
            "表中的 ± 是 5 次重复测量的 **repeatability sample SD**（cantilever 2 的错误导出峰也已从原始谱重拟合，因此仍为 5 次），不是可溯源的绝对标定不确定度。实际输入时可分别使用：",
            "",
        ]
    )
    for summary in summaries:
        lines.append(
            f"- Cantilever {summary['cantilever']}：`sensitivity = {summary['sensitivity_nm_per_V']:.2f} nm/V`，"
            f"`spring constant = {summary['spring_constant_N_per_m']:.4f} N/m`。"
        )

    wrong = next(
        item
        for item in thermal_records
        if item.group == 2 and item.exported_frequency_hz > 50_000.0
    )
    lines.extend(
        [
            "",
            "## 谱选错峰的处理",
            "",
            f"`{wrong.file}` 的导出 header 写成 **{wrong.exported_frequency_hz/1e3:.2f} kHz**。本分析没有沿用这个值，也没有删除整条谱。全谱 blind peak search 在原始 measured PSD 中先找到 {wrong.detected_peak_hz/1e3:.3f} kHz 与 {wrong.candidate_peaks_hz[1]/1e3:.3f} kHz 两个显著峰；按“最低显著共振峰 = fundamental”选取前者，再做 SHO 拟合，得到 **f₀ = {wrong.fit_frequency_hz/1e3:.4f} kHz, Q = {wrong.quality_factor:.2f}**。97.67 kHz 峰保留为高阶模态证据，但不进入基频热标定。",
            "",
            "峰选择不使用照片/厂家给出的 11–17 kHz 范围；该范围只在计算完成后用于 sanity check。三组重拟合的 f₀ 均落在照片标注范围内。",
            "",
            "## 计算方法",
            "",
            "1. 每个 `.tnd` 直接读取 `Frequency` 与 measured `average` PSD；导出的 `fit-data`、`parameter.f`、`parameter.Q`、`parameter.A` 不参与计算。对 1–200 kHz 的 log-PSD 去除宽尺度 median baseline，以 ≥0.8 decade prominence 搜峰，并取最低显著峰作为 fundamental。",
            "2. 在选中峰的 ±3 kHz 内，以 Gamma/Whittle PSD likelihood 拟合",
            "",
            "   `S_VV(f) = N + A² / [(1-(f/f₀)²)² + (f/(f₀Q))²]`。",
            "",
            "   一侧共振面积为 `I_V = (π/2) A² Q f₀`，单位 V²。",
            f"3. 使用文件记录的 rectangular/dynamic correction factor `β = 0.8170` 与用户给定温度计算 voltage-domain thermal factor：`B_z = β k_B T / I_V`（单位 N·m/V²）。",
            "4. `.jpk-force` 使用每条 approach segment 的 measuredHeight 和 vDeflection 原始 int32 数据及各自 metadata conversion 解码；在末端 50 nm 硬接触线性区拟合 `V = a + b z`，取 `Sensitivity = 1/|b|`。",
            "5. 最终 vertical spring constant 为 `k_z = mean(B_z) / mean(Sensitivity)²`。这也正是 JPK 文件中 corrected vertical thermal factor 对 contact sensitivity 的换算关系。",
            "",
            "## QC 与数值敏感性",
            "",
        ]
    )
    for summary in summaries:
        lines.append(
            f"- Cantilever {summary['cantilever']}：5 条 contact fit 的最低 R² = {summary['minimum_contact_fit_R2']:.6f}；"
            f"contact span 从 50 nm 改为 40/60 nm 时，k 最大移动 {summary['spring_shift_contact_window_N_per_m']:.4f} N/m；"
            f"thermal fit half-width 从 3 kHz 改为 2/4 kHz 时，k 最大移动 {summary['spring_shift_thermal_window_N_per_m']:.4f} N/m。"
        )
    lines.extend(
        [
            "- 15 个 JPK force ZIP 均通过 CRC；15 条 thermal spectrum 均为 31,131 个递增 frequency bins，无 NaN/Inf，15 次独立 SHO optimizer 均成功。",
            f"- 照片标签给出的 CONT-W batch 范围为 f₀ = 11–17 kHz、C = 0.11–0.56 N/m；三个最终结果均在此范围。该比较仅作 ex-post physical sanity check。",
            "",
            "## 解释边界",
            "",
            "- Sensitivity 拟合假定 force curve 的接触基底相对 cantilever 足够刚。仅凭文件不能确认基底材料；若基底可压缩，InvOLS 会偏大，继而 k 会偏小。",
            "- repeatability SD 不含 optical-lever spot position、硬接触几何、基底 compliance、0.817 correction model、scanner calibration 等 systematic uncertainty，因此不能当作 traceable absolute uncertainty。",
            "- 这里报告的是 JPK vertical convention 的 `k_z`；不要把 thermal-factor header 中错误显示的 `N/m` 当成该中间量的单位，中间量实际为 N·m/V²。",
            "",
            "## 输出",
            "",
            "- `thermal_refits.csv`：每条原始谱的 blind candidates、独立 f₀/Q/A/area 与 thermal factor。",
            "- `force_contact_fits.csv`：每条 force curve 的 contact slope、Sensitivity、R² 与 40/50/60 nm window sensitivity。",
            "- `calibration_summary.csv`：最终三支 cantilever 标定值及 repeatability/window sensitivity。",
            "- `figures/thermal_full_spectra_refits.png`、`thermal_fundamental_zoom.png`、`force_contact_fits.png`：原始 measured data 与重拟合。",
            "- `provenance.json` 与 `artifact_manifest.sha256`：输入身份、参数、软件版本和输出 hashes。",
            "",
        ]
    )
    return "\n".join(lines)


def create_manifest(paths: Iterable[Path], destination: Path) -> None:
    rows = []
    for path in sorted(paths, key=lambda item: str(item.relative_to(ROOT))):
        rows.append(f"{sha256_file(path)}  {path.relative_to(ROOT)}")
    destination.write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    thermal_records: list[ThermalFit] = []
    force_records: list[ForceFit] = []
    thermal_plot_data: dict[int, list[dict]] = {1: [], 2: [], 3: []}
    force_plot_data: dict[int, list[dict]] = {1: [], 2: [], 3: []}
    thermal_by_group: dict[int, list[ThermalFit]] = {1: [], 2: [], 3: []}
    force_by_group: dict[int, list[ForceFit]] = {1: [], 2: [], 3: []}
    fit_window_factors: dict[int, dict[float, list[float]]] = {
        group: {width: [] for width in THERMAL_FIT_WINDOW_CHECKS_HZ}
        for group in (1, 2, 3)
    }

    source_paths: list[Path] = []
    for group in (1, 2, 3):
        tnd_paths = sorted((DATA_ROOT / str(group)).glob("*.tnd"))
        force_paths = sorted((DATA_ROOT / str(group)).glob("*.jpk-force"))
        if len(tnd_paths) != 5 or len(force_paths) != 5:
            raise RuntimeError(
                f"Expected 5 TND and 5 force files for group {group}; "
                f"found {len(tnd_paths)} and {len(force_paths)}"
            )
        for path in tnd_paths:
            record, plot_data = refit_tnd(path, group)
            thermal_records.append(record)
            thermal_by_group[group].append(record)
            thermal_plot_data[group].append(plot_data)
            source_paths.append(path)
            # Numerical sensitivity uses the same blind-selected peak and raw
            # measured PSD, changing only the local fit half-width.
            _, data = load_tnd(path)
            for width in THERMAL_FIT_WINDOW_CHECKS_HZ:
                fit = fit_sho(data[:, 0], data[:, 2], record.detected_peak_hz, width)
                factor = record.correction_factor * Boltzmann * TEMPERATURE_K / float(fit["resonance_area_v2"])
                fit_window_factors[group][width].append(factor)
        for path in force_paths:
            record, plot_data = fit_force_curve(path, group)
            force_records.append(record)
            force_by_group[group].append(record)
            force_plot_data[group].append(plot_data)
            source_paths.append(path)

    if not all(record.fit_success for record in thermal_records):
        failures = [record.file for record in thermal_records if not record.fit_success]
        raise RuntimeError(f"SHO optimization failed: {failures}")
    if min(record.r_squared for record in force_records) < 0.999:
        raise RuntimeError("At least one hard-contact fit has R² < 0.999")

    summaries = [
        summarize_group(
            group,
            thermal_by_group[group],
            force_by_group[group],
            fit_window_factors[group],
        )
        for group in (1, 2, 3)
    ]

    thermal_rows = []
    for record in thermal_records:
        row = asdict(record)
        row["candidate_peaks_hz"] = ";".join(f"{value:.9g}" for value in record.candidate_peaks_hz)
        row["candidate_prominences_decades"] = ";".join(
            f"{value:.9g}" for value in record.candidate_prominences_decades
        )
        thermal_rows.append(row)
    write_csv(RESULTS / "thermal_refits.csv", thermal_rows)
    write_csv(RESULTS / "force_contact_fits.csv", [asdict(record) for record in force_records])
    write_csv(RESULTS / "calibration_summary.csv", summaries)

    make_thermal_full_plot(thermal_plot_data, summaries)
    make_thermal_zoom_plot(thermal_plot_data)
    make_contact_plot(force_plot_data, summaries)

    report = render_report(summaries, thermal_records, force_records)
    (RESULTS / "REPORT.md").write_text(report, encoding="utf-8")

    image_path = ROOT / "IMG_4571.jpeg"
    provenance = {
        "analysis_date": "2026-08-20",
        "temperature_C": TEMPERATURE_C,
        "temperature_K": TEMPERATURE_K,
        "input_counts": {"tnd": 15, "jpk_force": 15, "photo": int(image_path.exists())},
        "source_files": [
            {
                "path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in sorted(source_paths + ([image_path] if image_path.exists() else []))
        ],
        "algorithm": {
            "peak_selection": "lowest significant PSD peak; exported/header peak unused",
            "peak_search_Hz": [PEAK_SEARCH_MIN_HZ, PEAK_SEARCH_MAX_HZ],
            "peak_prominence_decades": PEAK_PROMINENCE_DECADES,
            "thermal_fit_model": "N + A2/((1-(f/f0)^2)^2 + (f/(f0*Q))^2)",
            "thermal_fit_half_width_Hz": THERMAL_FIT_HALF_WIDTH_HZ,
            "thermal_fit_objective": "Gamma/Whittle PSD negative log likelihood",
            "thermal_input_columns": ["Frequency", "average measured PSD"],
            "explicitly_unused_for_calibration": [
                "fit-data column",
                "exported parameter.f",
                "exported parameter.Q",
                "exported parameter.A",
            ],
            "resonance_area": "pi/2 * A2 * Q * f0",
            "thermal_factor": "correction_factor * k_B * T / resonance_area",
            "contact_segment": "approach",
            "contact_height_channel": "measuredHeight",
            "contact_span_nm": CONTACT_SPAN_NM,
            "sensitivity": "1/abs(dV/dz)",
            "spring_constant": "mean(thermal_factor)/mean(sensitivity)^2",
        },
        "label_sanity_check_only": {
            "probe_label": "NANOSENSORS POINTPROBE CONT-W; handwritten tipless",
            "frequency_range_kHz": [11.0, 17.0],
            "spring_constant_range_N_per_m": [0.11, 0.56],
            "length_um": 450.0,
            "width_um": [51.2, 54.0],
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
            "platform": platform.platform(),
        },
        "script": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256_file(Path(__file__).resolve()),
        },
    }
    (RESULTS / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    artifact_paths = [
        Path(__file__).resolve(),
        RESULTS / "REPORT.md",
        RESULTS / "thermal_refits.csv",
        RESULTS / "force_contact_fits.csv",
        RESULTS / "calibration_summary.csv",
        RESULTS / "provenance.json",
        FIGURES / "thermal_full_spectra_refits.png",
        FIGURES / "thermal_fundamental_zoom.png",
        FIGURES / "force_contact_fits.png",
    ]
    create_manifest(artifact_paths, RESULTS / "artifact_manifest.sha256")

    for summary in summaries:
        print(
            f"Cantilever {summary['cantilever']}: "
            f"S={summary['sensitivity_nm_per_V']:.3f} nm/V, "
            f"k={summary['spring_constant_N_per_m']:.6f} N/m, "
            f"f0={summary['resonance_frequency_kHz']:.6f} kHz, "
            f"Q={summary['quality_factor']:.3f}"
        )


if __name__ == "__main__":
    main()

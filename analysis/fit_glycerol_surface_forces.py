#!/usr/bin/env python3
"""Raw-only calibration and sphere-plane surface-force fitting for 20-08-26.

The JPK deflection-distance and force conversions are deliberately not used.
Vertical deflection is decoded in volts, InvOLS is re-estimated from hard
contact, and force is rebuilt with the independently calibrated cantilever-1
spring constant.  A robust linear far-field photodiode drift is fitted and
subtracted per approach curve, with per-pixel 16x16 slope maps retained for
spatial QC.  Equilibrium force curves are obtained by extrapolating the 1/2/4
um/s approach data to zero speed before fitting a same-potential silica-sphere/
silica-plane Poisson-Boltzmann/Derjaguin model.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
import zipfile
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.constants import Boltzmann, elementary_charge, epsilon_0
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT / "20-08-26"
DEFAULT_RESULTS = ROOT / "analysis" / "surface_force_results"

TEMPERATURE_C = 25.6
TEMPERATURE_K = TEMPERATURE_C + 273.15
SPRING_CONSTANT_N_PER_M = 0.2969899086921243
SPRING_CONSTANT_REPEAT_SD_N_PER_M = 0.0092768833
PROBE_RADIUS_M = 4.546848945303745e-6
PROBE_RADIUS_RANGE_M = (3.980307e-6, 4.546849e-6)
HAMAKER_J = 2.4e-21

# Rounded from the 25 C, 0.57 MHz glycerol-water measurements tabulated at
# 0.00/9.88/20.33/30.19/39.67 wt% as
# 78.48/75.98/73.86/71.44/68.93.  These are also the values used by the
# reference repository.  Concentration is glycerol mass fraction.
EPSILON_R = {0: 78.5, 10: 76.0, 20: 73.9, 30: 71.4, 40: 68.9}
SPEED_LABELS_UM_PER_S = np.array([1.0, 2.0, 4.0], dtype=np.float64)
BIN_CENTERS_NM = np.arange(5.0, 905.0, 5.0, dtype=np.float64)
BIN_HALF_WIDTH_NM = 2.5

CONTACT_WINDOW_POINTS = 12
CONTACT_SEARCH_POINTS = 150
CONTACT_R2_WINDOW_MIN = 0.997
CONTACT_R2_COMBINED_MIN = 0.995
CONTACT_SENSITIVITY_RANGE_M_PER_V = (35e-9, 120e-9)
FAR_FIELD_FRACTION = 0.20
FAR_FIELD_MIN_POINTS = 80

_GL_NODES, _GL_WEIGHTS = leggauss(32)
_GL_T = (_GL_NODES + 1.0) / 2.0
_GL_W = _GL_WEIGHTS / 2.0


@dataclass
class ContactFit:
    sensitivity_m_per_V: float
    start: int
    stop: int
    slope_V_per_m: float
    intercept_V: float
    r2: float


@dataclass
class FarFieldDriftFit:
    slope_V_per_m: float
    intercept_V: float
    r2: float
    n_points: int
    n_inliers: int
    height_min_m: float
    height_max_m: float
    residual_mad_V: float


@dataclass
class RawCurve:
    point_index: int | None
    measured_height_m: np.ndarray
    deflection_V: np.ndarray
    duration_s: float
    contact_fit: ContactFit | None = None
    far_field_fit: FarFieldDriftFit | None = None


@dataclass
class SourceData:
    path: Path
    concentration_wt_percent: int
    source_type: str
    timestamp: datetime
    sha256: str
    stored_sensitivity_m_per_V: float
    stored_spring_constant_N_per_m: float
    map_grid_i: int | None
    map_grid_j: int | None
    map_back_and_forth: bool
    map_ulength_m: float
    map_vlength_m: float
    curves: list[RawCurve]
    skipped_curves: int
    sensitivity_anchor_m_per_V: float = float("nan")
    sensitivity_anchor_mad_m_per_V: float = float("nan")
    sensitivity_valid_curves: int = 0
    sensitivity_used_m_per_V: float = float("nan")
    sensitivity_method: str = ""


@dataclass
class PreparedSource:
    source: SourceData
    speed_um_per_s: float
    speed_label_um_per_s: float
    force_pN: np.ndarray
    force_spread_pN: np.ndarray
    counts: np.ndarray
    valid_contact_curves: int
    endpoint_contact_curves: int
    far_offset_pN: float


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _ends_with_unescaped_backslash(line: str) -> bool:
    count = 0
    for character in reversed(line):
        if character != "\\":
            break
        count += 1
    return count % 2 == 1


def _unescape_java(value: str) -> str:
    output: list[str] = []
    index = 0
    translations = {"t": "\t", "n": "\n", "r": "\r", "f": "\f"}
    while index < len(value):
        character = value[index]
        if character != "\\" or index + 1 >= len(value):
            output.append(character)
            index += 1
            continue
        escaped = value[index + 1]
        if escaped == "u" and index + 5 < len(value):
            token = value[index + 2 : index + 6]
            try:
                output.append(chr(int(token, 16)))
                index += 6
                continue
            except ValueError:
                pass
        output.append(translations.get(escaped, escaped))
        index += 2
    return "".join(output)


def _split_property(line: str) -> tuple[str, str]:
    escaped = False
    for index, character in enumerate(line):
        if escaped:
            escaped = False
            continue
        if character == "\\":
            escaped = True
            continue
        if character in "=:":
            return line[:index], line[index + 1 :]
        if character.isspace():
            remainder = line[index:].lstrip()
            if remainder.startswith(("=", ":")):
                remainder = remainder[1:].lstrip()
            return line[:index], remainder
    return line, ""


def parse_properties(raw: bytes) -> dict[str, str]:
    physical_lines = raw.decode("iso-8859-1").splitlines()
    logical_lines: list[str] = []
    pending = ""
    for physical in physical_lines:
        line = pending + (physical.lstrip() if pending else physical)
        if _ends_with_unescaped_backslash(line):
            pending = line[:-1]
            continue
        logical_lines.append(line)
        pending = ""
    if pending:
        logical_lines.append(pending)
    properties: dict[str, str] = {}
    for line in logical_lines:
        stripped = line.lstrip()
        if not stripped or stripped.startswith(("#", "!")):
            continue
        key, value = _split_property(stripped)
        properties[_unescape_java(key.strip())] = _unescape_java(value.strip())
    return properties


def _required(properties: Mapping[str, str], key: str) -> str:
    try:
        return properties[key]
    except KeyError as exc:
        raise ValueError(f"missing JPK property: {key}") from exc


def _raw_dtype(encoder_type: str) -> np.dtype:
    mapping = {
        "signedinteger": np.dtype(">i4"),
        "unsignedinteger": np.dtype(">u4"),
        "signedshort": np.dtype(">i2"),
        "unsignedshort": np.dtype(">u2"),
        "float-data": np.dtype(">f4"),
        "double-data": np.dtype(">f8"),
    }
    try:
        return mapping[encoder_type.strip().lower()]
    except KeyError as exc:
        raise ValueError(f"unsupported JPK encoder: {encoder_type}") from exc


def decode_base_channel(
    archive: zipfile.ZipFile,
    shared: Mapping[str, str],
    segment: Mapping[str, str],
    segment_prefix: str,
    channel: str,
) -> tuple[np.ndarray, int]:
    lcd_index = int(_required(segment, f"channel.{channel}.lcd-info.*"))
    lcd = f"lcd-info.{lcd_index}"
    member = f"{segment_prefix}/{_required(segment, f'channel.{channel}.data.file.name')}"
    payload = archive.read(member)
    dtype = _raw_dtype(_required(shared, f"{lcd}.encoder.type"))
    if len(payload) % dtype.itemsize:
        raise ValueError(f"{member} has invalid byte length")
    raw = np.frombuffer(payload, dtype=dtype).astype(np.float64)
    offset = float(_required(shared, f"{lcd}.encoder.scaling.offset"))
    multiplier = float(_required(shared, f"{lcd}.encoder.scaling.multiplier"))
    values = offset + multiplier * raw
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{member} decoded non-finite values")
    return values, lcd_index


def apply_linear_conversion(
    values: np.ndarray,
    shared: Mapping[str, str],
    lcd_index: int,
    slot: str,
) -> np.ndarray:
    lcd = f"lcd-info.{lcd_index}"
    prefix = f"{lcd}.conversion-set.conversion.{slot}"
    if _required(shared, f"{prefix}.defined").lower() != "true":
        raise ValueError(f"conversion {slot} is undefined for LCD {lcd_index}")
    if _required(shared, f"{prefix}.scaling.type") != "linear":
        raise ValueError(f"conversion {slot} is not linear")
    offset = float(_required(shared, f"{prefix}.scaling.offset"))
    multiplier = float(_required(shared, f"{prefix}.scaling.multiplier"))
    converted = offset + multiplier * values
    if not np.all(np.isfinite(converted)):
        raise ValueError(f"conversion {slot} produced non-finite values")
    return converted


def _segment_style(segment: Mapping[str, str], shared: Mapping[str, str]) -> str:
    direct = segment.get("force-segment-header.settings.style")
    if direct:
        return direct.lower()
    reference = _required(segment, "force-segment-header.force-segment-header-info.*")
    return _required(
        shared, f"force-segment-header-info.{reference}.settings.style"
    ).lower()


def _source_timestamp(path: Path) -> datetime:
    match = re.search(
        r"(\d{4}\.\d{2}\.\d{2}-\d{2}\.\d{2}\.\d{2}\.\d{3})", path.name
    )
    if not match:
        raise ValueError(f"timestamp is absent from {path.name}")
    return datetime.strptime(match.group(1), "%Y.%m.%d-%H.%M.%S.%f")


def _segment_prefixes(
    names: Iterable[str], shared: Mapping[str, str], archive: zipfile.ZipFile, source_type: str
) -> list[tuple[int | None, str]]:
    if source_type == "force":
        candidates: list[tuple[int | None, int, str]] = []
        pattern = re.compile(r"^segments/(\d+)/segment-header\.properties$")
        for name in names:
            match = pattern.match(name)
            if match:
                candidates.append((None, int(match.group(1)), name.rsplit("/", 1)[0]))
    else:
        candidates = []
        pattern = re.compile(
            r"^index/(\d+)/segments/(\d+)/segment-header\.properties$"
        )
        for name in names:
            match = pattern.match(name)
            if match:
                candidates.append(
                    (int(match.group(1)), int(match.group(2)), name.rsplit("/", 1)[0])
                )
    output: list[tuple[int | None, str]] = []
    for point_index, segment_index, prefix in sorted(
        candidates, key=lambda item: (-1 if item[0] is None else item[0], item[1])
    ):
        segment = parse_properties(archive.read(f"{prefix}/segment-header.properties"))
        if _segment_style(segment, shared) == "extend":
            output.append((point_index, prefix))
    return output


def load_source(path: Path, concentration: int) -> SourceData:
    source_type = "map" if path.suffix == ".jpk-force-map" else "force"
    curves: list[RawCurve] = []
    skipped = 0
    with zipfile.ZipFile(path, "r") as archive:
        corrupt = archive.testzip()
        if corrupt is not None:
            raise ValueError(f"CRC failure in {path}: {corrupt}")
        header = parse_properties(archive.read("header.properties"))
        shared = parse_properties(archive.read("shared-data/header.properties"))
        if source_type == "map":
            map_grid_i = int(
                _required(header, "force-scan-map.position-pattern.grid.ilength")
            )
            map_grid_j = int(
                _required(header, "force-scan-map.position-pattern.grid.jlength")
            )
            map_back_and_forth = (
                _required(header, "force-scan-map.position-pattern.back-and-forth")
                .strip()
                .lower()
                == "true"
            )
            map_ulength = float(
                _required(header, "force-scan-map.position-pattern.grid.ulength")
            )
            map_vlength = float(
                _required(header, "force-scan-map.position-pattern.grid.vlength")
            )
            if (
                map_grid_i <= 0
                or map_grid_j <= 0
                or not np.isfinite(map_ulength)
                or not np.isfinite(map_vlength)
                or map_ulength <= 0.0
                or map_vlength <= 0.0
            ):
                raise ValueError(f"invalid map grid metadata in {path}")
        else:
            map_grid_i = None
            map_grid_j = None
            map_back_and_forth = False
            map_ulength = float("nan")
            map_vlength = float("nan")
        prefixes = _segment_prefixes(archive.namelist(), shared, archive, source_type)
        stored_sensitivity = float("nan")
        stored_spring = float("nan")
        for point_index, prefix in prefixes:
            try:
                segment = parse_properties(
                    archive.read(f"{prefix}/segment-header.properties")
                )
                deflection, v_lcd = decode_base_channel(
                    archive, shared, segment, prefix, "vDeflection"
                )
                measured_base, h_lcd = decode_base_channel(
                    archive, shared, segment, prefix, "measuredHeight"
                )
                measured_height = apply_linear_conversion(
                    measured_base, shared, h_lcd, "nominal"
                )
                duration = float(_required(segment, "force-segment-header.duration"))
                if (
                    deflection.shape != measured_height.shape
                    or deflection.size < 200
                    or not np.isfinite(duration)
                    or duration <= 0.0
                ):
                    raise ValueError("short or inconsistent approach segment")
                if measured_height[0] < measured_height[-1]:
                    measured_height = measured_height[::-1].copy()
                    deflection = deflection[::-1].copy()
                curves.append(
                    RawCurve(
                        point_index=point_index,
                        measured_height_m=measured_height,
                        deflection_V=deflection,
                        duration_s=duration,
                    )
                )
                if not np.isfinite(stored_sensitivity):
                    stored_sensitivity = float(
                        _required(
                            shared,
                            f"lcd-info.{v_lcd}.conversion-set.conversion.distance.scaling.multiplier",
                        )
                    )
                    stored_spring = float(
                        _required(
                            shared,
                            f"lcd-info.{v_lcd}.conversion-set.conversion.force.scaling.multiplier",
                        )
                    )
            except (KeyError, ValueError, zipfile.BadZipFile):
                skipped += 1
    if not curves:
        raise ValueError(f"no usable approach curves in {path}")
    return SourceData(
        path=path,
        concentration_wt_percent=concentration,
        source_type=source_type,
        timestamp=_source_timestamp(path),
        sha256=sha256_file(path),
        stored_sensitivity_m_per_V=stored_sensitivity,
        stored_spring_constant_N_per_m=stored_spring,
        map_grid_i=map_grid_i,
        map_grid_j=map_grid_j,
        map_back_and_forth=map_back_and_forth,
        map_ulength_m=map_ulength,
        map_vlength_m=map_vlength,
        curves=curves,
        skipped_curves=skipped,
    )


def robust_mad(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if not values.size:
        return float("nan")
    median = float(np.median(values))
    return float(1.4826 * np.median(np.abs(values - median)))


def robust_line(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, int]:
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    y = np.asarray(y, dtype=np.float64).reshape(-1)
    if x.size != y.size or x.size < 8 or not np.all(np.isfinite(x + y)):
        raise ValueError("line fit requires at least eight finite pairs")
    keep = np.ones(x.size, dtype=bool)
    for _ in range(4):
        selected_x = x[keep]
        selected_y = y[keep]
        center = float(np.mean(selected_x))
        scale = float(np.std(selected_x))
        if scale <= np.finfo(np.float64).tiny:
            raise ValueError("singular x range")
        design = np.column_stack(
            [(selected_x - center) / scale, np.ones(selected_x.size)]
        )
        coefficients, _, rank, _ = np.linalg.lstsq(design, selected_y, rcond=None)
        if rank != 2:
            raise ValueError("rank-deficient line fit")
        slope = float(coefficients[0] / scale)
        intercept = float(coefficients[1] - slope * center)
        residual = y - (slope * x + intercept)
        sigma = robust_mad(residual[keep])
        if not np.isfinite(sigma) or sigma <= np.finfo(np.float64).tiny:
            break
        new_keep = np.abs(residual - np.median(residual[keep])) <= 3.5 * sigma
        if np.count_nonzero(new_keep) < 8 or np.array_equal(new_keep, keep):
            break
        keep = new_keep
    prediction = slope * x[keep] + intercept
    ss_res = float(np.sum((y[keep] - prediction) ** 2))
    ss_tot = float(np.sum((y[keep] - np.mean(y[keep])) ** 2))
    r2 = float("nan") if ss_tot <= 0.0 else 1.0 - ss_res / ss_tot
    return slope, intercept, r2, int(np.count_nonzero(keep))


def fit_far_field_drift(
    height_m: np.ndarray, deflection_V: np.ndarray
) -> FarFieldDriftFit:
    """Fit the per-approach linear photodiode drift in the initial far field.

    JPK approach arrays are normalized above so that acquisition starts at the
    largest scanner height, approximately 0.8--1.0 um from contact.  The same
    first-20-percent definition was already used by this analysis; this helper
    makes the fitted line and its diagnostics explicit without changing the
    subtraction operator.
    """

    height = np.asarray(height_m, dtype=np.float64).reshape(-1)
    voltage = np.asarray(deflection_V, dtype=np.float64).reshape(-1)
    if height.shape != voltage.shape or height.size < FAR_FIELD_MIN_POINTS:
        raise ValueError("far-field drift fit needs equal, sufficiently long arrays")
    points = min(
        height.size,
        max(FAR_FIELD_MIN_POINTS, int(math.ceil(FAR_FIELD_FRACTION * height.size))),
    )
    far_height = height[:points]
    far_voltage = voltage[:points]
    slope, intercept, r2, n_inliers = robust_line(far_height, far_voltage)
    residual = far_voltage - (slope * far_height + intercept)
    return FarFieldDriftFit(
        slope_V_per_m=float(slope),
        intercept_V=float(intercept),
        r2=float(r2),
        n_points=int(points),
        n_inliers=int(n_inliers),
        height_min_m=float(np.min(far_height)),
        height_max_m=float(np.max(far_height)),
        residual_mad_V=robust_mad(residual),
    )


def baseline_voltage(
    height_m: np.ndarray,
    deflection_V: np.ndarray,
    fit: FarFieldDriftFit | None = None,
) -> np.ndarray:
    selected = fit if fit is not None else fit_far_field_drift(height_m, deflection_V)
    return selected.slope_V_per_m * height_m + selected.intercept_V


def find_contact_fit(
    height_m: np.ndarray,
    deflection_V: np.ndarray,
    far_field_fit: FarFieldDriftFit | None = None,
) -> ContactFit | None:
    try:
        corrected = deflection_V - baseline_voltage(
            height_m, deflection_V, far_field_fit
        )
    except (ValueError, np.linalg.LinAlgError):
        return None
    n = height_m.size
    width = CONTACT_WINDOW_POINTS
    first = max(0, n - CONTACT_SEARCH_POINTS)
    candidates: list[int] = []
    for start in range(first, n - width + 1):
        try:
            slope, _, r2, _ = robust_line(
                height_m[start : start + width], corrected[start : start + width]
            )
        except (ValueError, np.linalg.LinAlgError):
            continue
        if slope >= 0.0 or not np.isfinite(r2) or r2 < CONTACT_R2_WINDOW_MIN:
            continue
        sensitivity = 1.0 / abs(slope)
        if CONTACT_SENSITIVITY_RANGE_M_PER_V[0] <= sensitivity <= CONTACT_SENSITIVITY_RANGE_M_PER_V[1]:
            candidates.append(start)
    if not candidates:
        return None
    runs: list[list[int]] = []
    current = [candidates[0]]
    for start in candidates[1:]:
        if start == current[-1] + 1:
            current.append(start)
        else:
            runs.append(current)
            current = [start]
    runs.append(current)
    # Hard contact must occur at the terminal end of approach.  Allow a small
    # terminal margin for setpoint/feedback glitches.
    terminal_runs = [run for run in runs if run[-1] + width >= n - 15]
    if not terminal_runs:
        return None
    for run in sorted(terminal_runs, key=len, reverse=True):
        start = run[0]
        stop = run[-1] + width
        try:
            slope, intercept, r2, _ = robust_line(
                height_m[start:stop], corrected[start:stop]
            )
        except (ValueError, np.linalg.LinAlgError):
            continue
        sensitivity = 1.0 / abs(slope) if slope < 0.0 else float("inf")
        if (
            r2 >= CONTACT_R2_COMBINED_MIN
            and CONTACT_SENSITIVITY_RANGE_M_PER_V[0]
            <= sensitivity
            <= CONTACT_SENSITIVITY_RANGE_M_PER_V[1]
        ):
            return ContactFit(
                sensitivity_m_per_V=float(sensitivity),
                start=start,
                stop=stop,
                slope_V_per_m=float(slope),
                intercept_V=float(intercept),
                r2=float(r2),
            )
    return None


def calibrate_sensitivity(sources: list[SourceData]) -> None:
    for source in sources:
        values: list[float] = []
        for curve in source.curves:
            try:
                curve.far_field_fit = fit_far_field_drift(
                    curve.measured_height_m, curve.deflection_V
                )
            except (ValueError, np.linalg.LinAlgError):
                curve.far_field_fit = None
            curve.contact_fit = find_contact_fit(
                curve.measured_height_m,
                curve.deflection_V,
                curve.far_field_fit,
            )
            if curve.contact_fit is not None:
                values.append(curve.contact_fit.sensitivity_m_per_V)
        array = np.asarray(values, dtype=np.float64)
        if array.size:
            median = float(np.median(array))
            mad = robust_mad(array)
            tolerance = max(3.5 * mad if np.isfinite(mad) else 0.0, 5e-9)
            retained = array[np.abs(array - median) <= tolerance]
        else:
            retained = array
        required = 1 if source.source_type == "force" else max(
            10, int(math.ceil(0.03 * len(source.curves)))
        )
        source.sensitivity_valid_curves = int(retained.size)
        if retained.size >= required:
            source.sensitivity_anchor_m_per_V = float(np.median(retained))
            source.sensitivity_anchor_mad_m_per_V = robust_mad(retained)

    # InvOLS is an optical-lever calibration, not a per-curve fit parameter.
    # Use every accepted hard-contact interval to establish one consensus per
    # liquid condition.  Equal source weighting prevents a 256-point map from
    # overwhelming the independent force curves.  A lone short-contact anchor
    # more than 15 nm/V from the concentration median is not allowed to set the
    # force scale (the rejected anchor remains visible in the output table).
    for concentration in sorted(EPSILON_R):
        subset = [
            source
            for source in sources
            if source.concentration_wt_percent == concentration
            and np.isfinite(source.sensitivity_anchor_m_per_V)
        ]
        if not subset:
            raise RuntimeError(
                f"no raw hard-contact sensitivity anchor at {concentration} wt%"
            )
        values = np.array(
            [source.sensitivity_anchor_m_per_V for source in subset], dtype=float
        )
        center = float(np.median(values))
        retained = values[np.abs(values - center) <= 15e-9]
        if retained.size < 2:
            raise RuntimeError(
                f"fewer than two consistent sensitivity anchors at {concentration} wt%"
            )
        consensus = float(np.median(retained))
        for source in sources:
            if source.concentration_wt_percent == concentration:
                source.sensitivity_used_m_per_V = consensus
                source.sensitivity_method = "concentration_consensus_raw_contacts"


def _speed_label(speed_um_per_s: float) -> float:
    return float(SPEED_LABELS_UM_PER_S[np.argmin(np.abs(SPEED_LABELS_UM_PER_S - speed_um_per_s))])


def prepare_source(source: SourceData) -> PreparedSource:
    binned_curves: list[np.ndarray] = []
    speeds: list[float] = []
    valid_contact = 0
    endpoint_contact = 0
    sensitivity = source.sensitivity_used_m_per_V
    for curve in source.curves:
        height = curve.measured_height_m
        voltage = curve.deflection_V
        if curve.far_field_fit is None:
            try:
                curve.far_field_fit = fit_far_field_drift(height, voltage)
            except (ValueError, np.linalg.LinAlgError):
                continue
        try:
            corrected = voltage - baseline_voltage(
                height, voltage, curve.far_field_fit
            )
        except (ValueError, np.linalg.LinAlgError):
            continue
        delta_m = sensitivity * corrected
        if curve.contact_fit is not None:
            start = curve.contact_fit.start
            stop = curve.contact_fit.stop
            try:
                slope, intercept, _, _ = robust_line(
                    height[start:stop], corrected[start:stop]
                )
                contact_height = -intercept / slope
            except (ValueError, np.linalg.LinAlgError, ZeroDivisionError):
                continue
            precontact_index = np.arange(height.size) < start
            valid_contact += 1
        else:
            # The last acquired point is assumed to have reached the mechanical
            # contact plane.  h + delta is constant in rigid contact; its upper
            # terminal quantile is more stable than a single final sample.
            terminal = height[-20:] + delta_m[-20:]
            contact_height = float(np.quantile(terminal, 0.95))
            precontact_index = np.ones(height.size, dtype=bool)
            endpoint_contact += 1
        distance_nm = (height + delta_m - contact_height) * 1e9
        force_pN = SPRING_CONSTANT_N_PER_M * delta_m * 1e12
        usable = (
            precontact_index
            & np.isfinite(distance_nm)
            & np.isfinite(force_pN)
            & (distance_nm >= BIN_CENTERS_NM[0] - BIN_HALF_WIDTH_NM)
            & (distance_nm < BIN_CENTERS_NM[-1] + BIN_HALF_WIDTH_NM)
        )
        if np.count_nonzero(usable) < 30:
            continue
        row = np.full(BIN_CENTERS_NM.shape, np.nan, dtype=np.float64)
        indices = np.floor(
            (distance_nm[usable] - (BIN_CENTERS_NM[0] - BIN_HALF_WIDTH_NM)) / 5.0
        ).astype(int)
        values = force_pN[usable]
        for index in np.unique(indices):
            if 0 <= index < row.size:
                row[index] = float(np.median(values[indices == index]))
        binned_curves.append(row)
        sample_time = (
            np.arange(height.size, dtype=np.float64) + 0.5
        ) * curve.duration_s / height.size
        trim = max(5, int(math.ceil(0.10 * height.size)))
        slope, _, _, _ = robust_line(
            sample_time[trim:-trim], height[trim:-trim]
        )
        speeds.append(abs(slope) * 1e6)
    if not binned_curves:
        raise ValueError(f"no prepared curves for {source.path}")
    matrix = np.asarray(binned_curves, dtype=np.float64)
    force = np.full(BIN_CENTERS_NM.shape, np.nan, dtype=np.float64)
    counts = np.sum(np.isfinite(matrix), axis=0).astype(int)
    spread = np.full(force.shape, np.nan)
    for index in range(force.size):
        finite = matrix[:, index][np.isfinite(matrix[:, index])]
        if finite.size:
            force[index] = float(np.median(finite))
        if finite.size >= 2:
            spread[index] = robust_mad(finite)
    # Residual source offsets are referenced in a region much farther than the
    # expected EDL range; the raw photodiode baseline itself was fit at the
    # beginning of each approach (~0.8-1.0 um from contact).
    far = (BIN_CENTERS_NM >= 700.0) & (BIN_CENTERS_NM <= 880.0) & np.isfinite(force)
    far_offset = float(np.median(force[far])) if np.count_nonzero(far) >= 5 else 0.0
    force = force - far_offset
    speed = float(np.median(speeds))
    return PreparedSource(
        source=source,
        speed_um_per_s=speed,
        speed_label_um_per_s=_speed_label(speed),
        force_pN=force,
        force_spread_pN=spread,
        counts=counts,
        valid_contact_curves=valid_contact,
        endpoint_contact_curves=endpoint_contact,
        far_offset_pN=far_offset,
    )


def cheng_viscosity_mPa_s(glycerol_mass_fraction: float, temperature_C: float) -> float:
    concentration = float(glycerol_mass_fraction)
    temperature = float(temperature_C)
    if not 0.0 <= concentration <= 1.0 or not 0.0 <= temperature <= 100.0:
        raise ValueError("Cheng correlation requires mass fraction [0,1], T [0,100 C]")
    a = 0.705 - 0.0017 * temperature
    b = (4.9 + 0.036 * temperature) * a**2.5
    alpha = 1.0 - concentration
    if 0.0 < concentration < 1.0:
        alpha += (
            a
            * b
            * concentration
            * (1.0 - concentration)
            / (a * concentration + b * (1.0 - concentration))
        )
    mu_water = 1.790 * np.exp(
        ((-1230.0 - temperature) * temperature)
        / (36100.0 + 360.0 * temperature)
    )
    mu_glycerol = 12100.0 * np.exp(
        ((-1233.0 + temperature) * temperature)
        / (9900.0 + 70.0 * temperature)
    )
    return float(mu_water**alpha * mu_glycerol ** (1.0 - alpha))


def zero_speed_curve(
    prepared: list[PreparedSource],
    concentration: int,
    source_filter: str = "all",
) -> list[dict[str, float | int | str]]:
    selected = [
        item
        for item in prepared
        if item.source.concentration_wt_percent == concentration
        and (source_filter == "all" or item.source.source_type == source_filter)
    ]
    records: list[dict[str, float | int | str]] = []
    for bin_index, distance_nm in enumerate(BIN_CENTERS_NM):
        group_speed: list[float] = []
        group_force: list[float] = []
        group_sigma: list[float] = []
        group_sources: list[int] = []
        for label in SPEED_LABELS_UM_PER_S:
            members = [
                item
                for item in selected
                if item.speed_label_um_per_s == label
                and np.isfinite(item.force_pN[bin_index])
                and item.counts[bin_index] > 0
            ]
            if not members:
                continue
            values = np.array([item.force_pN[bin_index] for item in members])
            speeds = np.array([item.speed_um_per_s for item in members])
            spreads = np.array(
                [item.force_spread_pN[bin_index] for item in members], dtype=float
            )
            group_speed.append(float(np.median(speeds)))
            group_force.append(float(np.median(values)))
            between = robust_mad(values)
            within = float(np.nanmedian(spreads)) if np.any(np.isfinite(spreads)) else 0.0
            sigma = max(
                between if np.isfinite(between) else 0.0,
                within / math.sqrt(max(1, int(np.median([m.counts[bin_index] for m in members]))))
                if np.isfinite(within)
                else 0.0,
                3.0,
            )
            group_sigma.append(float(sigma))
            group_sources.append(len(members))
        if len(group_speed) < 2:
            continue
        x = np.asarray(group_speed, dtype=float)
        y = np.asarray(group_force, dtype=float)
        sigma = np.asarray(group_sigma, dtype=float)
        design = np.column_stack([np.ones(x.size), x])
        root_weight = 1.0 / sigma
        weighted_design = design * root_weight[:, None]
        weighted_y = y * root_weight
        coefficients, _, rank, singular = np.linalg.lstsq(
            weighted_design, weighted_y, rcond=None
        )
        if rank != 2:
            continue
        prediction = design @ coefficients
        residual = y - prediction
        information = weighted_design.T @ weighted_design
        covariance = np.linalg.inv(information)
        if x.size > 2:
            reduced = float(np.sum((residual / sigma) ** 2) / (x.size - 2))
            covariance *= max(reduced, 1.0)
        intercept_sigma = float(math.sqrt(max(0.0, covariance[0, 0])))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        speed_r2 = float("nan") if ss_tot <= 0.0 else 1.0 - float(np.sum(residual**2)) / ss_tot
        records.append(
            {
                "concentration_wt_percent": concentration,
                "source_filter": source_filter,
                "distance_nm": float(distance_nm),
                "equilibrium_force_pN": float(coefficients[0]),
                "equilibrium_sigma_pN": max(intercept_sigma, 3.0),
                "speed_slope_pN_per_um_s": float(coefficients[1]),
                "speed_regression_r2": speed_r2,
                "speed_groups": len(x),
                "source_count": int(sum(group_sources)),
                "speed_1_force_pN": float(y[np.argmin(abs(x - 1.0))]) if np.any(abs(x - 1.0) < 0.4) else float("nan"),
                "speed_2_force_pN": float(y[np.argmin(abs(x - 2.0))]) if np.any(abs(x - 2.0) < 0.4) else float("nan"),
                "speed_4_force_pN": float(y[np.argmin(abs(x - 4.0))]) if np.any(abs(x - 4.0) < 0.5) else float("nan"),
            }
        )
    return records


@lru_cache(maxsize=192)
def _pb_table(surface: float) -> tuple[np.ndarray, np.ndarray, float]:
    """Build a reusable nonlinear-PB pressure-integral table for one u_s."""

    grid = np.unique(
        np.concatenate(
            [
                np.geomspace(1e-4, 1.0, 160),
                np.linspace(1.0, 40.0, 420),
            ]
        )
    )
    target = grid / 2.0
    low = np.full(grid.shape, max(surface * 1e-14, 1e-15))
    high = np.full(grid.shape, surface * (1.0 - 1e-13))

    def half_gap(midplane: np.ndarray) -> np.ndarray:
        difference = surface - midplane
        t2 = _GL_T[None, :] ** 2
        increment = difference[:, None] * t2
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            denominator = np.sqrt(
                np.sinh(midplane[:, None] + 0.5 * increment)
                * np.sinh(0.5 * increment)
            )
            integrand = difference[:, None] * _GL_T[None, :] / denominator
        return integrand @ _GL_W

    for _ in range(55):
        middle = 0.5 * (low + high)
        value = half_gap(middle)
        # The half-gap decreases monotonically as the midplane potential rises.
        low = np.where(value > target, middle, low)
        high = np.where(value > target, high, middle)
    midplane = 0.5 * (low + high)
    pressure = 2.0 * np.sinh(0.5 * midplane) ** 2
    gamma = math.tanh(surface / 4.0)
    tail = 32.0 * gamma**2 * math.exp(-grid[-1])
    integral = np.empty_like(grid)
    integral[-1] = tail
    increments = np.diff(grid)
    for index in range(grid.size - 2, -1, -1):
        integral[index] = integral[index + 1] + 0.5 * (
            pressure[index] + pressure[index + 1]
        ) * increments[index]
    return grid, integral, gamma


def pb_dimensionless_G(H: np.ndarray, surface_potential_dimensionless: float) -> np.ndarray:
    """Return integral of dimensionless PB disjoining pressure.

    For equal constant-potential plates in a symmetric 1:1 electrolyte,
    ``G(H,u_s) = integral_H^infinity [cosh(u_m(s))-1] ds``.  ``u_m`` is
    obtained from the first integral of nonlinear PB.  The sphere-plane force
    is ``2*pi*R*eps*(kT/e)^2*kappa*G``.
    """

    query = np.asarray(H, dtype=np.float64)
    if not np.all(np.isfinite(query)) or np.any(query <= 0.0):
        raise ValueError("dimensionless separations must be finite and positive")
    surface = abs(float(surface_potential_dimensionless))
    if not np.isfinite(surface) or surface < 0.0:
        raise ValueError("surface potential must be finite")
    if surface == 0.0:
        return np.zeros_like(query)
    flat = query.reshape(-1)
    grid, integral, gamma = _pb_table(surface)
    result = np.empty_like(flat)
    # At very large H the midplane potential is so small that subtracting
    # cosh(u_m)-1 loses relative precision in float64.  The nonlinear
    # superposition asymptote is exponentially accurate there.
    inside = flat <= 25.0
    result[inside] = np.interp(flat[inside], grid, integral)
    result[~inside] = 32.0 * gamma**2 * np.exp(-flat[~inside])
    if not np.all(np.isfinite(result)) or np.any(result < 0.0):
        raise FloatingPointError("nonlinear PB evaluator produced invalid values")
    return result.reshape(query.shape)


def edl_force_pN(
    distance_nm: np.ndarray,
    lambda_nm: float,
    zeta_mV: float,
    eps_r: float,
    model: str,
) -> np.ndarray:
    distance_m = np.asarray(distance_nm, dtype=float) * 1e-9
    kappa = 1.0 / (float(lambda_nm) * 1e-9)
    zeta_V = abs(float(zeta_mV)) * 1e-3
    absolute_permittivity = eps_r * epsilon_0
    H = kappa * distance_m
    if model == "nonlinear_pb_derjaguin":
        u_surface = elementary_charge * zeta_V / (Boltzmann * TEMPERATURE_K)
        G = pb_dimensionless_G(H, u_surface)
        force_N = (
            2.0
            * np.pi
            * PROBE_RADIUS_M
            * absolute_permittivity
            * (Boltzmann * TEMPERATURE_K / elementary_charge) ** 2
            * kappa
            * G
        )
    elif model == "linear_hhf_equal_potential":
        reciprocal = np.where(H > 40.0, np.exp(-H), 1.0 / (np.exp(H) + 1.0))
        force_N = (
            4.0
            * np.pi
            * PROBE_RADIUS_M
            * absolute_permittivity
            * kappa
            * zeta_V**2
            * reciprocal
        )
    else:
        raise ValueError(f"unknown EDL model: {model}")
    return force_N * 1e12


def total_equilibrium_force_pN(
    distance_nm: np.ndarray,
    lambda_nm: float,
    zeta_mV: float,
    baseline_pN: float,
    eps_r: float,
    model: str,
) -> np.ndarray:
    distance_m = np.asarray(distance_nm, dtype=float) * 1e-9
    vdw_pN = -HAMAKER_J * PROBE_RADIUS_M / (6.0 * distance_m**2) * 1e12
    return (
        edl_force_pN(distance_nm, lambda_nm, zeta_mV, eps_r, model)
        + vdw_pN
        + baseline_pN
    )


def fit_equilibrium(
    records: list[dict[str, float | int | str]],
    concentration: int,
    model: str,
    dmin_nm: float,
    dmax_nm: float,
    variant: str,
) -> dict[str, float | int | str | bool]:
    selected = [
        record
        for record in records
        if dmin_nm <= float(record["distance_nm"]) <= dmax_nm
        and int(record["speed_groups"]) >= 3
    ]
    if len(selected) < 12:
        return {
            "concentration_wt_percent": concentration,
            "variant": variant,
            "model": model,
            "fit_valid": False,
            "invalid_reasons": "fewer_than_12_three_speed_bins",
            "n_points": len(selected),
            "dmin_nm": dmin_nm,
            "dmax_nm": dmax_nm,
        }
    distance = np.array([float(item["distance_nm"]) for item in selected])
    force = np.array([float(item["equilibrium_force_pN"]) for item in selected])
    sigma = np.array([float(item["equilibrium_sigma_pN"]) for item in selected])
    far = [
        float(item["equilibrium_force_pN"])
        for item in records
        if 500.0 <= float(item["distance_nm"]) <= 850.0
    ]
    noise_floor = max(3.0, robust_mad(np.asarray(far, dtype=float))) if far else 3.0
    sigma = np.maximum(sigma, noise_floor)
    eps_r = EPSILON_R[concentration]

    def residual(parameters: np.ndarray) -> np.ndarray:
        lambda_nm = math.exp(float(parameters[0]))
        zeta_mV = math.exp(float(parameters[1]))
        baseline = float(parameters[2])
        prediction = total_equilibrium_force_pN(
            distance, lambda_nm, zeta_mV, baseline, eps_r, model
        )
        return (prediction - force) / sigma

    lower = np.array([math.log(1.0), math.log(0.1), -500.0])
    upper = np.array([math.log(1000.0), math.log(250.0), 500.0])
    best = None
    for initial_lambda, initial_zeta in ((20.0, 35.0), (100.0, 100.0)):
        result = least_squares(
            residual,
            x0=np.array(
                [math.log(initial_lambda), math.log(initial_zeta), 0.0]
            ),
            bounds=(lower, upper),
            method="trf",
            loss="soft_l1",
            f_scale=1.0,
            xtol=1e-8,
            ftol=1e-8,
            gtol=1e-8,
            max_nfev=250,
        )
        objective = float(np.sum(residual(result.x) ** 2))
        if best is None or objective < best[0]:
            best = (objective, result)
    assert best is not None
    objective, result = best
    lambda_nm = math.exp(float(result.x[0]))
    zeta_mV = math.exp(float(result.x[1]))
    baseline = float(result.x[2])
    prediction = total_equilibrium_force_pN(
        distance, lambda_nm, zeta_mV, baseline, eps_r, model
    )
    raw_residual = force - prediction
    ss_res = float(np.sum(raw_residual**2))
    ss_tot = float(np.sum((force - np.mean(force)) ** 2))
    r2 = float("nan") if ss_tot <= 0.0 else 1.0 - ss_res / ss_tot
    jacobian = np.asarray(result.jac, dtype=float)
    column_norm = np.linalg.norm(jacobian, axis=0)
    rank = 0
    condition = float("inf")
    se = np.full(3, np.nan)
    if np.all(np.isfinite(column_norm)) and np.all(column_norm > 0.0):
        normalized = jacobian / column_norm
        singular = np.linalg.svd(normalized, compute_uv=False)
        tolerance = max(normalized.shape) * np.finfo(float).eps * singular[0]
        rank = int(np.count_nonzero(singular > tolerance))
        if singular[-1] > 0.0:
            condition = float(singular[0] / singular[-1])
        if rank == 3 and len(distance) > 3:
            try:
                covariance = np.linalg.inv(jacobian.T @ jacobian)
                covariance *= objective / (len(distance) - 3)
                se = np.sqrt(np.maximum(0.0, np.diag(covariance)))
            except np.linalg.LinAlgError:
                pass
    lambda_ci = (
        lambda_nm * math.exp(-1.96 * se[0]),
        lambda_nm * math.exp(1.96 * se[0]),
    ) if np.isfinite(se[0]) and se[0] < 5.0 else (float("nan"), float("nan"))
    zeta_ci = (
        zeta_mV * math.exp(-1.96 * se[1]),
        zeta_mV * math.exp(1.96 * se[1]),
    ) if np.isfinite(se[1]) and se[1] < 5.0 else (float("nan"), float("nan"))
    boundary_fraction = np.minimum(
        (result.x - lower) / (upper - lower), (upper - result.x) / (upper - lower)
    )
    reasons: list[str] = []
    if not result.success:
        reasons.append("optimizer")
    if rank != 3:
        reasons.append("jacobian_rank")
    if not np.isfinite(condition) or condition > 1e6:
        reasons.append("jacobian_condition")
    if np.any(boundary_fraction < 1e-4):
        reasons.append("parameter_boundary")
    if not np.isfinite(r2) or r2 < 0.70:
        reasons.append("r2")
    signal = float(np.max(force) - np.median(force[-max(3, len(force) // 5) :]))
    if signal < 5.0 * noise_floor:
        reasons.append("weak_signal")
    return {
        "concentration_wt_percent": concentration,
        "variant": variant,
        "model": model,
        "fit_valid": not reasons,
        "invalid_reasons": ";".join(reasons),
        "n_points": len(distance),
        "dmin_nm": dmin_nm,
        "dmax_nm": dmax_nm,
        "lambda_D_nm": lambda_nm,
        "lambda_ci95_low_nm": lambda_ci[0],
        "lambda_ci95_high_nm": lambda_ci[1],
        "zeta_magnitude_mV": zeta_mV,
        "zeta_signed_silica_mV": -zeta_mV,
        "zeta_ci95_low_mV": zeta_ci[0],
        "zeta_ci95_high_mV": zeta_ci[1],
        "baseline_pN": baseline,
        "r2": r2,
        "weighted_objective": objective,
        "noise_floor_pN": noise_floor,
        "jacobian_rank": rank,
        "jacobian_condition": condition,
        "optimizer_success": bool(result.success),
        "optimizer_status": int(result.status),
        "optimizer_nfev": int(result.nfev),
        "optimizer_message": str(result.message),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def map_pixel_from_index(
    point_index: int, columns: int, back_and_forth: bool
) -> tuple[int, int]:
    """Convert acquisition index to physical raster row and column."""

    index = int(point_index)
    width = int(columns)
    if index < 0 or width <= 0:
        raise ValueError("map index and grid width must be positive")
    row = index // width
    scan_column = index % width
    column = (
        width - 1 - scan_column
        if back_and_forth and row % 2 == 1
        else scan_column
    )
    return row, column


def _far_field_curve_record(
    source: SourceData,
    curve: RawCurve,
    acquisition_order: int,
    speed_um_per_s: float,
    speed_label_um_per_s: float,
) -> dict:
    fit = curve.far_field_fit
    row: int | None = None
    column: int | None = None
    if (
        source.source_type == "map"
        and curve.point_index is not None
        and source.map_grid_i is not None
    ):
        row, column = map_pixel_from_index(
            curve.point_index, source.map_grid_i, source.map_back_and_forth
        )
    record = {
        "concentration_wt_percent": source.concentration_wt_percent,
        "source": str(source.path.relative_to(ROOT)),
        "source_type": source.source_type,
        "timestamp": source.timestamp.isoformat(),
        "source_speed_um_per_s": speed_um_per_s,
        "speed_label_um_per_s": speed_label_um_per_s,
        "acquisition_order": acquisition_order,
        "point_index": curve.point_index,
        "grid_row": row,
        "grid_column": column,
        "grid_i": source.map_grid_i,
        "grid_j": source.map_grid_j,
        "far_field_fit_valid": fit is not None,
        "baseline_applied": "robust_linear_always" if fit is not None else "none",
    }
    if fit is None:
        return record
    sensitivity = source.sensitivity_used_m_per_V
    force_slope_N_per_m = (
        SPRING_CONSTANT_N_PER_M * sensitivity * fit.slope_V_per_m
    )
    span_m = fit.height_max_m - fit.height_min_m
    residual_noise_pN = (
        SPRING_CONSTANT_N_PER_M
        * sensitivity
        * fit.residual_mad_V
        * 1e12
    )
    drift_span_pN = abs(force_slope_N_per_m) * span_m * 1e12
    if np.isfinite(residual_noise_pN) and residual_noise_pN > 0.0:
        drift_to_noise = drift_span_pN / residual_noise_pN
    else:
        drift_to_noise = float("nan")
    record.update(
        {
            "far_field_fraction": FAR_FIELD_FRACTION,
            "far_field_n_points": fit.n_points,
            "far_field_n_inliers": fit.n_inliers,
            "far_field_height_min_um": fit.height_min_m * 1e6,
            "far_field_height_max_um": fit.height_max_m * 1e6,
            "far_field_height_span_nm": span_m * 1e9,
            "slope_V_per_m": fit.slope_V_per_m,
            "slope_mV_per_um": fit.slope_V_per_m * 1e-3,
            "slope_force_N_per_m": force_slope_N_per_m,
            "slope_pN_per_100nm": force_slope_N_per_m * 1e5,
            "far_field_intercept_V": fit.intercept_V,
            "far_field_r2": fit.r2,
            "residual_noise_mV": fit.residual_mad_V * 1e3,
            "residual_noise_pN": residual_noise_pN,
            "absolute_drift_span_pN": drift_span_pN,
            "drift_to_noise_ratio": drift_to_noise,
            "linear_drift_resolved_above_2noise": bool(
                np.isfinite(drift_to_noise) and drift_to_noise > 2.0
            ),
        }
    )
    return record


def _finite_correlation(first: np.ndarray, second: np.ndarray) -> float:
    x = np.asarray(first, dtype=np.float64)
    y = np.asarray(second, dtype=np.float64)
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    if x.size < 3 or np.std(x) <= 0.0 or np.std(y) <= 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _robust_spatial_plane(
    columns: np.ndarray, rows: np.ndarray, values: np.ndarray
) -> tuple[float, float, float, float, int]:
    x = np.asarray(columns, dtype=np.float64)
    y = np.asarray(rows, dtype=np.float64)
    z = np.asarray(values, dtype=np.float64)
    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x = x[finite]
    y = y[finite]
    z = z[finite]
    if z.size < 8:
        return (float("nan"),) * 4 + (0,)
    x_center = float(np.mean(x))
    y_center = float(np.mean(y))
    design = np.column_stack(
        [np.ones(z.size), x - x_center, y - y_center]
    )
    keep = np.ones(z.size, dtype=bool)
    coefficients = np.full(3, np.nan)
    rank = 0
    for _ in range(4):
        coefficients, _, rank, _ = np.linalg.lstsq(
            design[keep], z[keep], rcond=None
        )
        if rank != 3 or not np.all(np.isfinite(coefficients)):
            return (float("nan"),) * 4 + (int(np.count_nonzero(keep)),)
        residual = z - design @ coefficients
        sigma = robust_mad(residual[keep])
        if not np.isfinite(sigma) or sigma <= np.finfo(np.float64).tiny:
            break
        new_keep = np.abs(residual - np.median(residual[keep])) <= 3.5 * sigma
        if np.count_nonzero(new_keep) < 8 or np.array_equal(new_keep, keep):
            break
        keep = new_keep
    prediction = design[keep] @ coefficients
    ss_res = float(np.sum((z[keep] - prediction) ** 2))
    ss_tot = float(np.sum((z[keep] - np.mean(z[keep])) ** 2))
    r2 = float("nan") if ss_tot <= 0.0 else 1.0 - ss_res / ss_tot
    return (
        float(coefficients[0]),
        float(coefficients[1]),
        float(coefficients[2]),
        r2,
        int(np.count_nonzero(keep)),
    )


def build_far_field_drift_diagnostics(
    sources: list[SourceData], prepared: list[PreparedSource]
) -> tuple[list[dict], list[dict]]:
    prepared_by_path = {item.source.path: item for item in prepared}
    curve_rows: list[dict] = []
    for source in sources:
        item = prepared_by_path[source.path]
        for acquisition_order, curve in enumerate(source.curves):
            curve_rows.append(
                _far_field_curve_record(
                    source,
                    curve,
                    acquisition_order,
                    item.speed_um_per_s,
                    item.speed_label_um_per_s,
                )
            )

    summaries: list[dict] = []
    for source in sources:
        if source.source_type != "map":
            continue
        relative = str(source.path.relative_to(ROOT))
        subset = [
            row
            for row in curve_rows
            if row["source"] == relative
            and row.get("far_field_fit_valid")
            and np.isfinite(float(row.get("slope_pN_per_100nm", np.nan)))
        ]
        values = np.asarray(
            [float(row["slope_pN_per_100nm"]) for row in subset], dtype=float
        )
        grid_rows = np.asarray([int(row["grid_row"]) for row in subset], dtype=float)
        grid_columns = np.asarray(
            [int(row["grid_column"]) for row in subset], dtype=float
        )
        acquisition = np.asarray(
            [int(row["acquisition_order"]) for row in subset], dtype=float
        )
        fit_r2 = np.asarray([float(row["far_field_r2"]) for row in subset])
        resolved = np.asarray(
            [bool(row["linear_drift_resolved_above_2noise"]) for row in subset]
        )
        _, gradient_column, gradient_row, plane_r2, plane_inliers = (
            _robust_spatial_plane(grid_columns, grid_rows, values)
        )
        even_values = values[(grid_rows.astype(int) % 2) == 0]
        odd_values = values[(grid_rows.astype(int) % 2) == 1]
        even_odd_difference = (
            float(np.median(even_values) - np.median(odd_values))
            if even_values.size and odd_values.size
            else float("nan")
        )
        item = prepared_by_path[source.path]
        expected = int(source.map_grid_i * source.map_grid_j)
        summaries.append(
            {
                "concentration_wt_percent": source.concentration_wt_percent,
                "source": relative,
                "timestamp": source.timestamp.isoformat(),
                "speed_um_per_s": item.speed_um_per_s,
                "speed_label_um_per_s": item.speed_label_um_per_s,
                "grid_i": source.map_grid_i,
                "grid_j": source.map_grid_j,
                "expected_pixels": expected,
                "valid_slope_pixels": int(values.size),
                "coverage_fraction": values.size / expected,
                "median_slope_pN_per_100nm": float(np.median(values)),
                "slope_mad_pN_per_100nm": robust_mad(values),
                "slope_p05_pN_per_100nm": float(np.quantile(values, 0.05)),
                "slope_p95_pN_per_100nm": float(np.quantile(values, 0.95)),
                "median_far_field_r2": float(np.nanmedian(fit_r2)),
                "resolved_drift_fraction": float(np.mean(resolved)),
                "correlation_with_acquisition_order": _finite_correlation(
                    values, acquisition
                ),
                "correlation_with_grid_row": _finite_correlation(
                    values, grid_rows
                ),
                "correlation_with_grid_column": _finite_correlation(
                    values, grid_columns
                ),
                "spatial_plane_gradient_column_pN_per_100nm_per_pixel": gradient_column,
                "spatial_plane_gradient_row_pN_per_100nm_per_pixel": gradient_row,
                "spatial_plane_gradient_magnitude_pN_per_100nm_per_pixel": float(
                    np.hypot(gradient_column, gradient_row)
                ),
                "spatial_plane_r2": plane_r2,
                "spatial_plane_inliers": plane_inliers,
                "even_minus_odd_row_median_pN_per_100nm": even_odd_difference,
            }
        )
    return curve_rows, summaries


def _nice_symmetric_limit(values: np.ndarray, quantile: float = 0.99) -> float:
    finite = np.abs(np.asarray(values, dtype=np.float64))
    finite = finite[np.isfinite(finite)]
    if not finite.size:
        return 1.0
    raw = max(float(np.quantile(finite, quantile)), np.finfo(float).eps)
    decade = 10.0 ** math.floor(math.log10(raw))
    return float(math.ceil(raw / decade * 2.0) / 2.0 * decade)


def plot_far_field_slope_maps(
    sources: list[SourceData],
    curve_rows: list[dict],
    summaries: list[dict],
    figures_dir: Path,
) -> float:
    map_sources = sorted(
        [source for source in sources if source.source_type == "map"],
        key=lambda source: source.timestamp,
    )
    rows_by_source: dict[str, list[dict]] = {}
    for row in curve_rows:
        if row["source_type"] == "map" and row.get("far_field_fit_valid"):
            rows_by_source.setdefault(str(row["source"]), []).append(row)
    summary_by_source = {str(row["source"]): row for row in summaries}
    all_slopes = np.asarray(
        [
            float(row["slope_pN_per_100nm"])
            for row in curve_rows
            if row["source_type"] == "map"
            and row.get("far_field_fit_valid")
            and np.isfinite(float(row.get("slope_pN_per_100nm", np.nan)))
        ],
        dtype=float,
    )
    global_limit = _nice_symmetric_limit(all_slopes, 0.99)
    cmap = plt.colormaps["coolwarm"].copy()
    cmap.set_bad("#303030")

    individual_dir = figures_dir / "far_field_slope_maps"
    individual_dir.mkdir(parents=True, exist_ok=True)
    matrices: dict[str, np.ndarray] = {}
    for source in map_sources:
        relative = str(source.path.relative_to(ROOT))
        matrix = np.full(
            (int(source.map_grid_j), int(source.map_grid_i)),
            np.nan,
            dtype=float,
        )
        for row in rows_by_source.get(relative, []):
            grid_row = int(row["grid_row"])
            grid_column = int(row["grid_column"])
            if 0 <= grid_row < matrix.shape[0] and 0 <= grid_column < matrix.shape[1]:
                matrix[grid_row, grid_column] = float(row["slope_pN_per_100nm"])
        matrices[relative] = matrix
        local_limit = _nice_symmetric_limit(matrix, 0.99)
        summary = summary_by_source[relative]
        figure, axis = plt.subplots(figsize=(6.4, 5.5), constrained_layout=True)
        image = axis.imshow(
            matrix,
            origin="lower",
            interpolation="nearest",
            cmap=cmap,
            vmin=-local_limit,
            vmax=local_limit,
            extent=(-0.5, matrix.shape[1] - 0.5, -0.5, matrix.shape[0] - 0.5),
        )
        axis.set_xticks([0, 5, 10, 15])
        axis.set_yticks([0, 5, 10, 15])
        axis.set_xlabel("physical x pixel")
        axis.set_ylabel("physical y pixel")
        axis.set_title(
            f"{source.concentration_wt_percent} wt% | "
            f"{summary['speed_um_per_s']:.3g} um/s | "
            f"{source.timestamp.strftime('%H:%M:%S')}\n"
            "far-field dF/d(measuredHeight); dark gray = unacquired"
        )
        figure.colorbar(
            image,
            ax=axis,
            extend="both",
            label="slope (pN per 100 nm scanner motion)",
        )
        axis.text(
            0.01,
            -0.16,
            f"median={summary['median_slope_pN_per_100nm']:.2f}, "
            f"MAD={summary['slope_mad_pN_per_100nm']:.2f}, "
            f"coverage={summary['coverage_fraction']:.1%}, "
            f"plane R2={summary['spatial_plane_r2']:.3f}",
            transform=axis.transAxes,
            fontsize=8,
        )
        filename = (
            f"{source.concentration_wt_percent:02d}wt_"
            f"{source.timestamp.strftime('%H%M%S_%f')}_far_field_slope.png"
        )
        figure.savefig(individual_dir / filename, dpi=220)
        plt.close(figure)

    concentrations = sorted(EPSILON_R)
    maximum_maps = max(
        sum(source.concentration_wt_percent == concentration for source in map_sources)
        for concentration in concentrations
    )
    figure, axes = plt.subplots(
        len(concentrations),
        maximum_maps,
        figsize=(3.15 * maximum_maps, 3.05 * len(concentrations)),
        constrained_layout=True,
        squeeze=False,
    )
    last_image = None
    for row_index, concentration in enumerate(concentrations):
        subset = [
            source
            for source in map_sources
            if source.concentration_wt_percent == concentration
        ]
        for column_index, axis in enumerate(axes[row_index]):
            if column_index >= len(subset):
                axis.axis("off")
                continue
            source = subset[column_index]
            relative = str(source.path.relative_to(ROOT))
            matrix = matrices[relative]
            summary = summary_by_source[relative]
            last_image = axis.imshow(
                matrix,
                origin="lower",
                interpolation="nearest",
                cmap=cmap,
                vmin=-global_limit,
                vmax=global_limit,
            )
            axis.set_xticks([0, 5, 10, 15])
            axis.set_yticks([0, 5, 10, 15])
            axis.tick_params(labelsize=7)
            axis.set_title(
                f"{source.timestamp.strftime('%H:%M:%S')} | "
                f"{summary['speed_um_per_s']:.2g} um/s\n"
                f"med {summary['median_slope_pN_per_100nm']:.1f}; "
                f"plane R2 {summary['spatial_plane_r2']:.2f}",
                fontsize=8,
            )
            if column_index == 0:
                axis.set_ylabel(f"{concentration} wt%\ny pixel")
            if row_index == len(concentrations) - 1:
                axis.set_xlabel("x pixel")
    if last_image is not None:
        figure.colorbar(
            last_image,
            ax=axes.ravel().tolist(),
            shrink=0.55,
            extend="both",
            label="dF/d(measuredHeight) (pN per 100 nm); common 99% scale",
        )
    figure.suptitle(
        "Per-pixel far-field linear drift slopes; serpentine rows restored to physical x",
        fontsize=14,
    )
    figure.savefig(figures_dir / "far_field_slope_map_montage.png", dpi=220)
    plt.close(figure)
    return global_limit


def plot_sensitivity(sources: list[SourceData], path: Path) -> None:
    origin = min(source.timestamp for source in sources)
    figure, axis = plt.subplots(figsize=(10.5, 5.5), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, 5))
    for color, concentration in zip(colors, sorted(EPSILON_R)):
        subset = [s for s in sources if s.concentration_wt_percent == concentration]
        time = np.array([(s.timestamp - origin).total_seconds() / 60 for s in subset])
        used = np.array([s.sensitivity_used_m_per_V * 1e9 for s in subset])
        stored = np.array([s.stored_sensitivity_m_per_V * 1e9 for s in subset])
        axis.plot(time, used, "o-", color=color, label=f"{concentration} wt% raw consensus")
        axis.scatter(time, stored, marker="x", color=color, alpha=0.65)
        for x, y, source in zip(time, used, subset):
            if source.sensitivity_method.startswith("time_interpolation"):
                axis.scatter([x], [y], facecolors="none", edgecolors="black", s=65)
    axis.set_xlabel("Minutes since first acquisition")
    axis.set_ylabel("InvOLS / sensitivity (nm V$^{-1}$)")
    axis.set_title("Concentration-wide raw hard-contact sensitivity; crosses are ignored JPK values")
    axis.grid(alpha=0.25)
    axis.legend(ncol=2, fontsize=8)
    figure.savefig(path, dpi=220)
    plt.close(figure)


def plot_force_fits(
    primary_equilibrium: dict[int, list[dict]],
    primary_fits: dict[int, dict],
    path: Path,
) -> None:
    figure, axes = plt.subplots(3, 2, figsize=(12, 13), constrained_layout=True)
    axes_flat = axes.ravel()
    for axis, concentration in zip(axes_flat, sorted(EPSILON_R)):
        records = primary_equilibrium[concentration]
        distance = np.array([float(row["distance_nm"]) for row in records])
        equilibrium = np.array([float(row["equilibrium_force_pN"]) for row in records])
        sigma = np.array([float(row["equilibrium_sigma_pN"]) for row in records])
        for speed, key, color in (
            (1, "speed_1_force_pN", "#4c78a8"),
            (2, "speed_2_force_pN", "#f58518"),
            (4, "speed_4_force_pN", "#54a24b"),
        ):
            values = np.array([float(row[key]) for row in records])
            axis.plot(distance, values, color=color, alpha=0.45, lw=1.0, label=f"{speed} um/s")
        axis.errorbar(
            distance,
            equilibrium,
            yerr=sigma,
            fmt=".",
            ms=3,
            color="black",
            alpha=0.75,
            label="v -> 0 intercept",
        )
        fit = primary_fits.get(concentration)
        if fit and "lambda_D_nm" in fit:
            fit_distance = np.linspace(5.0, 300.0, 400)
            prediction = total_equilibrium_force_pN(
                fit_distance,
                float(fit["lambda_D_nm"]),
                float(fit["zeta_magnitude_mV"]),
                float(fit["baseline_pN"]),
                EPSILON_R[concentration],
                "nonlinear_pb_derjaguin",
            )
            axis.plot(fit_distance, prediction, color="#b22222", lw=2.0, label="nonlinear PB fit")
        axis.axhline(0.0, color="0.6", lw=0.8)
        axis.set_xlim(0, 300)
        axis.set_title(f"{concentration} wt% glycerol")
        axis.set_xlabel("Contact-aligned separation (nm)")
        axis.set_ylabel("Force (pN)")
        axis.grid(alpha=0.2)
        axis.legend(fontsize=7)
    axes_flat[-1].axis("off")
    figure.savefig(path, dpi=220)
    plt.close(figure)


def plot_summary(primary_fits: dict[int, dict], path: Path) -> None:
    concentration = np.array(sorted(primary_fits), dtype=float)
    lambda_D = np.array([primary_fits[int(c)].get("lambda_D_nm", np.nan) for c in concentration], dtype=float)
    zeta = np.array([primary_fits[int(c)].get("zeta_signed_silica_mV", np.nan) for c in concentration], dtype=float)
    lambda_low = np.array([primary_fits[int(c)].get("lambda_ci95_low_nm", np.nan) for c in concentration], dtype=float)
    lambda_high = np.array([primary_fits[int(c)].get("lambda_ci95_high_nm", np.nan) for c in concentration], dtype=float)
    zeta_low_mag = np.array([primary_fits[int(c)].get("zeta_ci95_low_mV", np.nan) for c in concentration], dtype=float)
    zeta_high_mag = np.array([primary_fits[int(c)].get("zeta_ci95_high_mV", np.nan) for c in concentration], dtype=float)
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), constrained_layout=True)
    axes[0].errorbar(
        concentration,
        lambda_D,
        yerr=np.vstack([lambda_D - lambda_low, lambda_high - lambda_D]),
        fmt="o-",
        capsize=3,
    )
    axes[0].set_xlabel("Glycerol mass concentration (wt%)")
    axes[0].set_ylabel("Debye length (nm)")
    axes[0].grid(alpha=0.25)
    axes[1].errorbar(
        concentration,
        zeta,
        yerr=np.vstack([zeta_high_mag + zeta, -zeta - zeta_low_mag]),
        fmt="o-",
        capsize=3,
        color="#b22222",
    )
    axes[1].set_xlabel("Glycerol mass concentration (wt%)")
    axes[1].set_ylabel("Silica-assigned zeta (mV)")
    axes[1].grid(alpha=0.25)
    figure.savefig(path, dpi=220)
    plt.close(figure)


def self_test() -> dict[str, float]:
    synthetic_height = np.linspace(1.0e-6, 0.0, 1000)
    synthetic_slope = 2.5e4
    synthetic_voltage = synthetic_slope * synthetic_height + 0.031
    synthetic_fit = fit_far_field_drift(synthetic_height, synthetic_voltage)
    drift_slope_error = abs(
        synthetic_fit.slope_V_per_m / synthetic_slope - 1.0
    )
    corrected_far = synthetic_voltage - baseline_voltage(
        synthetic_height, synthetic_voltage, synthetic_fit
    )
    drift_subtraction_residual = float(
        np.max(np.abs(corrected_far[: synthetic_fit.n_points]))
    )
    if drift_slope_error > 1e-12 or drift_subtraction_residual > 1e-12:
        raise AssertionError("far-field linear drift recovery failed")
    expected_pixels = {
        0: (0, 0),
        15: (0, 15),
        16: (1, 15),
        31: (1, 0),
    }
    for point_index, expected in expected_pixels.items():
        if map_pixel_from_index(point_index, 16, True) != expected:
            raise AssertionError("serpentine map reconstruction failed")
    H = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0])
    surface = 0.01
    nonlinear = pb_dimensionless_G(H, surface)
    linear = 2.0 * surface**2 / (np.exp(H) + 1.0)
    pb_linear_limit_error = float(np.max(np.abs(nonlinear / linear - 1.0)))
    if pb_linear_limit_error > 2e-3:
        raise AssertionError(f"PB linear limit error {pb_linear_limit_error}")
    large_H = np.array([20.0, 30.0, 50.0])
    surface = 3.0
    nonlinear_far = pb_dimensionless_G(large_H, surface)
    far = 32.0 * math.tanh(surface / 4.0) ** 2 * np.exp(-large_H)
    pb_far_error = float(np.max(np.abs(nonlinear_far / far - 1.0)))
    if pb_far_error > 0.03:
        raise AssertionError(f"PB far-field error {pb_far_error}")
    viscosity_water = cheng_viscosity_mPa_s(0.0, TEMPERATURE_C)
    viscosity_glycerol = cheng_viscosity_mPa_s(1.0, TEMPERATURE_C)
    if not (0.8 < viscosity_water < 1.0 and 700.0 < viscosity_glycerol < 1000.0):
        raise AssertionError("Cheng endpoint viscosities are implausible")
    force = edl_force_pN(np.array([20.0, 100.0]), 50.0, 60.0, 78.5, "nonlinear_pb_derjaguin")
    if not np.all(np.isfinite(force)) or not np.all(force > 0.0) or force[0] <= force[1]:
        raise AssertionError("EDL force positivity/monotonicity failed")
    return {
        "far_field_synthetic_slope_relative_error": drift_slope_error,
        "far_field_synthetic_subtraction_max_abs_V": drift_subtraction_residual,
        "pb_linear_limit_max_relative_error": pb_linear_limit_error,
        "pb_far_asymptote_max_relative_error": pb_far_error,
        "cheng_water_viscosity_mPa_s": viscosity_water,
        "cheng_glycerol_viscosity_mPa_s": viscosity_glycerol,
    }


def build_report(
    sources: list[SourceData],
    prepared: list[PreparedSource],
    fits: list[dict],
    primary_fits: dict[int, dict],
    mixture: list[dict],
    drift_summaries: list[dict],
    self_checks: dict[str, float],
) -> str:
    lines = [
        "# 20-08-26 silica sphere-plane surface-force fit",
        "",
        "## Result status",
        "",
        "This analysis decodes raw JPK `vDeflection` in volts and raw scanner readback in metres. The JPK header InvOLS and force channels were retained only for comparison and were never used to calculate deflection or force. Force was rebuilt with the cantilever-1 calibration `k = 0.2969899087 N/m` at 25.6 C.",
        "",
        "The operational primary result uses the maps acquired over the same area: at each separation, the measured 1/2/4 um/s map medians are extrapolated to zero approach speed. The 25 independent curves form a separate formal cohort and a cross-check, rather than being allowed to outvote a map merely because they are stored as five files. The zero-speed map intercept is fitted to a nonlinear, equal-constant-potential, symmetric 1:1 Poisson-Boltzmann planar solution and converted to sphere-plane force with the Derjaguin approximation. This avoids using the reference repository's equal-sphere prefactor.",
        "",
        "Important experimental consistency finding: the three sequential maps show force decreasing with approach speed over 20-100 nm, opposite to the positive no-slip lubrication slope. Therefore the numerical PB fits below are reproducible operational fits, but the speed extrapolation cannot be claimed as a validated hydrodynamic correction or a unique equilibrium measurement. Sequential surface-state/contact-alignment changes remain confounded with speed.",
        "",
        "| glycerol wt% | Debye length (nm) | common-potential magnitude (mV) | silica-assigned sign (mV) | R2 | validity |",
        "|---:|---:|---:|---:|---:|:---|",
    ]
    for concentration in sorted(primary_fits):
        fit = primary_fits[concentration]
        if "lambda_D_nm" not in fit:
            lines.append(f"| {concentration} | n/a | n/a | n/a | n/a | invalid: {fit.get('invalid_reasons','')} |")
            continue
        numerical = "numerically valid" if fit.get("fit_valid") else f"numerically invalid: {fit.get('invalid_reasons','')}"
        speed = "speed consistent" if fit.get("map_speed_consistency_pass") else "map speed inconsistency"
        validity = f"{numerical}; {speed}"
        lines.append(
            f"| {concentration} | {fit['lambda_D_nm']:.3g} "
            f"[{fit['lambda_ci95_low_nm']:.3g}, {fit['lambda_ci95_high_nm']:.3g}] | "
            f"{fit['zeta_magnitude_mV']:.3g} [{fit['zeta_ci95_low_mV']:.3g}, {fit['zeta_ci95_high_mV']:.3g}] | "
            f"{fit['zeta_signed_silica_mV']:.3g} | {fit['r2']:.4f} | {validity} |"
        )
    lines += [
        "",
        "The confidence intervals above are local Jacobian intervals conditional on the chosen contact alignment, common-potential boundary condition, dielectric constants, radius and Hamaker constant. They are much narrower than the experimental systematic spread and must not be quoted alone.",
        "",
        "| wt% | Debye-length systematic range (nm) | potential-magnitude systematic range (mV) |",
        "|---:|---:|---:|",
    ]
    for concentration in sorted(primary_fits):
        fit = primary_fits[concentration]
        if "systematic_lambda_min_nm" in fit:
            lines.append(
                f"| {concentration} | {fit['systematic_lambda_min_nm']:.3g}-{fit['systematic_lambda_max_nm']:.3g} | "
                f"{fit['systematic_zeta_min_mV']:.3g}-{fit['systematic_zeta_max_mV']:.3g} |"
            )
    lines += [
        "",
        "These ranges span nonlinear map fit-window, all-source and independent-force-only variants where available. They are sensitivity ranges, not probability intervals. The 10 wt% Debye length is especially weak because many traces lack a resolved hard-contact interval and the force signal is small.",
        "",
        "Force between identical surfaces is even in the potential, so AFM normal-force data determine `|zeta|`, not its sign. The negative sign in the table is assigned from silica surface chemistry, not inferred from the force curve. More precisely, the fitted quantity is the constant-potential PB boundary potential; calling it zeta additionally assumes that this potential represents the electrokinetic slipping-plane potential.",
        "",
        "## Model and derivation",
        "",
        "For dimensionless potential `u=e psi/(k_B T)` and dimensionless plate spacing `H=kappa D`, the equal-potential planar boundary-value problem is",
        "",
        "`u'' = sinh(u)`, `u(+/-H/2)=u_s`, and `u'(0)=0`.",
        "",
        "Its first integral gives",
        "",
        "`H/2 = integral_(u_m)^(u_s) du / sqrt(2[cosh(u)-cosh(u_m)])`.",
        "",
        "At the symmetry plane the Maxwell-field term is zero, so `Pi = eps (k_B T/e)^2 kappa^2 [cosh(u_m)-1]`. The sphere-plane Derjaguin force used here is therefore",
        "",
        "`F_EDL(D)=2 pi R eps (k_B T/e)^2 kappa integral_(kappa D)^infinity [cosh(u_m(H))-1] dH`.",
        "",
        "The fixed van der Waals term is `F_vdW=-A_H R/(6D^2)` with `A_H=2.4e-21 J`. The fitted common potential and inverse screening length are shared by the silica sphere and silica plane within each concentration.",
        "",
        "In the low-potential limit this implementation reduces to the sphere-plane same-potential HHF expression",
        "",
        "`F_EDL = 4 pi R eps kappa zeta^2 / [exp(kappa D)+1]`.",
        "",
        "For two equal spheres the Derjaguin radius would be `R/2`; therefore copying the reference paper's equal-sphere formula would understate the present sphere-plane force by a factor of two for the same physical sphere radius.",
        "",
        "Primary sources: [Hogg-Healy-Fuerstenau, linear Debye-Huckel sphere interactions](https://pubs.rsc.org/en/content/articlehtml/1966/tf/tf9666201638); [Stankovich and Carnie, nonlinear PB sphere/plate and Derjaguin accuracy](https://pubs.acs.org/doi/10.1021/la950384k); [Polat and Polat, nonlinear PB for arbitrary same-sign plate potentials](https://doi.org/10.1016/j.jcis.2009.09.008); [Liu and Li, AFM zeta workflow used by the reference repository](https://doi.org/10.1016/j.jcis.2020.05.061).",
        "",
        "## Sensitivity reconstruction",
        "",
        "Each usable approach trace was searched near its terminal end for a contiguous rigid-contact interval. A robust line `V=a+b h` gives `InvOLS=1/|b|`. First, each source supplies a robust contact-slope anchor. Then all consistent source anchors at the same glycerol concentration are combined with equal source weighting into one concentration-wide sensitivity. This keeps the force calibration independent of approach speed and uses all available raw hard-contact data; the embedded JPK sensitivity is never substituted.",
        "",
        "| wt% | sources | raw curves | source anchors | retained anchors | common sensitivity used (nm/V) |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for concentration in sorted(EPSILON_R):
        subset = [s for s in sources if s.concentration_wt_percent == concentration]
        values = [s.sensitivity_used_m_per_V * 1e9 for s in subset]
        retained = sum(
            np.isfinite(s.sensitivity_anchor_m_per_V)
            and abs(s.sensitivity_anchor_m_per_V - s.sensitivity_used_m_per_V) <= 15e-9
            for s in subset
        )
        lines.append(
            f"| {concentration} | {len(subset)} | {sum(len(s.curves) for s in subset)} | "
            f"{sum(np.isfinite(s.sensitivity_anchor_m_per_V) for s in subset)} | "
            f"{retained} | {np.median(values):.3f} |"
        )
    lines += [
        "",
        "A curve lacking a resolvable hard-contact interval cannot independently establish its contact plane. For those curves, the terminal `h+delta` value is used as the mechanical-contact estimate and the curve is flagged as endpoint-aligned. This limitation is particularly important for 10 wt% and is visible in the source-subset/window sensitivity results.",
        "",
        "## Far-field linear drift removal and spatial QC",
        "",
        "Before contact finding or force conversion, every approach trace uses its initial 20% of samples (200 of 1000 for the present data; minimum 80) as the far field. A four-pass robust line `V_ff(h)=a h+b` is fitted in raw `vDeflection` versus raw `measuredHeight`, with 3.5-MAD residual rejection. The full fitted line is subtracted from the entire trace: `V_corrected=V-V_ff`. This is the linear drift-removal step used by all sensitivity, force, speed-extrapolation and PB results; the later 700-880 nm source correction removes only a residual constant offset.",
        "",
        "The slope below is reported after applying the concentration-wide calculated sensitivity and cantilever-1 spring constant, as `dF/d(measuredHeight)` in pN per 100 nm scanner motion. The sign is therefore tied to the scanner-height coordinate, not acquisition time. `resolved` means the fitted change across the far window exceeds twice the post-fit residual MAD; it is diagnostic only, because the declared linear operator is applied to every valid trace.",
        "",
        "Every map declares a 16x16, back-and-forth raster. Odd acquisition rows are reversed back into physical x before plotting. Dark-gray pixels in partial 30 wt% maps are unacquired records, not interpolated values.",
        "",
        "| wt% | start time | speed (um/s) | pixels | median slope +/- MAD (pN/100 nm) | resolved | corr. acquisition order | spatial-plane R2 |",
        "|---:|:---|---:|---:|---:|---:|---:|---:|",
    ]
    for summary in drift_summaries:
        lines.append(
            f"| {summary['concentration_wt_percent']} | "
            f"{str(summary['timestamp'])[11:19]} | "
            f"{summary['speed_um_per_s']:.3g} | "
            f"{summary['valid_slope_pixels']}/{summary['expected_pixels']} | "
            f"{summary['median_slope_pN_per_100nm']:.3g} +/- "
            f"{summary['slope_mad_pN_per_100nm']:.3g} | "
            f"{summary['resolved_drift_fraction']:.1%} | "
            f"{summary['correlation_with_acquisition_order']:.3f} | "
            f"{summary['spatial_plane_r2']:.3f} |"
        )
    lines += [
        "",
        "Per-pixel values, fit quality, residual noise, drift/noise ratio, raster coordinates and plane-gradient summaries are retained in `far_field_drift_by_curve.csv` and `far_field_drift_by_map.csv`; no visual clipping changes the saved numerical values.",
        "",
        "## Speed dependence and liquid properties",
        "",
        "At each 5 nm bin, the same-area map medians at the measured approximately 1, 2 and 4 um/s speeds are regressed linearly to zero speed. This is an operational intercept, not a validated lubrication subtraction. The measured slope is checked against the positive no-slip prediction `6 pi eta R^2/D` using the Cheng glycerol-water viscosity correlation.",
        "",
        "| wt% | epsilon_r | viscosity (mPa s) | map slope / lubrication slope | map positive-slope fraction | independent-force slope / lubrication slope |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in mixture:
        lines.append(
            f"| {row['concentration_wt_percent']} | {row['epsilon_r']:.3f} | "
            f"{row['viscosity_mPa_s']:.4f} | "
            f"{row['map_median_observed_to_lubrication_slope_ratio_20_100nm']:.3g} | "
            f"{row['map_positive_speed_slope_fraction_20_100nm']:.3f} | "
            f"{row['force_median_observed_to_lubrication_slope_ratio_20_100nm']:.3g} |"
        )
    lines += [
        "",
        "All map medians fail the expected positive lubrication direction. Some independent-force cohorts instead have large positive slopes, so the two acquisition modes are not mutually interchangeable replicates. Viscosity source: [Cheng, Formula for the Viscosity of a Glycerol-Water Mixture](https://doi.org/10.1021/ie071349z). The dielectric constants (78.5, 76.0, 73.9, 71.4 and 68.9) are rounded from the 25 C, 0.57 MHz measurements tabulated at 0.00, 9.88, 20.33, 30.19 and 39.67 wt% in [Physical Properties of Glycerine and Its Solutions, Table 27](https://www.cleaninginstitute.org/sites/default/files/research-pdfs/Physical_properties_of_glycerine_and_its_solutions.pdf), and agree with the reference repository. Modern primary measurements over composition and 10-50 C are reported by [Behrends et al.](https://pubmed.ncbi.nlm.nih.gov/16626219/). The PB temperature itself is the measured 25.6 C; the closest tabulated permittivity temperature is 25 C.",
        "",
        "## Numerical and provenance checks",
        "",
        f"- Parsed {len(sources)} raw sources and {sum(len(s.curves) for s in sources)} usable approach curves; partial/short terminal records skipped: {sum(s.skipped_curves for s in sources)}.",
        f"- Far-field linear-drift synthetic slope relative error: `{self_checks['far_field_synthetic_slope_relative_error']:.3e}`; maximum residual after subtraction: `{self_checks['far_field_synthetic_subtraction_max_abs_V']:.3e} V`.",
        f"- Nonlinear PB to analytic linear-PB limit maximum relative error: `{self_checks['pb_linear_limit_max_relative_error']:.3e}`.",
        f"- Nonlinear PB far-field asymptote maximum relative error: `{self_checks['pb_far_asymptote_max_relative_error']:.3e}`.",
        f"- Cheng endpoint check at 25.6 C: water `{self_checks['cheng_water_viscosity_mPa_s']:.6g} mPa s`, glycerol `{self_checks['cheng_glycerol_viscosity_mPa_s']:.6g} mPa s`.",
        "- Every input ZIP passed CRC during decoding; SHA-256 values are in `input_manifest.csv`.",
        "- Optimizer termination alone is not accepted: fit validity also checks Jacobian rank/condition, parameter bounds, R2 and signal-to-noise.",
        "",
        "## Files",
        "",
        "- `sensitivity_by_source.csv`: raw source anchors, retained-anchor flags, common calculated sensitivity and ignored embedded calibration.",
        "- `far_field_drift_by_curve.csv`: fitted slope, R2, residual noise and drift/noise for every approach trace.",
        "- `far_field_drift_by_map.csv`: 16x16 map-level spatial-gradient and acquisition-order diagnostics.",
        "- `source_binned_force_curves.csv`: reconstructed source-median force curves by actual approach speed.",
        "- `zero_speed_equilibrium.csv`: three-speed intercepts and speed slopes.",
        "- `fit_results.csv`: nonlinear primary fit, fit-window/source-subset checks and linear-HHF diagnostic.",
        "- `mixture_properties.csv`, `input_manifest.csv`, `provenance.json`, figures and SHA-256 manifest.",
    ]
    return "\n".join(lines) + "\n"


def run(results_dir: Path) -> None:
    self_checks = self_test()
    sources: list[SourceData] = []
    for concentration in sorted(EPSILON_R):
        directory = DATA_ROOT / str(concentration)
        paths = sorted(directory.glob("*.jpk-force")) + sorted(
            directory.glob("*.jpk-force-map")
        )
        for path in paths:
            sources.append(load_source(path, concentration))
    sources.sort(key=lambda item: item.timestamp)
    calibrate_sensitivity(sources)
    prepared = [prepare_source(source) for source in sources]

    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    sensitivity_rows: list[dict] = []
    manifest_rows: list[dict] = []
    for source in sources:
        sensitivity_rows.append(
            {
                "concentration_wt_percent": source.concentration_wt_percent,
                "source": str(source.path.relative_to(ROOT)),
                "source_type": source.source_type,
                "timestamp": source.timestamp.isoformat(),
                "raw_curves": len(source.curves),
                "skipped_curves": source.skipped_curves,
                "valid_contact_curves": source.sensitivity_valid_curves,
                "calculated_anchor_nm_per_V": source.sensitivity_anchor_m_per_V * 1e9,
                "calculated_anchor_mad_nm_per_V": source.sensitivity_anchor_mad_m_per_V * 1e9,
                "anchor_used_in_consensus": bool(
                    np.isfinite(source.sensitivity_anchor_m_per_V)
                    and abs(
                        source.sensitivity_anchor_m_per_V
                        - source.sensitivity_used_m_per_V
                    )
                    <= 15e-9
                ),
                "sensitivity_used_nm_per_V": source.sensitivity_used_m_per_V * 1e9,
                "sensitivity_method": source.sensitivity_method,
                "stored_ignored_sensitivity_nm_per_V": source.stored_sensitivity_m_per_V * 1e9,
                "stored_ignored_spring_constant_N_per_m": source.stored_spring_constant_N_per_m,
                "calibrated_spring_constant_used_N_per_m": SPRING_CONSTANT_N_PER_M,
            }
        )
        manifest_rows.append(
            {
                "source": str(source.path.relative_to(ROOT)),
                "sha256": source.sha256,
                "size_bytes": source.path.stat().st_size,
                "zip_crc_pass": True,
            }
        )
    write_csv(results_dir / "sensitivity_by_source.csv", sensitivity_rows)
    write_csv(results_dir / "input_manifest.csv", manifest_rows)

    drift_curve_rows, drift_summaries = build_far_field_drift_diagnostics(
        sources, prepared
    )
    write_csv(results_dir / "far_field_drift_by_curve.csv", drift_curve_rows)
    write_csv(results_dir / "far_field_drift_by_map.csv", drift_summaries)
    far_field_plot_limit = plot_far_field_slope_maps(
        sources, drift_curve_rows, drift_summaries, figures_dir
    )

    source_curve_rows: list[dict] = []
    for item in prepared:
        for index, distance in enumerate(BIN_CENTERS_NM):
            if not np.isfinite(item.force_pN[index]):
                continue
            source_curve_rows.append(
                {
                    "concentration_wt_percent": item.source.concentration_wt_percent,
                    "source": str(item.source.path.relative_to(ROOT)),
                    "source_type": item.source.source_type,
                    "speed_um_per_s": item.speed_um_per_s,
                    "speed_label_um_per_s": item.speed_label_um_per_s,
                    "distance_nm": distance,
                    "force_pN": item.force_pN[index],
                    "force_spread_pN": item.force_spread_pN[index],
                    "curve_count": item.counts[index],
                    "sensitivity_used_nm_per_V": item.source.sensitivity_used_m_per_V * 1e9,
                    "valid_contact_curves": item.valid_contact_curves,
                    "endpoint_contact_curves": item.endpoint_contact_curves,
                    "far_offset_pN": item.far_offset_pN,
                }
            )
    write_csv(results_dir / "source_binned_force_curves.csv", source_curve_rows)

    equilibrium_by_filter: dict[tuple[int, str], list[dict]] = {}
    equilibrium_rows: list[dict] = []
    for concentration in sorted(EPSILON_R):
        for source_filter in ("all", "map", "force"):
            records = zero_speed_curve(prepared, concentration, source_filter)
            equilibrium_by_filter[(concentration, source_filter)] = records
            equilibrium_rows.extend(records)
    write_csv(results_dir / "zero_speed_equilibrium.csv", equilibrium_rows)

    fits: list[dict] = []
    variants = (
        ("map", 10.0, 250.0, "primary"),
        ("map", 15.0, 250.0, "dmin15"),
        ("map", 20.0, 250.0, "dmin20"),
        ("map", 10.0, 150.0, "dmax150"),
        ("map", 10.0, 350.0, "dmax350"),
        ("all", 10.0, 250.0, "all_sources"),
        ("force", 10.0, 250.0, "independent_force_only"),
    )
    for concentration in sorted(EPSILON_R):
        for source_filter, dmin, dmax, variant in variants:
            fits.append(
                fit_equilibrium(
                    equilibrium_by_filter[(concentration, source_filter)],
                    concentration,
                    "nonlinear_pb_derjaguin",
                    dmin,
                    dmax,
                    variant,
                )
            )
        fits.append(
            fit_equilibrium(
                equilibrium_by_filter[(concentration, "map")],
                concentration,
                "linear_hhf_equal_potential",
                10.0,
                250.0,
                "linear_hhf_diagnostic",
            )
        )
    primary_fits = {
        int(row["concentration_wt_percent"]): row
        for row in fits
        if row.get("variant") == "primary"
        and row.get("model") == "nonlinear_pb_derjaguin"
    }

    mixture_rows: list[dict] = []
    for concentration in sorted(EPSILON_R):
        viscosity = cheng_viscosity_mPa_s(concentration / 100.0, TEMPERATURE_C)
        summaries: dict[str, tuple[float, float, float]] = {}
        for source_filter in ("all", "map", "force"):
            ratios: list[float] = []
            slopes: list[float] = []
            speed_r2: list[float] = []
            for record in equilibrium_by_filter[(concentration, source_filter)]:
                distance = float(record["distance_nm"])
                if 20.0 <= distance <= 100.0:
                    predicted = (
                        6.0
                        * np.pi
                        * (viscosity * 1e-3)
                        * PROBE_RADIUS_M**2
                        / (distance * 1e-9)
                        * 1e6
                    )
                    observed = float(record["speed_slope_pN_per_um_s"])
                    if predicted > 0.0 and np.isfinite(observed):
                        ratios.append(observed / predicted)
                        slopes.append(observed)
                        speed_r2.append(float(record["speed_regression_r2"]))
            summaries[source_filter] = (
                float(np.median(ratios)) if ratios else float("nan"),
                float(np.mean(np.asarray(slopes) > 0.0)) if slopes else float("nan"),
                float(np.nanmedian(speed_r2)) if speed_r2 else float("nan"),
            )
        mixture_rows.append(
            {
                "concentration_wt_percent": concentration,
                "temperature_C": TEMPERATURE_C,
                "epsilon_r": EPSILON_R[concentration],
                "viscosity_mPa_s": viscosity,
                "all_median_observed_to_lubrication_slope_ratio_20_100nm": summaries["all"][0],
                "all_positive_speed_slope_fraction_20_100nm": summaries["all"][1],
                "map_median_observed_to_lubrication_slope_ratio_20_100nm": summaries["map"][0],
                "map_positive_speed_slope_fraction_20_100nm": summaries["map"][1],
                "map_median_speed_regression_r2_20_100nm": summaries["map"][2],
                "force_median_observed_to_lubrication_slope_ratio_20_100nm": summaries["force"][0],
                "force_positive_speed_slope_fraction_20_100nm": summaries["force"][1],
            }
        )

    mixture_by_concentration = {
        int(row["concentration_wt_percent"]): row for row in mixture_rows
    }
    for concentration, primary in primary_fits.items():
        mixture_row = mixture_by_concentration[concentration]
        map_ratio = float(
            mixture_row[
                "map_median_observed_to_lubrication_slope_ratio_20_100nm"
            ]
        )
        map_positive_fraction = float(
            mixture_row["map_positive_speed_slope_fraction_20_100nm"]
        )
        speed_consistency = map_ratio > 0.0 and map_positive_fraction >= 0.70
        primary["map_speed_consistency_pass"] = speed_consistency
        primary["map_observed_to_lubrication_slope_ratio"] = map_ratio
        primary["map_positive_speed_slope_fraction"] = map_positive_fraction
        primary["scientific_status"] = (
            "quantitative_fit_with_speed_consistency"
            if speed_consistency
            else "operational_fit_speed_trend_not_lubrication_consistent"
        )

        nonlinear_variants = [
            row
            for row in fits
            if int(row["concentration_wt_percent"]) == concentration
            and row.get("model") == "nonlinear_pb_derjaguin"
            and "lambda_D_nm" in row
            and np.isfinite(float(row["lambda_D_nm"]))
        ]
        lambda_values = np.array(
            [float(row["lambda_D_nm"]) for row in nonlinear_variants]
        )
        zeta_values = np.array(
            [float(row["zeta_magnitude_mV"]) for row in nonlinear_variants]
        )
        primary["systematic_lambda_min_nm"] = float(np.min(lambda_values))
        primary["systematic_lambda_max_nm"] = float(np.max(lambda_values))
        primary["systematic_zeta_min_mV"] = float(np.min(zeta_values))
        primary["systematic_zeta_max_mV"] = float(np.max(zeta_values))

    write_csv(results_dir / "fit_results.csv", fits)
    write_csv(results_dir / "mixture_properties.csv", mixture_rows)

    plot_sensitivity(sources, figures_dir / "sensitivity_vs_time.png")
    plot_force_fits(
        {c: equilibrium_by_filter[(c, "map")] for c in sorted(EPSILON_R)},
        primary_fits,
        figures_dir / "speed_extrapolation_and_pb_fits.png",
    )
    plot_summary(primary_fits, figures_dir / "zeta_and_debye_summary.png")

    provenance = {
        "analysis_script": str(Path(__file__).relative_to(ROOT)),
        "analysis_script_sha256": sha256_file(Path(__file__)),
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "temperature_C": TEMPERATURE_C,
        "temperature_K": TEMPERATURE_K,
        "spring_constant_N_per_m": SPRING_CONSTANT_N_PER_M,
        "spring_constant_repeatability_sd_N_per_m": SPRING_CONSTANT_REPEAT_SD_N_PER_M,
        "probe_radius_m": PROBE_RADIUS_M,
        "probe_radius_range_m": PROBE_RADIUS_RANGE_M,
        "hamaker_J": HAMAKER_J,
        "epsilon_r": EPSILON_R,
        "epsilon_r_source": {
            "table_url": "https://www.cleaninginstitute.org/sites/default/files/research-pdfs/Physical_properties_of_glycerine_and_its_solutions.pdf",
            "table": 27,
            "temperature_C": 25.0,
            "measurement_frequency_Hz": 0.57e6,
            "tabulated_glycerol_wt_percent": [0.00, 9.88, 20.33, 30.19, 39.67],
            "tabulated_epsilon_r": [78.48, 75.98, 73.86, 71.44, 68.93],
            "modern_primary_measurement_doi": "10.1063/1.2188391",
        },
        "randomness": "none",
        "jpk_embedded_sensitivity_used": False,
        "jpk_embedded_force_used": False,
        "sensitivity_aggregation": "one equal-source-weighted raw hard-contact consensus per glycerol concentration",
        "far_field_linear_drift": {
            "coordinate": "raw measuredHeight in m",
            "signal": "raw vDeflection in V",
            "window": "initial 20 percent of approach, minimum 80 points",
            "fit": "four-pass robust line with 3.5-MAD residual rejection",
            "operator": "subtract fitted slope and intercept from the full trace before contact and force conversion",
            "slope_output_unit": "pN per 100 nm measuredHeight",
            "resolved_drift_diagnostic": "absolute fitted change across far window greater than two residual MAD",
            "map_geometry": "16x16 physical raster with serpentine acquisition rows restored",
            "montage_global_symmetric_color_limit_pN_per_100nm": far_field_plot_limit,
        },
        "speed_extrapolation_interpretation": "operational only; map speed slopes oppose no-slip lubrication",
        "self_checks": self_checks,
    }
    (results_dir / "provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    report = build_report(
        sources,
        prepared,
        fits,
        primary_fits,
        mixture_rows,
        drift_summaries,
        self_checks,
    )
    (results_dir / "REPORT.md").write_text(report, encoding="utf-8")

    artifacts = sorted(
        path
        for path in results_dir.rglob("*")
        if path.is_file() and path.name != "artifact_manifest.sha256"
    )
    manifest = "\n".join(
        f"{sha256_file(path)}  {path.relative_to(results_dir)}" for path in artifacts
    )
    (results_dir / "artifact_manifest.sha256").write_text(manifest + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--self-test", action="store_true")
    arguments = parser.parse_args()
    if arguments.self_test:
        print(json.dumps(self_test(), indent=2))
        return
    run(arguments.results_dir.resolve())


if __name__ == "__main__":
    main()

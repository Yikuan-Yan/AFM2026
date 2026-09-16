# Numerical audit: 16-09-26 low-speed analysis

## Data flow and units

`JPK int/float encoder -> raw vDeflection [V] and measuredHeight [m] -> per-curve constant baseline [V] -> contact InvOLS [nm/V] -> deflection [m] -> force [pN] and separation D [nm] -> local U_gap=-dD/dt [um/s] -> pixel bins -> map medians -> palindrome means -> double-centered model comparison`.

Force uses `k=0.174367775352 N/m`; the resulting scale is 10.1559363 nN/V. `K` has units pN nm/(um/s)^alpha; for actual-gap models alpha=1.

## Array and missing-data checks

- Force matrix shape/dtype: `(18, 64, 60)`, `float64`; finite=68215, NaN=905, Inf=0`.
- Gap-speed matrix shape/dtype: `(18, 64, 60)`, `float64`; finite=68212, NaN=908, Inf=0`.
- The 15 rejected curves contain 0 finite force bins and 0 finite gap-speed bins; both must be zero. Missing curves remain NaN and are never imputed.
- Every map retains finite map medians at 5-300 nm; each target bin has at least 61 finite pixels.

## Sensitivity and conditioning

- Contact-window batch InvOLS at 30/40/50 nm: 58.147362/58.244342/58.393970 nm/V (maximum relative shift 0.257%).
- Alternative constant-baseline windows change absolute force by a map-dependent constant. Separation and the double-centered velocity-distance interaction are invariant to that constant; absolute force-time panels remain finite-window-gauge quantities.
- Hydrodynamic fit-range sensitivity is tabulated below. Correlated 5 nm bins are not counted as independent replicates and no bin-level p-value is reported.

| range | model | K | alpha | D0 (nm) | normalized RMSE | Jacobian condition |
|---|---|---:|---:|---:|---:|---:|
| primary_20_200 | fitted_v_over_D | 32675.5 | 1.00000 | 0.00000 | 0.18015 | 1 |
| primary_20_200 | shifted_v_over_D | 61584.5 | 1.00000 | 17.81608 | 0.09031 | 7878.96 |
| primary_20_200 | empirical_power_shift | 57512.3 | 0.78390 | 17.76076 | 0.08190 | 69486.5 |
| primary_20_200 | actual_gap_shifted | 52803.5 | 1.00000 | 6.38685 | 0.08692 | 7154.6 |
| sensitivity_30_200 | fitted_v_over_D | 37653.5 | 1.00000 | 0.00000 | 0.14707 | 1 |
| sensitivity_30_200 | shifted_v_over_D | 61016 | 1.00000 | 17.30050 | 0.10508 | 7461.54 |
| sensitivity_30_200 | empirical_power_shift | 57380.1 | 0.81537 | 17.26650 | 0.09998 | 88411.2 |
| sensitivity_30_200 | actual_gap_shifted | 53232.6 | 1.00000 | 6.63768 | 0.10478 | 6648.29 |
| sensitivity_20_250 | fitted_v_over_D | 33805.2 | 1.00000 | 0.00000 | 0.18733 | 1 |
| sensitivity_20_250 | shifted_v_over_D | 62227.7 | 1.00000 | 18.14269 | 0.09288 | 6851 |
| sensitivity_20_250 | empirical_power_shift | 58065.1 | 0.77500 | 18.11020 | 0.08399 | 61974.3 |
| sensitivity_20_250 | actual_gap_shifted | 53900.9 | 1.00000 | 6.84138 | 0.09007 | 6351.21 |

All parameterized optimizations report success: failures=0. The empirical three-parameter Jacobian is full-rank but ill-conditioned, so alpha and D0 are model-shape descriptors rather than independently identified material parameters.

## Synthetic closure

- `synthetic_double_center_max_abs_error_pN` = `1.01252e-13`
- `synthetic_classical_coefficient_relative_error` = `1.47346e-16`
- `synthetic_actual_gap_coefficient_relative_error` = `1.47346e-16`
- `double_center_row_mean_max_abs_pN` = `3.78956e-14`
- `double_center_column_mean_max_abs_pN` = `5.53071e-14`

All synthetic errors are below 1e-10. Input ZIP CRC is exercised by the reused parser for every map; input and output SHA-256 manifests are written separately.

## Claim boundary

The numerics establish a reproducible velocity-distance interaction and quantify which curves are usable. They do not identify the interaction uniquely as hydrodynamic, validate the bulk no-slip prefactor, establish low-speed equivalence, or justify equilibrium PB parameters after subtraction.

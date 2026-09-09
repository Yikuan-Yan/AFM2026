# 31-08-26: 256 independent colloid-colloid fits

## Scope

Map2–map5 each contribute all 64 original signed F(D) curves. Each curve is reconstructed from raw vDeflection/measuredHeight, binned every 5 nm, and fitted independently over 20–200 nm. No map-level median or four-point |F| time curve is used as a fit input.

Force reconstruction uses global liquid InvOLS `68.289604 nm/V` and independently calibrated `k=0.311340173 N/m`.

## Requested colloid-colloid model

This is a deliberate geometry-mismatch diagnostic because the experiment is silica sphere–plane, while this fit assumes two identical silica spheres with `R1=R2=4.546849 µm` and `R_eff=R/2=2.273424 µm`.

The equal-potential linearized HHF and two-sphere Derjaguin terms are

`F_EDL = 4π R_eff ε κ ψ² / [exp(κD)+1]`,

`F_vdW = -A_H R_eff/(6D²)`,

`F_hyd = 6πη R_eff² v_gap/D`,

with a fitted additive baseline. Primary fits include the fixed pure-water no-slip hydrodynamic term; no-hydrodynamics and far-constant branches are saved as sensitivity variants. Free parameters are apparent `lambda_D`, equal-surface `|psi|`, and baseline.

## Primary per-map distributions

Values are medians [IQR] over 64 independent curve fits.

| map | time (min) | apparent lambda_D (nm) | apparent |psi| (mV) | baseline (pN) | RMSE (pN) | R2 | interior fits | <=25 mV |
|---:|---:|:---|:---|:---|:---|:---|---:|---:|
| 2 | 0.00 | 23.65 [20.06, 26.71] | 38.5 [31.1, 47.4] | -39.0 [-88.0, +60.9] | 11.3 [8.7, 12.9] | 0.987 [0.980, 0.993] | 64/64 | 1.6% |
| 3 | 27.16 | 17.40 [14.06, 21.61] | 36.1 [31.2, 45.7] | -22.4 [-85.8, +69.4] | 14.9 [11.7, 18.4] | 0.964 [0.937, 0.989] | 64/64 | 6.2% |
| 4 | 53.33 | 13.67 [11.62, 18.98] | 31.1 [24.9, 39.0] | -27.3 [-58.0, +42.1] | 13.9 [11.6, 17.2] | 0.937 [0.859, 0.980] | 64/64 | 26.6% |
| 5 | 78.66 | 12.10 [11.04, 14.00] | 23.0 [17.6, 31.3] | -10.7 [-47.2, +37.9] | 7.8 [6.9, 9.1] | 0.930 [0.836, 0.981] | 64/64 | 54.7% |

## Interpretation boundary

- These are per-curve apparent parameters under the user-requested wrong geometry, not equilibrium sphere–plane zeta potentials or bulk Debye lengths.
- Equal-sphere conservative-force amplitudes are half the actual sphere–plane Derjaguin amplitudes at the same physical radius. The fitted potential therefore compensates for the geometry choice.
- The HHF expression is linearized PB. Fits above the 25 mV reference ceiling are retained numerically but flagged as outside the reference implementation's operational potential range.
- Adjacent 5 nm bins from one curve are correlated, and hydrodynamics/history/baseline remain model-conditioned. Across-curve IQR is more relevant than local optimizer precision.

## Outputs

- `per_curve_colloid_colloid_fits.csv`: all fit parameters and diagnostics.
- `per_curve_colloid_colloid_predictions.csv`: observed and fitted values at every included distance bin.
- `map_colloid_colloid_fit_summary.csv`: map-level medians and IQRs for all variants.
- `figures/all_256_colloid_colloid_curve_fits.png` and `figures/apparent_colloid_colloid_parameters_vs_time.png`.
- `provenance.json` and `artifact_manifest.sha256`.

# 20-08-26 silica sphere-plane surface-force fit

## Result status

This analysis decodes raw JPK `vDeflection` in volts and raw scanner readback in metres. The JPK header InvOLS and force channels were retained only for comparison and were never used to calculate deflection or force. Force was rebuilt with the cantilever-1 calibration `k = 0.2969899087 N/m` at 25.6 C.

The operational primary result uses the maps acquired over the same area: at each separation, the measured 1/2/4 um/s map medians are extrapolated to zero approach speed. The 25 independent curves form a separate formal cohort and a cross-check, rather than being allowed to outvote a map merely because they are stored as five files. The zero-speed map intercept is fitted to a nonlinear, equal-constant-potential, symmetric 1:1 Poisson-Boltzmann planar solution and converted to sphere-plane force with the Derjaguin approximation. This avoids using the reference repository's equal-sphere prefactor.

Important experimental consistency finding: the three sequential maps show force decreasing with approach speed over 20-100 nm, opposite to the positive no-slip lubrication slope. Therefore the numerical PB fits below are reproducible operational fits, but the speed extrapolation cannot be claimed as a validated hydrodynamic correction or a unique equilibrium measurement. Sequential surface-state/contact-alignment changes remain confounded with speed.

| glycerol wt% | Debye length (nm) | common-potential magnitude (mV) | silica-assigned sign (mV) | R2 | validity |
|---:|---:|---:|---:|---:|:---|
| 0 | 19.2 [18.2, 20.2] | 76 [75.7, 76.3] | -76 | 0.9849 | numerically valid; map speed inconsistency |
| 10 | 15.9 [14.9, 17] | 12 [11.7, 12.2] | -12 | 0.9879 | numerically valid; map speed inconsistency |
| 20 | 22.1 [21.5, 22.7] | 74.8 [74, 75.5] | -74.8 | 0.9949 | numerically valid; map speed inconsistency |
| 30 | 22 [20.1, 24.2] | 72.7 [70.2, 75.3] | -72.7 | 0.9683 | numerically valid; map speed inconsistency |
| 40 | 25.7 [25.2, 26.2] | 81.4 [80.7, 82.1] | -81.4 | 0.9989 | numerically valid; map speed inconsistency |

The confidence intervals above are local Jacobian intervals conditional on the chosen contact alignment, common-potential boundary condition, dielectric constants, radius and Hamaker constant. They are much narrower than the experimental systematic spread and must not be quoted alone.

| wt% | Debye-length systematic range (nm) | potential-magnitude systematic range (mV) |
|---:|---:|---:|
| 0 | 18.5-19.6 | 76-76.2 |
| 10 | 5.55-16.5 | 11.8-15.2 |
| 20 | 14.4-22.3 | 70.9-75.2 |
| 30 | 11.2-23.7 | 72.6-75.8 |
| 40 | 14.6-26.5 | 66-82.6 |

These ranges span nonlinear map fit-window, all-source and independent-force-only variants where available. They are sensitivity ranges, not probability intervals. The 10 wt% Debye length is especially weak because many traces lack a resolved hard-contact interval and the force signal is small.

Force between identical surfaces is even in the potential, so AFM normal-force data determine `|zeta|`, not its sign. The negative sign in the table is assigned from silica surface chemistry, not inferred from the force curve. More precisely, the fitted quantity is the constant-potential PB boundary potential; calling it zeta additionally assumes that this potential represents the electrokinetic slipping-plane potential.

## Model and derivation

For dimensionless potential `u=e psi/(k_B T)` and dimensionless plate spacing `H=kappa D`, the equal-potential planar boundary-value problem is

`u'' = sinh(u)`, `u(+/-H/2)=u_s`, and `u'(0)=0`.

Its first integral gives

`H/2 = integral_(u_m)^(u_s) du / sqrt(2[cosh(u)-cosh(u_m)])`.

At the symmetry plane the Maxwell-field term is zero, so `Pi = eps (k_B T/e)^2 kappa^2 [cosh(u_m)-1]`. The sphere-plane Derjaguin force used here is therefore

`F_EDL(D)=2 pi R eps (k_B T/e)^2 kappa integral_(kappa D)^infinity [cosh(u_m(H))-1] dH`.

The fixed van der Waals term is `F_vdW=-A_H R/(6D^2)` with `A_H=2.4e-21 J`. The fitted common potential and inverse screening length are shared by the silica sphere and silica plane within each concentration.

In the low-potential limit this implementation reduces to the sphere-plane same-potential HHF expression

`F_EDL = 4 pi R eps kappa zeta^2 / [exp(kappa D)+1]`.

For two equal spheres the Derjaguin radius would be `R/2`; therefore copying the reference paper's equal-sphere formula would understate the present sphere-plane force by a factor of two for the same physical sphere radius.

Primary sources: [Hogg-Healy-Fuerstenau, linear Debye-Huckel sphere interactions](https://pubs.rsc.org/en/content/articlehtml/1966/tf/tf9666201638); [Stankovich and Carnie, nonlinear PB sphere/plate and Derjaguin accuracy](https://pubs.acs.org/doi/10.1021/la950384k); [Polat and Polat, nonlinear PB for arbitrary same-sign plate potentials](https://doi.org/10.1016/j.jcis.2009.09.008); [Liu and Li, AFM zeta workflow used by the reference repository](https://doi.org/10.1016/j.jcis.2020.05.061).

## Sensitivity reconstruction

Each usable approach trace was searched near its terminal end for a contiguous rigid-contact interval. A robust line `V=a+b h` gives `InvOLS=1/|b|`. First, each source supplies a robust contact-slope anchor. Then all consistent source anchors at the same glycerol concentration are combined with equal source weighting into one concentration-wide sensitivity. This keeps the force calibration independent of approach speed and uses all available raw hard-contact data; the embedded JPK sensitivity is never substituted.

| wt% | sources | raw curves | source anchors | retained anchors | common sensitivity used (nm/V) |
|---:|---:|---:|---:|---:|---:|
| 0 | 8 | 773 | 8 | 8 | 66.609 |
| 10 | 8 | 773 | 6 | 6 | 64.177 |
| 20 | 8 | 773 | 8 | 7 | 60.718 |
| 30 | 11 | 1011 | 11 | 11 | 58.673 |
| 40 | 8 | 773 | 8 | 8 | 57.903 |

A curve lacking a resolvable hard-contact interval cannot independently establish its contact plane. For those curves, the terminal `h+delta` value is used as the mechanical-contact estimate and the curve is flagged as endpoint-aligned. This limitation is particularly important for 10 wt% and is visible in the source-subset/window sensitivity results.

## Far-field linear drift removal and spatial QC

Before contact finding or force conversion, every approach trace uses its initial 20% of samples (200 of 1000 for the present data; minimum 80) as the far field. A four-pass robust line `V_ff(h)=a h+b` is fitted in raw `vDeflection` versus raw `measuredHeight`, with 3.5-MAD residual rejection. The full fitted line is subtracted from the entire trace: `V_corrected=V-V_ff`. This is the linear drift-removal step used by all sensitivity, force, speed-extrapolation and PB results; the later 700-880 nm source correction removes only a residual constant offset.

The slope below is reported after applying the concentration-wide calculated sensitivity and cantilever-1 spring constant, as `dF/d(measuredHeight)` in pN per 100 nm scanner motion. The sign is therefore tied to the scanner-height coordinate, not acquisition time. `resolved` means the fitted change across the far window exceeds twice the post-fit residual MAD; it is diagnostic only, because the declared linear operator is applied to every valid trace.

Every map declares a 16x16, back-and-forth raster. Odd acquisition rows are reversed back into physical x before plotting. Dark-gray pixels in partial 30 wt% maps are unacquired records, not interpolated values.

| wt% | start time | speed (um/s) | pixels | median slope +/- MAD (pN/100 nm) | resolved | corr. acquisition order | spatial-plane R2 |
|---:|:---|---:|---:|---:|---:|---:|---:|
| 0 | 15:01:47 | 2 | 256/256 | -3.21 +/- 4.88 | 8.2% | -0.259 | 0.086 |
| 0 | 15:06:26 | 1 | 256/256 | -4.02 +/- 4.12 | 15.6% | -0.277 | 0.080 |
| 0 | 15:15:57 | 4 | 256/256 | -7.01 +/- 4.41 | 9.4% | -0.027 | 0.012 |
| 10 | 15:29:51 | 2 | 256/256 | -0.958 +/- 6.95 | 9.0% | -0.082 | 0.008 |
| 10 | 15:35:17 | 1 | 256/256 | 17.5 +/- 13 | 70.3% | 0.432 | 0.251 |
| 10 | 15:44:18 | 4 | 256/256 | 10.8 +/- 11.6 | 30.1% | 0.378 | 0.182 |
| 20 | 15:54:06 | 2 | 256/256 | 12.1 +/- 7.29 | 50.0% | 0.250 | 0.305 |
| 20 | 15:59:55 | 1 | 256/256 | 33.6 +/- 12.3 | 94.1% | -0.389 | 0.171 |
| 20 | 16:09:19 | 4 | 256/256 | 4.5 +/- 12.1 | 12.1% | 0.013 | 0.005 |
| 30 | 16:20:49 | 2 | 149/256 | 0.656 +/- 4.35 | 1.3% | 0.030 | 0.001 |
| 30 | 16:24:44 | 2 | 111/256 | 4.34 +/- 7.64 | 10.8% | 0.421 | 0.167 |
| 30 | 16:27:34 | 2 | 208/256 | 8.49 +/- 12.7 | 34.6% | 0.381 | 0.131 |
| 30 | 16:37:53 | 1 | 101/256 | 23.2 +/- 17.8 | 80.2% | -0.046 | 0.001 |
| 30 | 16:42:02 | 1 | 181/256 | 34.4 +/- 18.4 | 83.4% | 0.044 | 0.003 |
| 30 | 16:53:04 | 4 | 256/256 | 20.1 +/- 12.1 | 39.5% | -0.041 | 0.006 |
| 40 | 17:03:03 | 2 | 256/256 | 1.8 +/- 5.2 | 3.5% | 0.005 | 0.001 |
| 40 | 17:07:48 | 1 | 256/256 | 19.7 +/- 16.1 | 67.2% | 0.297 | 0.108 |
| 40 | 17:18:39 | 4 | 256/256 | 8.91 +/- 13.9 | 18.0% | -0.043 | 0.022 |

Per-pixel values, fit quality, residual noise, drift/noise ratio, raster coordinates and plane-gradient summaries are retained in `far_field_drift_by_curve.csv` and `far_field_drift_by_map.csv`; no visual clipping changes the saved numerical values.

## Speed dependence and liquid properties

At each 5 nm bin, the same-area map medians at the measured approximately 1, 2 and 4 um/s speeds are regressed linearly to zero speed. This is an operational intercept, not a validated lubrication subtraction. The measured slope is checked against the positive no-slip prediction `6 pi eta R^2/D` using the Cheng glycerol-water viscosity correlation.

| wt% | epsilon_r | viscosity (mPa s) | map slope / lubrication slope | map positive-slope fraction | independent-force slope / lubrication slope |
|---:|---:|---:|---:|---:|---:|
| 0 | 78.500 | 0.8806 | -3.49 | 0.000 | 5.94 |
| 10 | 76.000 | 1.1292 | -1.71 | 0.000 | -0.738 |
| 20 | 73.900 | 1.5013 | -3.17 | 0.059 | 11.7 |
| 30 | 71.400 | 2.0886 | -3.58 | 0.000 | 8.09 |
| 40 | 68.900 | 3.0796 | -2.43 | 0.000 | 17.6 |

All map medians fail the expected positive lubrication direction. Some independent-force cohorts instead have large positive slopes, so the two acquisition modes are not mutually interchangeable replicates. Viscosity source: [Cheng, Formula for the Viscosity of a Glycerol-Water Mixture](https://doi.org/10.1021/ie071349z). The dielectric constants (78.5, 76.0, 73.9, 71.4 and 68.9) are rounded from the 25 C, 0.57 MHz measurements tabulated at 0.00, 9.88, 20.33, 30.19 and 39.67 wt% in [Physical Properties of Glycerine and Its Solutions, Table 27](https://www.cleaninginstitute.org/sites/default/files/research-pdfs/Physical_properties_of_glycerine_and_its_solutions.pdf), and agree with the reference repository. Modern primary measurements over composition and 10-50 C are reported by [Behrends et al.](https://pubmed.ncbi.nlm.nih.gov/16626219/). The PB temperature itself is the measured 25.6 C; the closest tabulated permittivity temperature is 25 C.

## Numerical and provenance checks

- Parsed 43 raw sources and 4103 usable approach curves; partial/short terminal records skipped: 5.
- Far-field linear-drift synthetic slope relative error: `5.995e-15`; maximum residual after subtraction: `4.163e-17 V`.
- Nonlinear PB to analytic linear-PB limit maximum relative error: `1.652e-03`.
- Nonlinear PB far-field asymptote maximum relative error: `1.054e-03`.
- Cheng endpoint check at 25.6 C: water `0.880645 mPa s`, glycerol `860.334 mPa s`.
- Every input ZIP passed CRC during decoding; SHA-256 values are in `input_manifest.csv`.
- Optimizer termination alone is not accepted: fit validity also checks Jacobian rank/condition, parameter bounds, R2 and signal-to-noise.

## Files

- `sensitivity_by_source.csv`: raw source anchors, retained-anchor flags, common calculated sensitivity and ignored embedded calibration.
- `far_field_drift_by_curve.csv`: fitted slope, R2, residual noise and drift/noise for every approach trace.
- `far_field_drift_by_map.csv`: 16x16 map-level spatial-gradient and acquisition-order diagnostics.
- `source_binned_force_curves.csv`: reconstructed source-median force curves by actual approach speed.
- `zero_speed_equilibrium.csv`: three-speed intercepts and speed slopes.
- `fit_results.csv`: nonlinear primary fit, fit-window/source-subset checks and linear-HHF diagnostic.
- `mixture_properties.csv`, `input_manifest.csv`, `provenance.json`, figures and SHA-256 manifest.

# AFM2026

Reproducible analysis code, derived results, quality-control figures and experimental planning for silica colloidal-probe AFM measurements in water–glycerol mixtures at 25.6 °C.

## Repository contents

- `analysis/calibrate_cantilevers.py`: raw thermal-noise and hard-contact calibration workflow.
- `analysis/fit_glycerol_surface_forces.py`: raw reconstruction and silica sphere–plane surface-force fitting.
- `analysis/analyze_velocity_systematics.py`: contact, snap-in, pull-off, raster and acquisition-history diagnostics.
- `analysis/analyze_velocity_joint_fit.py`: same-pixel velocity comparisons, hydrodynamic shape tests and joint PB fits.
- `analysis/plot_zero_wt_all_fd_curves.py`: all measured 0 wt% approach curves, split into three speed-highlight figures with pointwise median and interquartile range.
- `analysis/analyze_zero_wt_distribution_separation.py`: spatial-block bootstrap, simultaneous confidence bands and exact functional sign-flip tests for the three 0 wt% maps.
- `analysis/analyze_zero_wt_classical_distribution_tests.py`: paired t/Wilcoxon, marginal Welch/Mann-Whitney/KS and three-group repeated-measures tests for the same maps.
- `analysis/*_results/`: derived CSV data, figures, reports, provenance and SHA-256 manifests.
- `analysis/FOLLOWUP_8X8_PALINDROME_FORCE_MAPPING_PROTOCOL.md`: detailed follow-up experiment using balanced 8×8 palindromic force mapping, including the 2026-08-24 extension to 0–99.5 wt% glycerol and equilibrium-force recovery over 20–200 nm.

## Current numerical snapshot

The three cantilevers were calibrated at 25.6 °C. Their recommended `(InvOLS, k, f0, Q)` values are respectively `(85.30 nm/V, 0.2970 N/m, 16.1101 kHz, 68.91)`, `(65.72 nm/V, 0.2736 N/m, 15.0774 kHz, 64.75)` and `(91.56 nm/V, 0.2384 N/m, 15.0675 kHz, 67.62)`. The cantilever-2 spectrum whose header selected 97.67 kHz was re-fitted from its measured PSD; the fundamental is 15.0777 kHz and the 97.67 kHz peak is retained only as a higher-mode observation.

The glycerol experiment uses cantilever 1 with `k = 0.2969899087 N/m`. Embedded JPK sensitivity and force are never used: raw hard-contact data give concentration-wide sensitivities of `66.609`, `64.177`, `60.718`, `58.673` and `57.903 nm/V` at 0, 10, 20, 30 and 40 wt%. Ten wt% is excluded from the primary interpretation because its terminal load and hard-contact validity define a different measurement regime.

For 0 wt%, the three map medians at 50 nm are `514.4`, `573.0` and `447.2 pN` at 1, 2 and 4 µm/s; at 100 nm they are `35.9`, `57.6` and `13.2 pN`. Thus the measured ordering is `F(2) > F(1) > F(4)`, not a monotonic increase with speed. Spatial-block tests show that the three measured map distributions are separated, but the acquisition order was `2 -> 1 -> 4 µm/s` and there is only one map per speed. The separation is therefore conditional evidence for a map/time/history difference, not a replicated causal velocity effect or a validated hydrodynamic subtraction.

The joint equal-silica sphere-plane PB fits, excluding 10 wt%, give model-conditioned primary values `lambda_D = 18.9, 18.1, 16.3, 21.5 nm` and `|zeta| = 75.2, 78.1, 79.3, 85.8 mV` at 0, 20, 30 and 40 wt%. Across the tested model/window/contact perturbations, the combined ranges expand to `lambda_D = 15.5-24.1 nm` and `|zeta| = 68.3-99.3 mV`. These are exploratory systematic ranges because the fitted hydrodynamic amplitude frequently has the wrong sign and remains confounded with map history.

The current scientific interpretation and its limitations are documented in:

- [`analysis/velocity_joint_fit_results/REPORT.md`](analysis/velocity_joint_fit_results/REPORT.md)
- [`analysis/velocity_joint_fit_results/ZERO_WT_VELOCITY_DISTRIBUTION_SUMMARY.md`](analysis/velocity_joint_fit_results/ZERO_WT_VELOCITY_DISTRIBUTION_SUMMARY.md)
- [`analysis/velocity_joint_fit_results/ZERO_WT_DISTRIBUTION_SEPARATION.md`](analysis/velocity_joint_fit_results/ZERO_WT_DISTRIBUTION_SEPARATION.md)
- [`analysis/velocity_joint_fit_results/ZERO_WT_CLASSICAL_DISTRIBUTION_TESTS.md`](analysis/velocity_joint_fit_results/ZERO_WT_CLASSICAL_DISTRIBUTION_TESTS.md)
- [`analysis/velocity_systematics_results/REPORT.md`](analysis/velocity_systematics_results/REPORT.md)
- [`analysis/surface_force_results/REPORT.md`](analysis/surface_force_results/REPORT.md)
- [`analysis/results/REPORT.md`](analysis/results/REPORT.md)

The future concentration range is 0–99.5 wt% glycerol. At 25.6 °C the Cheng viscosity estimate rises from `0.8806 mPa s` in water to `774.9 mPa s` at 99.5 wt%, so a fixed approach speed of 0.1 or 0.2 µm/s cannot make hydrodynamic drainage negligible across the range. The protocol records concentration-dependent `eta U` speed design, finite-distance zero-force semantics and the requirement that the reported 20–200 nm equilibrium force be obtained from a map/block/session-supported `U -> 0` intercept rather than from one finite-speed curve.

## Raw-data policy

Original JPK force maps, individual force curves, thermal-noise files and acquisition photographs are intentionally not stored in this Git repository. The `.gitignore` excludes the local `20-08-26/` and `calibration/` trees as well as the raw JPK/TND file extensions globally.

The derived provenance and input manifests retain source identities and hashes so that results can be checked against an authorized raw-data archive without placing the raw measurements in Git.

## Reproduction

With the excluded raw folders restored at the paths expected by the scripts, the analyses are run from the repository root:

```bash
python analysis/calibrate_cantilevers.py
python analysis/fit_glycerol_surface_forces.py
python analysis/analyze_velocity_systematics.py
python analysis/analyze_velocity_joint_fit.py
python analysis/plot_zero_wt_all_fd_curves.py
python analysis/analyze_zero_wt_distribution_separation.py
python analysis/analyze_zero_wt_classical_distribution_tests.py
```

Each result directory contains an `artifact_manifest.sha256` file for checking the primary committed derived artifacts. The supplemental 0 wt% analysis has its own `zero_wt_distribution_analysis_manifest.sha256` so it can be verified independently.

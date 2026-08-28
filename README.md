# AFM2026

Reproducible analysis code, derived results, quality-control figures and experimental planning for silica colloidal-probe AFM measurements in water–glycerol mixtures at 25.6 °C.

## Repository contents

- `analysis/calibrate_cantilevers.py`: raw thermal-noise and hard-contact calibration workflow.
- `analysis/fit_glycerol_surface_forces.py`: raw reconstruction and silica sphere–plane surface-force fitting.
- `analysis/analyze_velocity_systematics.py`: contact, snap-in, pull-off, raster and acquisition-history diagnostics.
- `analysis/analyze_velocity_joint_fit.py`: same-pixel velocity comparisons, hydrodynamic shape tests and joint PB fits.
- `analysis/analyze_27_08_26_palindrome_pilot.py`: D4–D6 air calibration plus raw D4 pure-water 8×8 palindrome pilot, including same-pixel early/late decomposition, contact/retract history QC and baseline sensitivity.
- `analysis/analyze_27_08_26_full_palindrome.py`: full 29-map reconstruction of the overlapping map1–3 liquid-refresh test and map3–5 no-refresh test, including time-aware velocity models, contact/retract-state diagnostics and explicit handling of the missing final 0.1 µm/s map.
- `analysis/analyze_27_08_26_palindrome_distribution_tests.py`: same-pixel palindrome-mean force distributions with 64-pixel, sixteen 2×2 spatial-tile and n=3 block-level paired t-tests, including clock-center and baseline sensitivity.
- `analysis/plot_27_08_26_fixed_pixel_chronology.py`: chronological force-distance overlay and 20/50/100/200 nm force slices for one fixed physical 8×8 map position.
- `analysis/plot_zero_wt_all_fd_curves.py`: all measured 0 wt% approach curves, split into three speed-highlight figures with pointwise median and interquartile range.
- `analysis/analyze_zero_wt_distribution_separation.py`: spatial-block bootstrap, simultaneous confidence bands and exact functional sign-flip tests for the three 0 wt% maps.
- `analysis/analyze_zero_wt_classical_distribution_tests.py`: paired t/Wilcoxon, marginal Welch/Mann-Whitney/KS and three-group repeated-measures tests for the same maps.
- `analysis/*_results/`: derived CSV data, figures, reports, provenance and SHA-256 manifests.
- `analysis/FOLLOWUP_8X8_PALINDROME_FORCE_MAPPING_PROTOCOL.md`: detailed follow-up experiment using balanced 8×8 palindromic force mapping, including the 2026-08-24 extension to 0–99.5 wt% glycerol and equilibrium-force recovery over 20–200 nm.

## Current numerical snapshot

The three cantilevers were calibrated at 25.6 °C. Their recommended `(InvOLS, k, f0, Q)` values are respectively `(85.30 nm/V, 0.2970 N/m, 16.1101 kHz, 68.91)`, `(65.72 nm/V, 0.2736 N/m, 15.0774 kHz, 64.75)` and `(91.56 nm/V, 0.2384 N/m, 15.0675 kHz, 67.62)`. The cantilever-2 spectrum whose header selected 97.67 kHz was re-fitted from its measured PSD; the fundamental is 15.0777 kHz and the 97.67 kHz peak is retained only as a higher-mode observation.

The glycerol experiment uses cantilever 1 with `k = 0.2969899087 N/m`. Embedded JPK sensitivity and force are never used: raw hard-contact data give concentration-wide sensitivities of `66.609`, `64.177`, `60.718`, `58.673` and `57.903 nm/V` at 0, 10, 20, 30 and 40 wt%. Ten wt% is excluded from the primary interpretation because its terminal load and hard-contact validity define a different measurement regime.

For 0 wt%, the three line-corrected map medians at 50 nm are `514.0`, `572.6` and `447.2 pN` at 1, 2 and 4 µm/s; at 100 nm they are `35.7`, `57.5` and `13.2 pN`. Thus the measured ordering is `F(2) > F(1) > F(4)`, not a monotonic increase with speed. Spatial-block tests show that the three measured map distributions are separated, but the acquisition order was `2 -> 1 -> 4 µm/s` and there is only one map per speed. The separation is therefore conditional evidence for a map/time/history difference, not a replicated causal velocity effect or a validated hydrodynamic subtraction.

The velocity joint-fit report now tabulates measured map-force slices at 25, 50, 100, 150 and 200 nm for 0, 20, 30 and 40 wt%. At 50 nm, the line-corrected 1/2/4 µm/s values are respectively `514.0/572.6/447.2`, `551.6/610.2/463.8`, `388.4/579.4/345.3` and `639.7/699.8/509.7 pN`; 30 wt% uses only the comparable high-load stratum. The corresponding constant-referenced values are `547.8/595.6/501.6`, `285.5/507.5/429.8`, `146.4/511.7/185.8` and `439.6/688.6/444.9 pN`. These branch differences show that far-field linear subtraction materially changes the reported force and cannot be treated as an innocuous definition of absolute zero.

The joint equal-silica sphere-plane PB fits, excluding 10 wt%, give model-conditioned primary values `lambda_D = 18.9, 18.1, 16.3, 21.5 nm` and `|zeta| = 75.2, 78.1, 79.3, 85.8 mV` at 0, 20, 30 and 40 wt%. Across the tested model/window/contact perturbations, the combined ranges expand to `lambda_D = 15.5-24.1 nm` and `|zeta| = 68.3-99.3 mV`. These are exploratory systematic ranges because the fitted hydrodynamic amplitude frequently has the wrong sign and remains confounded with map history.

The current scientific interpretation and its limitations are documented in:

- [`analysis/palindrome_27_08_26_full_results/REPORT.md`](analysis/palindrome_27_08_26_full_results/REPORT.md)
- [`analysis/palindrome_27_08_26_distribution_results/REPORT.md`](analysis/palindrome_27_08_26_distribution_results/REPORT.md)
- [`analysis/palindrome_27_08_26_pilot_results/REPORT.md`](analysis/palindrome_27_08_26_pilot_results/REPORT.md)
- [`analysis/velocity_joint_fit_results/REPORT.md`](analysis/velocity_joint_fit_results/REPORT.md)
- [`analysis/velocity_joint_fit_results/ZERO_WT_VELOCITY_DISTRIBUTION_SUMMARY.md`](analysis/velocity_joint_fit_results/ZERO_WT_VELOCITY_DISTRIBUTION_SUMMARY.md)
- [`analysis/velocity_joint_fit_results/ZERO_WT_DISTRIBUTION_SEPARATION.md`](analysis/velocity_joint_fit_results/ZERO_WT_DISTRIBUTION_SEPARATION.md)
- [`analysis/velocity_joint_fit_results/ZERO_WT_CLASSICAL_DISTRIBUTION_TESTS.md`](analysis/velocity_joint_fit_results/ZERO_WT_CLASSICAL_DISTRIBUTION_TESTS.md)
- [`analysis/velocity_systematics_results/REPORT.md`](analysis/velocity_systematics_results/REPORT.md)
- [`analysis/surface_force_results/REPORT.md`](analysis/surface_force_results/REPORT.md)
- [`analysis/results/REPORT.md`](analysis/results/REPORT.md)

The future concentration range is 0–99.5 wt% glycerol. At 25.6 °C the Cheng viscosity estimate rises from `0.8806 mPa s` in water to `774.9 mPa s` at 99.5 wt%, so a fixed approach speed of 0.1 or 0.2 µm/s cannot make hydrodynamic drainage negligible across the range. The protocol records concentration-dependent `eta U` speed design, finite-distance zero-force semantics and the requirement that the reported 20–200 nm equilibrium force be obtained from a map/block/session-supported `U -> 0` intercept rather than from one finite-speed curve.

The complete 27-08-26 follow-up contains 29 pure-water maps. Map1–3 is a liquid-refresh-affected palindrome test; map3–5 is the primary no-refresh test, with only the final 0.1 µm/s map missing from map5. Re-fitting the air calibrations at 25.6 °C gives D4/D5/D6 `(InvOLS, k)` values of `(72.76 nm/V, 0.164829 N/m)`, `(79.90 nm/V, 0.311340 N/m)` and `(70.90 nm/V, 0.181704 N/m)`. One of the five D4 thermal-noise files is an exact byte duplicate, so D4 uses four unique spectra. The force maps themselves give a separate global D4 water-contact InvOLS of `51.822 nm/V`; all map forces are rebuilt from this value and the calibrated D4 spring constant rather than the embedded JPK sensitivity or force.

In the no-refresh map3–5 chronology, the far-linear map-median force falls from `2578.5` to `495.9 pN` at 20 nm, `647.6` to `16.2 pN` at 50 nm and `77.1` to `8.2 pN` at 100 nm, while the local contact InvOLS changes by only about `-0.73%` and the terminal contact load remains near `21.8 nN`. All eight complete same-speed pairs have a negative 50 nm history residual. Endpoint velocity slopes change sign between blocks, and time-aware fits depend strongly on time model and far-field baseline; every tested HC3 95% interval includes zero. The present data therefore demonstrate strong acquisition-history relaxation but do not identify a hydrodynamic slope or support a defensible `U -> 0` force correction. Retract pull-off is also systematically less negative after the nominally faster approach, showing a speed/order-correlated contact-history effect that cannot yet be separated causally from time.

Palindromic early/late averaging does reveal conditional within-map velocity-distribution shifts at selected distances, but it does not establish replicated velocity separation. In primary Test B, the n=3 block-level paired t-test gives raw p-values `0.757/0.338/0.405/0.440` at 20/50/100/200 nm for the 0.2-minus-0.05 µm/s contrast, and no block-level result survives Holm correction over 20–200 nm. Pixel p-values are therefore retained only as within-area evidence; sixteen 2×2 spatial tiles and the three map blocks are reported separately to avoid treating 64 correlated pixels as independent experiments.

## Raw-data policy

Original JPK force maps, individual force curves, thermal-noise files and acquisition photographs are intentionally not stored in this Git repository. The `.gitignore` excludes the local `20-08-26/`, `27-08-26/` and `calibration/` trees as well as the raw JPK/TND file extensions globally.

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
python analysis/analyze_27_08_26_palindrome_pilot.py
python analysis/plot_27_08_26_fixed_pixel_chronology.py
python analysis/analyze_27_08_26_full_palindrome.py
python analysis/analyze_27_08_26_palindrome_distribution_tests.py
```

Each result directory contains an `artifact_manifest.sha256` file for checking the primary committed derived artifacts. The supplemental 0 wt% analysis has its own `zero_wt_distribution_analysis_manifest.sha256` so it can be verified independently.

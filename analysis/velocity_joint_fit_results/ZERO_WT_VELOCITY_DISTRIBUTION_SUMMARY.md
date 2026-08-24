# 0 wt% velocity-dependent F-D distributions: integrated summary

## Main finding

The three measured 0 wt% map distributions are statistically separated over much of `25-250 nm`, but they do **not** follow a monotonic hydrodynamic velocity law. In the force-bearing `50-100 nm` region the observed ordering is

`F(2 um/s) > F(1 um/s) > F(4 um/s)`.

The maps were acquired chronologically as `2 -> 1 -> 4 um/s`, with only one map at each speed. The force level therefore also decreases with acquisition history. The data establish separation among these three maps; they do not identify speed as its cause independently of time, surface conditioning, optical/contact response or residual baseline history.

## Data and force definition

- Environment: 0 wt% glycerol (water), 25.6 °C; silica sphere against a nominally identical silica plane.
- Cantilever: cantilever 1, calibrated spring constant `k = 0.2969899087 N/m`.
- In-experiment 0 wt% sensitivity: `66.608814 nm/V`, calculated from all eight accepted raw hard-contact source anchors. This differs from the separate calibration-session InvOLS `85.2963 nm/V`; force reconstruction deliberately uses the sensitivity measured from the experiment itself and the independently calibrated `k`.
- Reconstruction: raw `vDeflection` is decoded in volts. A robust line fitted to the initial 20% of each approach trace is subtracted; then `delta = S V_corrected`, `F = k delta`, and `D = h + delta - h_contact`. Embedded JPK sensitivity and force channels are not used.
- Zero-force convention: zero is the extrapolated per-curve far-field line, not an independently measured force at infinite separation. Any slowly varying physical contribution present in the reference window can therefore be partly absorbed into the constant/linear baseline.

The map inventory is:

| acquisition order | start time | nominal speed (um/s) | map pixels reconstructed | pixels common to all three maps |
|---:|:---|---:|---:|---:|
| 1 | 15:01:47 | 2 | 256 | 255 |
| 2 | 15:06:26 | 1 | 256 | 255 |
| 3 | 15:15:57 | 4 | 256, including one endpoint-aligned curve | 255 |

The all-curve figures additionally include independent force curves: `259` curves at 1 um/s (`256 map + 3 independent`), `258` at 2 um/s (`256 + 2`) and `256` at 4 um/s (`256 + 0`). Statistical comparisons use only the `255` same-position map pixels common to all speeds; independent curves are not mixed with map pixels as if they were exchangeable replicates.

## Directly measured force levels

The table gives map-pixel median `[25th, 75th percentile]` in pN on the line-corrected branch. It is the force at the stated sphere-plane separation, not one scalar assigned to an entire force curve.

| D (nm) | 1 um/s (pN) | 2 um/s (pN) | 4 um/s (pN) |
|---:|:---|:---|:---|
| 25 | 1902.5 [1838.6, 1957.7] | 1931.3 [1861.9, 1977.3] | 1819.7 [1739.4, 1881.1] |
| 50 | 514.4 [479.0, 544.9] | 573.0 [533.0, 597.0] | 447.2 [404.7, 476.0] |
| 100 | 35.9 [14.8, 57.4] | 57.6 [27.6, 79.7] | 13.2 [-8.0, 37.3] |
| 150 | -12.0 [-28.7, 9.1] | -3.5 [-27.8, 14.9] | -21.0 [-44.2, -1.3] |
| 200 | -15.1 [-30.7, 8.3] | -8.0 [-31.0, 9.7] | -25.5 [-45.6, -6.8] |
| 250 | -9.7 [-25.1, 8.3] | -5.8 [-25.7, 13.6] | -18.1 [-39.1, -1.8] |

The small negative far-range values are relative to the fitted far-field line. They are not an attractive-force absolute calibration and should not be interpreted without the baseline convention above.

## Distribution-separation evidence

The primary traditional comparison uses paired 4x4 spatial-block medians (`16` matched blocks). Pairwise p-values are Holm-corrected over `3 pairs x 46 distances`. The robust companion analysis uses `10,000` spatial-block bootstrap replicates and one studentized max-T band simultaneously covering all pairs and distances.

At the two most informative separations:

| D (nm) | comparison | block mean difference (pN) | 95% t-CI (pN) | Holm paired-t p | paired median difference (pN) | simultaneous 95% CI (pN) | P(F_fast > F_slow) |
|---:|:---|---:|:---|---:|---:|:---|---:|
| 50 | 2-1 | +56.66 | [45.94, 67.38] | 1.21e-6 | +52.95 | [27.54, 78.37] | 0.831 |
| 50 | 4-2 | -128.60 | [-142.41, -114.78] | 4.72e-10 | -128.53 | [-149.01, -108.05] | 0.043 |
| 50 | 4-1 | -71.94 | [-84.24, -59.63] | 3.12e-7 | -77.94 | [-95.93, -59.94] | 0.102 |
| 100 | 2-1 | +17.82 | [10.14, 25.51] | 0.0123 | +16.90 | [8.04, 25.76] | 0.659 |
| 100 | 4-2 | -39.23 | [-47.93, -30.52] | 9.29e-6 | -44.71 | [-61.72, -27.70] | 0.161 |
| 100 | 4-1 | -21.40 | [-28.75, -14.05] | 0.00155 | -23.20 | [-32.73, -13.68] | 0.286 |

Across the full distance window:

- The 4x4-block repeated-measures ANOVA rejects equality of the three maps at all `46/46` bins from `25` to `250 nm`; the largest Holm-adjusted p-value is `0.00591`. At 50 nm, `F(2,30) = 247.64`, raw `p = 2.24e-19`.
- The block Friedman test is significant at `43/46` bins; the exceptions are `215`, `225` and `250 nm`.
- The robust simultaneous band places `F2-F1 > 0` over `30-105 nm` plus the isolated `115 nm` bin, `F4-F2 < 0` over all `25-250 nm`, and `F4-F1 < 0` over `25-145`, `180-190` and `200-250 nm`.
- Each pair's 4x4-block exact functional max-T test has raw `p = 3.0518e-5` and Holm `p = 9.1553e-5`.
- Repeating the exact test with 8x8 blocks leaves only four blocks and gives `p = 0.125` for every pair. This is the minimum two-sided resolution of that four-block sign-flip test, not positive evidence that the distributions are equal.

Pixel-level t/Wilcoxon and marginal Welch/Mann-Whitney/KS tests are also saved, but their nominal p-values treat spatial pixels as independent and are therefore secondary. They confirm visible separation but should not replace the spatial-block result.

## Comparison with the no-slip hydrodynamic scale

For water at 25.6 °C (`eta = 0.880645 mPa s`) and sphere radius `R = 4.54685 um`, the no-slip lubrication estimate is

`F_hyd = 6 pi eta R^2 U / D`.

At 100 nm it is `3.4318 pN` per `1 um/s`; at 50 nm it is `6.8636 pN` per `1 um/s`. Consequently, no slip predicts:

| D (nm) | expected F2-F1 (pN) | observed paired median (pN) | expected F4-F2 (pN) | observed paired median (pN) |
|---:|---:|---:|---:|---:|
| 50 | +6.86 | +52.95 | +13.73 | -128.53 |
| 100 | +3.43 | +16.90 | +6.86 | -44.71 |

The `2-1` contrast has the expected sign but is roughly five to eight times the ideal no-slip increment; the `4-2` contrast reverses sign and is much larger in magnitude. A single positive hydrodynamic coefficient cannot produce both. This agrees with the earlier pairwise shape fits: all later-slower comparisons fitted positive `chi`, whereas every comparison ending at the chronologically last 4 um/s map fitted negative `chi`.

## Other measured changes that can shift F-D curves

In chronological `2 -> 1 -> 4 um/s` order at 0 wt%:

- median local contact InvOLS changes `66.5 -> 66.7 -> 67.1 nm/V`;
- terminal load remains close but changes `17.8 -> 17.8 -> 17.7 nN`;
- median far-field slopes are `-3.21 -> -4.02 -> -7.01 pN per 100 nm` before their per-curve subtraction;
- detected snap-in fractions are `5.1% -> 17.6% -> 7.8%`, with detected-event medians `10.4 -> 6.58 -> 5.40 nm`;
- median exact retract pull-off among uncensored/detected curves changes `-131 -> -105 -> -84.4 nN`.

These changes are not all mutually consistent with a single velocity mechanism. They demonstrate simultaneous evolution of contact/optical response, event detection and adhesion during the sequential maps. The far-field linear correction removes the fitted affine baseline, but it cannot remove nonlinear drift or a real change in surface/probe state.

## Claim boundary and experimental consequence

Supported claim: **the three measured map F-D distributions are conditionally separated within the mapped area, with a strong nonmonotonic ordering in the 50-100 nm force-bearing range.**

Unsupported claim: **the separation is a causal velocity dependence, a measured hydrodynamic coefficient, or a validated zero-speed equilibrium extrapolation.** There is one map per speed, so the experimental velocity sample size is `n = 1` per group regardless of the 255 paired pixels.

The appropriate follow-up is the documented balanced 8x8 palindromic protocol, in which whole maps rather than pixels are the velocity replicates and acquisition order is counterbalanced. See [`../FOLLOWUP_8X8_PALINDROME_FORCE_MAPPING_PROTOCOL.md`](../FOLLOWUP_8X8_PALINDROME_FORCE_MAPPING_PROTOCOL.md).

## Reproducible artifacts

- [`ZERO_WT_DISTRIBUTION_SEPARATION.md`](ZERO_WT_DISTRIBUTION_SEPARATION.md): robust effects, pointwise/simultaneous intervals and functional max-T tests.
- [`ZERO_WT_CLASSICAL_DISTRIBUTION_TESTS.md`](ZERO_WT_CLASSICAL_DISTRIBUTION_TESTS.md): classical pairwise and omnibus tests.
- [`zero_wt_distribution_separation_selected.csv`](zero_wt_distribution_separation_selected.csv) and [`zero_wt_classical_distribution_tests_selected.csv`](zero_wt_classical_distribution_tests_selected.csv): selected-distance numerical values quoted above.
- [`zero_wt_map_pixel_forces_25_250nm.csv`](zero_wt_map_pixel_forces_25_250nm.csv): all 255 matched pixels, three maps and 46 separations.
- [`figures/zero_wt_fd_highlight_1_um_s.png`](figures/zero_wt_fd_highlight_1_um_s.png), [`figures/zero_wt_fd_highlight_2_um_s.png`](figures/zero_wt_fd_highlight_2_um_s.png) and [`figures/zero_wt_fd_highlight_4_um_s.png`](figures/zero_wt_fd_highlight_4_um_s.png): all-curve views with one speed highlighted, other speeds gray, and pointwise median/interquartile range.
- [`figures/zero_wt_distribution_separation_confidence.png`](figures/zero_wt_distribution_separation_confidence.png) and [`figures/zero_wt_classical_distribution_test_pvalues.png`](figures/zero_wt_classical_distribution_test_pvalues.png): confidence and p-value summaries.
- `analysis/plot_zero_wt_all_fd_curves.py`, `analysis/analyze_zero_wt_distribution_separation.py` and `analysis/analyze_zero_wt_classical_distribution_tests.py`: generators, run in that order after the primary pipelines.
- [`zero_wt_distribution_analysis_manifest.sha256`](zero_wt_distribution_analysis_manifest.sha256): SHA-256 identities for all supplemental scripts, tables, reports, provenance records and figures listed here.

# Same-pixel velocity diagnostics and joint silica sphere-plane fit

## Result status

The raw-only pipeline decoded 43 JPK sources and retained 18 maps; 15 maps at 0, 20, 30 and 40 wt% form the primary analysis. Ten wt% is retained only as a measurement-regime QC control. Every force value uses a concentration-specific consensus calculated from all accepted raw hard-contact anchors and the calibrated cantilever-1 spring constant `k = 0.2969899087 N/m`; embedded JPK sensitivity and force are unused.

The main conclusion is evidence-bounded: separation-resolved velocity differences are real, but the current sequential-map design does not identify a unique hydrodynamic coefficient independently of map history. Same-speed 30 wt% controls retain a median linear-nuisance residual scale of `27.2 pN`, while only `5/14` nonzero-speed pair fits have a positive fitted no-slip amplitude ratio. The joint PB parameters below are therefore model-conditioned exploratory results and systematic ranges, not a validated zero-velocity extrapolation.

## Map inventory and load strata

| wt% | time | nominal speed (um/s) | retained pixels | load regime | terminal load (nN) | contact response |
|---:|:---|---:|---:|:---|---:|---:|
| 0 | 15:01:47 | 2 | 256 | main | 17.8 | 1.0018 |
| 0 | 15:06:26 | 1 | 256 | main | 17.8 | 0.9994 |
| 0 | 15:15:57 | 4 | 255 | main | 17.7 | 0.9924 |
| 10 | 15:29:51 | 2 | 63 | main | 0.875 | 0.9580 |
| 10 | 15:35:17 | 1 | 102 | main | 1.02 | 0.9838 |
| 10 | 15:44:18 | 4 | 112 | main | 0.878 | 1.0101 |
| 20 | 15:54:06 | 2 | 242 | main | 7.04 | 0.9672 |
| 20 | 15:59:55 | 1 | 251 | main | 7.2 | 1.0126 |
| 20 | 16:09:19 | 4 | 252 | main | 6.87 | 1.0269 |
| 30 | 16:20:49 | 2 | 149 | low_load | 6.67 | 0.9623 |
| 30 | 16:24:44 | 2 | 111 | low_load | 6.7 | 0.9650 |
| 30 | 16:27:34 | 2 | 206 | high_load | 18.4 | 0.9942 |
| 30 | 16:37:53 | 1 | 101 | high_load | 18.6 | 1.0133 |
| 30 | 16:42:02 | 1 | 179 | high_load | 18.6 | 1.0173 |
| 30 | 16:53:04 | 4 | 256 | high_load | 18.4 | 1.0205 |
| 40 | 17:03:03 | 2 | 255 | main | 18.1 | 1.0078 |
| 40 | 17:07:48 | 1 | 255 | main | 18.3 | 1.0040 |
| 40 | 17:18:39 | 4 | 256 | main | 18.1 | 0.9644 |

At 30 wt%, the first two 2 um/s maps are a roughly 6.7 nN low-load stratum. The later 2, 1, 1 and 4 um/s maps are a roughly 18.5 nN high-load stratum. Cross-stratum differences are not used as velocity contrasts.

## Same-pixel event and contact comparisons

All differences are `second map - first map` at the same physical raster coordinate after undoing serpentine row reversal. Spatial-block confidence intervals use 4x4 pixel blocks and quantify only within-area heterogeneity. Snap differences require a detected event in both maps, so their paired `n` and interval are detection-conditioned; missing/no-snap pixels are not silently assigned zero distance.

| wt% | load | speed pair (um/s) | metric | paired n | median difference | block CI95 | acquisition-order r | plane R2 |
|---:|:---|:---|:---|---:|---:|:---|---:|---:|
| 0 | main | 2->1 | contact_invOLS_nm_per_V | 256 | 0.196 nm/V | [0.0459, 0.26] | -0.0538 | 0.00827 |
| 0 | main | 2->1 | terminal_load_nN | 256 | 0.00627 nN | [-0.0261, 0.027] | -0.0742 | 0.0087 |
| 0 | main | 2->1 | snap_distance_nm | 1 | -2.09 nm | [n/a, n/a] | n/a | n/a |
| 0 | main | 2->4 | contact_invOLS_nm_per_V | 255 | 0.749 nm/V | [0.588, 0.896] | 0.061 | 0.00757 |
| 0 | main | 2->4 | terminal_load_nN | 255 | -0.0364 nN | [-0.0542, -0.019] | -0.0266 | 0.000902 |
| 0 | main | 2->4 | snap_distance_nm | 2 | -6.06 nm | [n/a, n/a] | n/a | n/a |
| 0 | main | 1->4 | contact_invOLS_nm_per_V | 255 | 0.578 nm/V | [0.413, 0.72] | 0.101 | 0.00458 |
| 0 | main | 1->4 | terminal_load_nN | 255 | -0.0353 nN | [-0.0708, -0.0174] | 0.0582 | 0.00649 |
| 0 | main | 1->4 | snap_distance_nm | 6 | -1.83 nm | [n/a, n/a] | 0.465 | n/a |
| 20 | main | 2->1 | contact_invOLS_nm_per_V | 237 | -2.6 nm/V | [-3, -2.15] | 0.038 | 0.064 |
| 20 | main | 2->1 | terminal_load_nN | 237 | 0.175 nN | [0.116, 0.25] | -0.321 | 0.104 |
| 20 | main | 2->1 | snap_distance_nm | 93 | 2.33 nm | [1.11, 3.05] | -0.196 | 0.0432 |
| 20 | main | 2->4 | contact_invOLS_nm_per_V | 238 | -3.44 nm/V | [-4.07, -2.82] | 0.0617 | 0.00482 |
| 20 | main | 2->4 | terminal_load_nN | 238 | -0.157 nN | [-0.227, -0.12] | -0.194 | 0.0347 |
| 20 | main | 2->4 | snap_distance_nm | 112 | 0.265 nm | [-0.722, 1.32] | -0.155 | 0.0703 |
| 20 | main | 1->4 | contact_invOLS_nm_per_V | 249 | -0.758 nm/V | [-1.26, -0.493] | 0.0319 | 0.0119 |
| 20 | main | 1->4 | terminal_load_nN | 249 | -0.336 nN | [-0.394, -0.301] | 0.221 | 0.0431 |
| 20 | main | 1->4 | snap_distance_nm | 158 | -1.69 nm | [-2.52, -0.989] | 0.0937 | 0.0369 |
| 30 | low_load | 2->2 | contact_invOLS_nm_per_V | 110 | -0.0911 nm/V | [-0.588, 0.362] | 0.304 | 0.0415 |
| 30 | low_load | 2->2 | terminal_load_nN | 110 | 0.0409 nN | [0.00295, 0.0681] | 0.325 | 0.0914 |
| 30 | low_load | 2->2 | snap_distance_nm | 100 | 1.5 nm | [-2, 3.85] | -0.181 | 0.0476 |
| 30 | high_load | 2->1 | contact_invOLS_nm_per_V | 101 | -1.65 nm/V | [-1.85, -1.45] | -0.245 | 0.0742 |
| 30 | high_load | 2->1 | terminal_load_nN | 101 | 0.213 nN | [0.159, 0.246] | 0.0512 | 0.000273 |
| 30 | high_load | 2->1 | snap_distance_nm | 73 | -1.88 nm | [-3.8, -0.567] | 0.0223 | 0.00111 |
| 30 | high_load | 2->1 | contact_invOLS_nm_per_V | 177 | -1.69 nm/V | [-1.93, -0.875] | 0.324 | 0.291 |
| 30 | high_load | 2->1 | terminal_load_nN | 177 | 0.201 nN | [0.131, 0.224] | -0.125 | 0.0223 |
| 30 | high_load | 2->1 | snap_distance_nm | 128 | -2.8 nm | [-3.76, -1.56] | -0.0896 | 0.0213 |
| 30 | high_load | 2->4 | contact_invOLS_nm_per_V | 206 | -1.34 nm/V | [-1.76, -0.709] | 0.643 | 0.673 |
| 30 | high_load | 2->4 | terminal_load_nN | 206 | 0.00726 nN | [-0.0251, 0.0549] | -0.156 | 0.0142 |
| 30 | high_load | 2->4 | snap_distance_nm | 154 | -5.79 nm | [-7.35, -4.53] | -0.0382 | 0.0163 |
| 30 | high_load | 1->1 | contact_invOLS_nm_per_V | 100 | -0.287 nm/V | [-0.362, -0.149] | -0.17 | 0.0867 |
| 30 | high_load | 1->1 | terminal_load_nN | 100 | 0.0565 nN | [-0.0389, 0.127] | 0.119 | 0.114 |
| 30 | high_load | 1->1 | snap_distance_nm | 65 | 0.709 nm | [-0.648, 1.91] | 0.152 | 0.0763 |
| 30 | high_load | 1->4 | contact_invOLS_nm_per_V | 101 | -0.206 nm/V | [-0.634, 0.0426] | 0.4 | 0.245 |
| 30 | high_load | 1->4 | terminal_load_nN | 101 | -0.127 nN | [-0.174, -0.0936] | 0.0101 | 0.0473 |
| 30 | high_load | 1->4 | snap_distance_nm | 63 | -3.34 nm | [-3.62, -2.33] | 0.089 | 0.0313 |
| 30 | high_load | 1->4 | contact_invOLS_nm_per_V | 179 | -0.0413 nm/V | [-0.215, 0.266] | 0.0992 | 0.0172 |
| 30 | high_load | 1->4 | terminal_load_nN | 179 | -0.188 nN | [-0.226, -0.143] | 0.0341 | 0.0139 |
| 30 | high_load | 1->4 | snap_distance_nm | 113 | -3.57 nm | [-4.49, -2.91] | 0.0813 | 0.0437 |
| 40 | main | 2->1 | contact_invOLS_nm_per_V | 255 | 0.169 nm/V | [0.0356, 0.315] | -0.0477 | 0.015 |
| 40 | main | 2->1 | terminal_load_nN | 255 | 0.207 nN | [0.164, 0.24] | 0.159 | 0.0641 |
| 40 | main | 2->1 | snap_distance_nm | 124 | 0.298 nm | [-0.239, 0.943] | 0.126 | 0.064 |
| 40 | main | 2->4 | contact_invOLS_nm_per_V | 255 | 2.57 nm/V | [2.43, 2.75] | 0.0349 | 0.124 |
| 40 | main | 2->4 | terminal_load_nN | 255 | 0.063 nN | [0.0305, 0.108] | -0.0905 | 0.0391 |
| 40 | main | 2->4 | snap_distance_nm | 94 | -3.23 nm | [-3.98, -2.54] | -0.143 | 0.0722 |
| 40 | main | 1->4 | contact_invOLS_nm_per_V | 255 | 2.33 nm/V | [2.22, 2.55] | 0.085 | 0.14 |
| 40 | main | 1->4 | terminal_load_nN | 255 | -0.151 nN | [-0.211, -0.0948] | -0.237 | 0.0666 |
| 40 | main | 1->4 | snap_distance_nm | 99 | -4.19 nm | [-4.6, -3.32] | -0.0996 | 0.0119 |

The strongest per-pair correlations between a snap-position change and a nuisance-metric change are diagnostic, not causal:

| wt% | load | speed pair | predictor | n | Spearman rho |
|---:|:---|:---|:---|---:|---:|
| 30 | high_load | 1->4 | delta_contact_height_nm | 113 | -0.201 |
| 30 | high_load | 2->1 | delta_terminal_load_nN | 128 | 0.199 |
| 40 | main | 2->4 | delta_contact_invOLS_nm_per_V | 94 | -0.191 |
| 30 | low_load | 2->2 | delta_contact_height_nm | 100 | 0.176 |
| 40 | main | 2->1 | delta_contact_invOLS_nm_per_V | 124 | -0.173 |
| 20 | main | 1->4 | delta_contact_height_nm | 158 | -0.151 |
| 30 | high_load | 1->1 | delta_contact_height_nm | 65 | 0.149 |
| 20 | main | 2->1 | delta_terminal_load_nN | 93 | 0.146 |

## Separation-resolved hydrodynamic shape test

For this test, each trace is referenced by only a far-field voltage constant; its far-field slope is retained. Each paired force difference is fitted as `Delta F = a + b(D-D0) + chi Delta[6 pi eta R^2 U/D]`. Thus the linear instrumental drift and the physical `1/D` shape compete in the same fit instead of the latter being removed beforehand.

| wt% | load | speed pair | chi | local CI95 | weighted RSS improvement | detrended sign agreement | full RMSE (pN) | status |
|---:|:---|:---|---:|:---|---:|---:|---:|:---|
| 0 | main | 2->1 | 6.45 | [3.93, 8.97] | 0.368 | 0.667 | 9.71 | descriptive_shape_test_speed_confounded_with_map_history |
| 0 | main | 2->4 | -14.6 | [-16.3, -12.9] | 0.868 | 0 | 19.8 | descriptive_shape_test_speed_confounded_with_map_history |
| 0 | main | 1->4 | -5.88 | [-6.3, -5.47] | 0.947 | 0 | 6.15 | descriptive_shape_test_speed_confounded_with_map_history |
| 20 | main | 2->1 | 6.85 | [5.08, 8.62] | 0.572 | 0.933 | 12.9 | descriptive_shape_test_speed_confounded_with_map_history |
| 20 | main | 2->4 | -4.52 | [-5.97, -3.07] | 0.463 | 0.167 | 22.6 | descriptive_shape_test_speed_confounded_with_map_history |
| 20 | main | 1->4 | -0.126 | [-0.804, 0.552] | 0.00306 | 0.5 | 12 | descriptive_shape_test_speed_confounded_with_map_history |
| 30 | low_load | 2->2 | n/a | [n/a, n/a] | 0 | n/a | 36.1 | same_speed_history_control |
| 30 | high_load | 2->1 | 14.8 | [12.6, 16.9] | 0.804 | 1 | 33.7 | descriptive_shape_test_speed_confounded_with_map_history |
| 30 | high_load | 2->1 | 25.7 | [23.4, 28] | 0.916 | 1 | 64.7 | descriptive_shape_test_speed_confounded_with_map_history |
| 30 | high_load | 2->4 | -16.8 | [-17.5, -16.1] | 0.982 | 0 | 20.2 | descriptive_shape_test_speed_confounded_with_map_history |
| 30 | high_load | 1->1 | n/a | [n/a, n/a] | 0 | n/a | 18.2 | same_speed_history_control |
| 30 | high_load | 1->4 | -4.41 | [-4.94, -3.89] | 0.863 | 0 | 9.94 | descriptive_shape_test_speed_confounded_with_map_history |
| 30 | high_load | 1->4 | -1.31 | [-1.53, -1.09] | 0.756 | 0 | 4.96 | descriptive_shape_test_speed_confounded_with_map_history |
| 40 | main | 2->1 | 3.14 | [2.31, 3.96] | 0.565 | 1 | 10.7 | descriptive_shape_test_speed_confounded_with_map_history |
| 40 | main | 2->4 | -5.21 | [-5.85, -4.57] | 0.856 | 0 | 21.1 | descriptive_shape_test_speed_confounded_with_map_history |
| 40 | main | 1->4 | -1.68 | [-1.9, -1.46] | 0.841 | 0 | 10.2 | descriptive_shape_test_speed_confounded_with_map_history |

The sign pattern is locked to acquisition order rather than a common velocity law: `5/5` comparisons in which the later map is slower have positive `chi`, whereas `9/9` comparisons ending at the chronologically later 4 um/s map have negative `chi`. A no-slip force cannot reverse its fitted sign merely because the speed order is changed; this is direct evidence for a superposed time/surface-history component.

A physical no-slip result would require positive `chi` near one, a positive nested weighted-RSS improvement beyond the linear drift alone, viscosity/velocity scaling, and an approach/retract sign reversal. Sign agreement is evaluated only after projecting both the measured difference and hydrodynamic basis away from the same intercept/linear-distance nuisance. These pairwise local CIs treat separation bins as independent and are deliberately not promoted to experimental confidence intervals.

## Approach/retract odd-even feasibility

Only `0` primary maps retain at least 12 bins with 50% verified free-retract support inside 25-250 nm. The adhered retract branch is never substituted for a free non-contact branch; encoder-floor values and right-censored detachments remain excluded.

| wt% | time | speed | minimum D with 50% free support (nm) | qualified bins 25-250 | usable for primary window |
|---:|:---|---:|---:|---:|:---|
| 0 | 15:01:47 | 2 | 440 | 0 | False |
| 0 | 15:06:26 | 1 | 345 | 0 | False |
| 0 | 15:15:57 | 4 | 285 | 0 | False |
| 20 | 15:54:06 | 2 | 620 | 0 | False |
| 20 | 15:59:55 | 1 | 580 | 0 | False |
| 20 | 16:09:19 | 4 | 580 | 0 | False |
| 30 | 16:20:49 | 2 | n/a | 0 | False |
| 30 | 16:24:44 | 2 | n/a | 0 | False |
| 30 | 16:27:34 | 2 | n/a | 0 | False |
| 30 | 16:37:53 | 1 | n/a | 0 | False |
| 30 | 16:42:02 | 1 | n/a | 0 | False |
| 30 | 16:53:04 | 4 | n/a | 0 | False |
| 40 | 17:03:03 | 2 | 605 | 0 | False |
| 40 | 17:07:48 | 1 | 645 | 0 | False |
| 40 | 17:18:39 | 4 | 640 | 0 | False |

## Joint same-surface PB fits

The fitted observation model is `Fobs = gm [Feq(D; |zeta|, lambda_D) + chi Fhyd(D,U)] + am + bm(D-D0)`. `Feq` is nonlinear equal-constant-potential 1:1 PB plus sphere-plane van der Waals force, converted with Derjaguin. Each map's separate raw hard-contact response enters the joint likelihood as an observation of `gm` with a minimum 2% uncertainty; it is not a per-curve force recalibration. Zeta and Debye length are shared across speeds within a concentration/load stratum.

| wt% | model | lambda_D (nm) | |zeta| (mV) | chi | R2 | RMSE (pN) | Delta AICc | numerical status |
|---:|:---|---:|---:|---:|---:|---:|---:|:---|
| 0 | M0_no_hydrodynamic | 18.6 | 74 | 0 | 0.9994 | 10.1 | 16.6 | valid |
| 0 | M1_no_slip_chi1 | 18.5 | 73.7 | 1 | 0.9994 | 10.1 | 24.3 | valid |
| 0 | M2_fitted_chi | 18.9 | 75.2 | -4.05 | 0.9994 | 9.74 | 0 | valid |
| 20 | M0_no_hydrodynamic | 18.1 | 78.1 | 0 | 0.9988 | 14.2 | 0 | valid |
| 20 | M1_no_slip_chi1 | 17.9 | 77.5 | 1 | 0.9988 | 14.3 | 4.09 | valid |
| 20 | M2_fitted_chi | 18.1 | 78.2 | -0.15 | 0.9988 | 14.2 | 2.33 | valid |
| 30 | M0_no_hydrodynamic | 15.9 | 78.6 | 0 | 0.9961 | 22.4 | 11.1 | valid |
| 30 | M1_no_slip_chi1 | 15.7 | 78.1 | 1 | 0.9957 | 23.5 | 35.5 | valid |
| 30 | M2_fitted_chi | 16.3 | 79.3 | -1.55 | 0.9967 | 20.6 | 0 | valid |
| 40 | M0_no_hydrodynamic | 20.3 | 80.9 | 0 | 0.9992 | 11.9 | 47.4 | valid |
| 40 | M1_no_slip_chi1 | 20 | 79.7 | 1 | 0.9992 | 12.3 | 82.5 | valid |
| 40 | M2_fitted_chi | 21.5 | 85.8 | -3.42 | 0.9993 | 11.5 | 0 | valid |

The negative silica sign is assigned chemically; identical-surface normal-force data determine only `|zeta|`. These parameters are boundary-potential fit parameters. Equating them to electrokinetic zeta additionally assumes the boundary potential represents the slipping-plane potential.

## Model and perturbation systematic ranges

| wt% | preferred tested primary model | lambda all valid tests (nm) | |zeta| all valid tests (mV) | M2 chi range | claim boundary |
|---:|:---|:---|:---|:---|:---|
| 0 | M2_fitted_chi | 18.1-20.1 | 68.3-82.3 | -12.6--0.0626 | exploratory_systematic_range_not_velocity_identified_parameter |
| 20 | M0_no_hydrodynamic | 17.3-18.7 | 70.9-86.3 | -1.88-18.8 | exploratory_systematic_range_not_velocity_identified_parameter |
| 30 | M2_fitted_chi | 15.5-16.8 | 71.4-88.3 | -3.54-9.04 | exploratory_systematic_range_not_velocity_identified_parameter |
| 40 | M2_fitted_chi | 18.9-24.1 | 70.8-99.3 | -9.03-11 | exploratory_systematic_range_not_velocity_identified_parameter |

The ranges above include M0/M1/M2, fit-window changes, fixed-unity versus contact-constrained force scale, a ±2.5 nm contact-zero shift, and leave-one-nominal-speed-out M2 fits, retaining only numerically valid results. The contact shift is half one 5 nm analysis bin and is a resolution-scale sensitivity test, not a calibrated contact-position confidence interval.

### Contact-zero and leave-one-speed-out detail

| wt% | perturbation | excluded speed (um/s) | contact shift (nm) | lambda_D (nm) | |zeta| (mV) | chi | status |
|---:|:---|---:|---:|---:|---:|---:|:---|
| 0 | M2_contact_minus2p5nm | n/a | -2.5 | 19.4 | 68.3 | -4.71 | valid |
| 0 | M2_contact_plus2p5nm | n/a | 2.5 | 19.2 | 82.3 | -6.38 | valid |
| 0 | M2_leave_out_1um_s | 1 | 0 | 20.1 | 79.1 | -12.6 | valid |
| 0 | M2_leave_out_2um_s | 2 | 0 | 18.1 | 75.1 | -2.88 | valid |
| 0 | M2_leave_out_4um_s | 4 | 0 | 17.9 | 69.8 | 20 | invalid: parameter_boundary |
| 20 | M2_contact_minus2p5nm | n/a | -2.5 | 18.5 | 71.8 | -1.31 | valid |
| 20 | M2_contact_plus2p5nm | n/a | 2.5 | 18.4 | 86.3 | -1.13 | valid |
| 20 | M2_leave_out_1um_s | 1 | 0 | 17.4 | 77.2 | 1.88 | valid |
| 20 | M2_leave_out_2um_s | 2 | 0 | 17.9 | 77.8 | -0.148 | valid |
| 20 | M2_leave_out_4um_s | 4 | 0 | 17.9 | 70.9 | 18.8 | valid |
| 30 | M2_contact_minus2p5nm | n/a | -2.5 | 16.3 | 71.4 | -1.39 | valid |
| 30 | M2_contact_plus2p5nm | n/a | 2.5 | 16.3 | 88.3 | -1.63 | valid |
| 30 | M2_leave_out_1um_s | 1 | 0 | 16.8 | 84.4 | -3.44 | valid |
| 30 | M2_leave_out_2um_s | 2 | 0 | 15.8 | 75.5 | -1.12 | valid |
| 30 | M2_leave_out_4um_s | 4 | 0 | 15.5 | 73.9 | 9.04 | valid |
| 40 | M2_contact_minus2p5nm | n/a | -2.5 | 21.6 | 77.6 | -2.78 | valid |
| 40 | M2_contact_plus2p5nm | n/a | 2.5 | 22 | 93.8 | -4.38 | valid |
| 40 | M2_leave_out_1um_s | 1 | 0 | 24.1 | 99.3 | -9.03 | valid |
| 40 | M2_leave_out_2um_s | 2 | 0 | 19.2 | 82.6 | -0.769 | valid |
| 40 | M2_leave_out_4um_s | 4 | 0 | 18.9 | 70.8 | 11 | valid |

AICc uses the Gaussian joint weighted likelihood of the force bins and separate map contact-scale observations: because all residuals are normalized by held-fixed uncertainty estimates, the model-dependent `-2 log L` term is the joint weighted chi-square. M2 is genuinely nested with M0 at `chi=0`. The comparison remains conditional on correlated map-median bins and does not remove the experimental speed/time confounding. A preferred model is therefore the best descriptor among M0/M1/M2, not proof that its velocity term is causal.

## Numerical semantics and checks

- Sphere radius: `4.54685 um`; primary maximum `D/R = 0.0550`, supporting the leading lubrication/Derjaguin small-gap approximation as a controlled diagnostic rather than an exact finite-separation law.
- Water no-slip force at 1 um/s and 100 nm: `3.43181 pN`; inverse-distance and linear-speed self-test errors are `0.000e+00` and `0.000e+00`.
- Synthetic gap-speed maximum error: `8.646e-12 um/s`; synthetic pairwise chi error: `4.441e-16`.
- Nonlinear PB linear-limit and far-asymptote maximum relative errors: `1.652e-03` and `1.054e-03`.
- Optimizer success is supplemented by Jacobian rank/condition, parameter-boundary and R2 checks. Local Jacobian intervals remain conditional because distance bins from one map are correlated.
- All random resampling uses the recorded local seed; every raw ZIP passes CRC and is SHA-256 listed.

## Files

- `paired_metric_differences.csv`, `paired_metric_summary.csv`: same-pixel event/contact comparisons and block intervals.
- `snap_nuisance_associations.csv`: paired snap changes versus contact, height, load and far-drift changes.
- `map_force_curves.csv`, `paired_force_differences.csv`: map medians and same-pixel full-curve contrasts.
- `pairwise_hydrodynamic_fits.csv`: `1/D` hydrodynamic shape test with simultaneous linear drift.
- `branch_odd_even.csv`, `branch_odd_even_summary.csv`: verified free-retract support and odd/even candidates.
- `joint_fit_results.csv`, `joint_fit_predictions.csv`, `joint_fit_map_nuisance.csv`: M0/M1/M2 results and components.
- `parameter_systematic_summary.csv`, `map_inventory.csv`, `input_manifest.csv`, figures, provenance and SHA-256 manifest.

# Velocity-correlated AFM systematics

## Scope and physical definitions

This analysis re-decodes raw `vDeflection` and `measuredHeight` for paired approach/retract branches. Embedded JPK sensitivity and force are not used. One raw hard-contact sensitivity consensus is recomputed per glycerol concentration and force is rebuilt with the calibrated cantilever-1 spring constant.

- **Approach snap-in distance**: contact-aligned `h + delta - h_contact` immediately before a statistically resolved downward force step. This is a tip-surface distance estimate.
- **Retract detachment position**: scanner/piezo travel `h_detach - h_contact` immediately before snap-off. It is not called a true gap because the probe remains adhered before release.
- **Pull-off force**: minimum reconstructed force on the adhered retract branch between snap-off and hard contact.
- **Contact response ratio**: common concentration sensitivity divided by the local contact InvOLS. A rigid, unchanged optical/contact response is near one; systematic changes flag calibration/contact mechanics bias.
- Velocity interpretation uses whole maps as experimental units. Pixel distributions characterize within-map heterogeneity only. The 10 wt% group is retained as anomaly QC and excluded from primary interpretation.

## Main findings

1. **Snap-in changes, but there is no concentration-independent monotonic velocity law.** The 4 um/s map has the smallest median detected snap-in distance in 0, 30 and 40 wt%, which is qualitatively compatible with hydrodynamic squeeze pressure delaying an instability. However, 20 wt% has nearly the same 2 and 4 um/s value and its 1 um/s value is largest; detection fractions also change strongly between maps. The result is therefore a possible dynamic contribution superposed on larger history/calibration effects, not an independently fitted velocity effect.

2. **A fixed-speed control directly shows time/history drift.** In the three consecutive 30 wt%, 2 um/s maps, median snap-in distance changes from 17.5 to 13.8 nm, while local contact InvOLS changes from 61 to 59 nm/V. The terminal load simultaneously changes from 6.67 to 18.4 nN. These changes cannot be caused by velocity because velocity is unchanged.

3. **The drift also occurs inside a single map.** The largest robust contact-InvOLS change across one primary map is -2.92 nm/V (30 wt%, 16:27:34), with spatial/acquisition plane R2 0.820. A fixed-speed map therefore contains a systematic optical/contact-response gradient.

4. **Retract adhesion is heavily censored by the raw encoder floor.** This is not a subtle statistical effect: many raw signed-int32 traces equal the lower encoding bound. Exact pull-off values are therefore available only for unclipped adhered branches. In particular, most 30 wt% traces do not contain a verified unsaturated free-cantilever prefix, so their detachment is not observed within the usable retract record; reporting near-zero pull-off for them would be wrong.

5. **Independent force curves also expose history dependence.** At 0 wt% and the same approximately 2 um/s speed, exact pull-off changes from -132 nN at 14:59:33 to -76.5 nN at 15:18:48; detachment travel changes from 431 to 248 nm. This same-speed change is point heterogeneity, surface/probe conditioning, or time drift—not hydrodynamic velocity dependence.

6. **The 10 wt% block is a different measurement regime.** Its approach hard-contact validity is low and terminal loads are about 1 nN rather than the roughly 7 or 18 nN regimes in the primary data. Retaining it as QC supports its exclusion from the primary surface-force fit.

7. **Separation zero has a large map plane.** The largest robust hard-contact-height change across a primary map is 84.4 nm (30 wt%, 16:20:49), with plane R2 0.982. The repeatable column component is compatible with sample tilt/topography, while changing row components can also contain piezo creep or time drift. Per-pixel contact alignment removes the first-order offset; a single common contact zero would not.

8. **Far-field photodiode drift is also map dependent and nonmonotonic.** At 20 wt% its chronological map medians are 12.1, 33.6 and 4.5 pN per 100 nm for the 2, 1 and 4 um/s maps. The robust line is subtracted independently from every approach trace, so this is retained as instrumental QC rather than mistaken for a surface force.

## Map-level results

| wt% | time | speed | contact valid | local InvOLS (nm/V) | snap detected | snap D (nm) | terminal load (nN) | retract clipped | free baseline | snap-off detected | exact pull-off (nN) | exact n / total | detach piezo (nm) |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 15:01:47 | 2 | 100.0% | 66.5 | 5.1% | 10.4 | 17.8 | 0.0% | 100.0% | 96.9% | -131 | 248/256 | 430 |
| 0 | 15:06:26 | 1 | 100.0% | 66.7 | 17.6% | 6.58 | 17.8 | 0.0% | 100.0% | 98.4% | -105 | 252/256 | 342 |
| 0 | 15:15:57 | 4 | 99.6% | 67.1 | 7.8% | 5.4 | 17.7 | 0.0% | 100.0% | 98.4% | -84.4 | 252/256 | 274 |
| 10 | 15:29:51 | 2 | 24.6% | 67 | 10.2% | 7.56 | 0.856 | 0.0% | 100.0% | 61.3% | -0.265 | 157/256 | 26.4 |
| 10 | 15:35:17 | 1 | 39.8% | 65.2 | 11.7% | 7.87 | 1 | 0.0% | 100.0% | 64.1% | -2.08 | 164/256 | 30.6 |
| 10 | 15:44:18 | 4 | 43.8% | 63.5 | 14.5% | 7.96 | 0.88 | 0.0% | 100.0% | 54.7% | -3.17 | 140/256 | 38.7 |
| 20 | 15:54:06 | 2 | 94.5% | 62.8 | 47.7% | 8.1 | 7.04 | 37.5% | 95.7% | 94.1% | -155 | 157/256 | 597 |
| 20 | 15:59:55 | 1 | 98.0% | 60 | 68.0% | 10.2 | 7.2 | 58.2% | 97.7% | 97.3% | -137 | 107/256 | 568 |
| 20 | 16:09:19 | 4 | 98.4% | 59.1 | 84.4% | 7.92 | 6.87 | 75.0% | 100.0% | 100.0% | -106 | 64/256 | 560 |
| 30 | 16:20:49 | 2 | 100.0% | 61 | 94.6% | 17.5 | 6.67 | 99.3% | 12.1% | 12.1% | -79.8 | 1/149 | 865 |
| 30 | 16:24:44 | 2 | 100.0% | 60.8 | 97.3% | 16.3 | 6.7 | 99.1% | 7.2% | 7.2% | -7.48 | 1/111 | 821 |
| 30 | 16:27:34 | 2 | 100.0% | 59 | 88.5% | 13.8 | 18.4 | 99.5% | 5.8% | 5.8% | -67.8 | 1/208 | 785 |
| 30 | 16:37:53 | 1 | 100.0% | 57.9 | 80.2% | 10 | 18.6 | 99.0% | 24.8% | 24.8% | -1.57 | 1/101 | 821 |
| 30 | 16:42:02 | 1 | 99.4% | 57.7 | 77.3% | 10.9 | 18.6 | 99.4% | 5.5% | 5.5% | -3.77 | 1/181 | 824 |
| 30 | 16:53:04 | 4 | 100.0% | 57.5 | 83.2% | 7.87 | 18.4 | 99.2% | 9.8% | 9.8% | -137 | 2/256 | 799 |
| 40 | 17:03:03 | 2 | 99.6% | 57.5 | 64.8% | 9.15 | 18.1 | 11.3% | 98.8% | 98.4% | -179 | 226/256 | 597 |
| 40 | 17:07:48 | 1 | 99.6% | 57.7 | 71.1% | 9.64 | 18.3 | 53.5% | 94.5% | 94.5% | -171 | 119/256 | 637 |
| 40 | 17:18:39 | 4 | 100.0% | 60 | 55.5% | 5.83 | 18.1 | 43.0% | 96.5% | 96.5% | -168 | 146/256 | 631 |

The table is descriptive. In conditions with one map at each speed, speed is confounded with acquisition time and surface history. The six 30 wt% maps are unbalanced (three at 2 um/s, two at 1 um/s and one at 4 um/s), so they provide fixed-speed drift controls but not balanced velocity replication.

`Retract clipped` means that raw signed-int32 `vDeflection` reached the JPK encoder floor. Such pull-off amplitudes are censored and are not included as exact force values. `Free baseline` requires at least 80 leading unsaturated samples after normalizing retract to far-to-contact orientation. If the physical retract ends at the encoder floor, detachment was not verifiably observed within the available travel; its position and pull-off force are right/left censored rather than zero.
For velocity summaries/regressions, retract event metrics are used only when at least 50% of a map has an uncensored/detected value. The per-source table retains smaller subsets visibly, with their exact numerator and total, for QC only.

## Physical interpretation

Quasi-statically, snap-in is a cantilever instability when the attractive interaction-force gradient reaches the restoring stiffness. A finite approach rate can shift that event because hydrodynamic pressure and cantilever dynamics alter the instantaneous force balance; dynamic drive-velocity effects are treated explicitly by [Bowen and Cheneler, Langmuir 2012](https://doi.org/10.1021/la304009c), and hydrodynamic pressure delaying jump-to-contact in liquid is analysed by [Dai, Journal of Fluid Mechanics 2025](https://doi.org/10.1017/jfm.2025.10290). Thus a velocity-dependent snap position is physically possible, but only after contact-zero, optical response, load and history are controlled.

Approach/retract comparison can isolate a hydrodynamic component when the non-hydrodynamic surface force is reproducible between branches and speeds; an example with silica sphere/plate is [McNamee et al., Colloids and Surfaces A 2007](https://doi.org/10.1016/j.colsurfa.2007.01.047). The present data violate that prerequisite through fixed-speed drift, load-regime changes and retract clipping. Large or gradual silica-silica adhesion in water has also been reported and attributed provisionally to bubbles/cavities by [Troncoso et al., Journal of Colloid and Interface Science 2014](https://doi.org/10.1016/j.jcis.2014.03.020); the current curves alone cannot identify that mechanism.

## Event-detection and numerical semantics

Approach events require a two-sample median-filtered drop larger than max(8 x far-field step MAD, 5 x far-field force residual MAD, 50 pN). Retract snap-off uses max(10 x step MAD, 5 x residual MAD, 25 pN). Candidate values are retained per curve, but detected-event summaries include only threshold-passing events.

Decoded 43 raw sources and 4103 approach curves. Synthetic step recovery: `True`; smooth no-event rejection: `True`; serpentine map-coordinate check: `True`.

## Files

- `curve_event_metrics.csv`: every paired branch, event candidates/detections, contact response, loads, drift/noise and explicit coordinates.
- `source_event_summary.csv`: robust per-source medians/MADs; map is the velocity-comparison unit.
- `within_map_metric_trends.csv`: acquisition-order correlations and spatial-plane gradients.
- `concentration_speed_summary.csv`: speed summaries of map medians, not pooled pixels.
- `velocity_time_models.csv`: descriptive map-level speed/time correlations and, where possible, a two-predictor model.
- `input_manifest.csv`, `provenance.json`, figures and SHA-256 artifact manifest.

# 31-08-26 五张 force map：map2–5 同地点四时点比较

## 直接结论

全部文件按timestamp依次定义为map1–map5。map2–map5的64个XY坐标逐点完全相同（最大差 `0 nm`），相对开始时间为 **0.00、27.16、53.33、78.66 min**。四图field、scan order和nominal speed相同，因此构成same-location四时点时间序列。

Force scale暂按metadata fingerprint最近的 **D5**、25.6 °C空气独立标定 `k=0.311340173 N/m`（repeatability SD `0.021740 N/m`）。五张map全部hard-contact共同得到global liquid InvOLS **68.290 nm/V**；文件写入的InvOLS与force conversion未用于最终force。

## Map inventory

| map | start time | relative time (min) | field | XY center (µm) | scanner/gap speed (µm/s) | local/global InvOLS (nm/V) | terminal load (nN) | far slope (pN/100 nm) |
|---:|:---|---:|:---|:---|:---|:---|---:|---:|
| 1 | 14:40:37 | — | 3.5×3.5 µm | (10.02, -18.94) | 1.000 / 0.999 | 63.137 / 68.290 | 14.130 | -0.64 |
| 2 | 15:01:59 | 0.00 | 2.5×2.5 µm | (11.89, -26.59) | 2.001 / 1.992 | 71.575 / 68.290 | 14.129 | -2.11 |
| 3 | 15:29:09 | 27.16 | 2.5×2.5 µm | (11.89, -26.59) | 2.001 / 1.995 | 70.877 / 68.290 | 14.121 | -2.01 |
| 4 | 15:55:19 | 53.33 | 2.5×2.5 µm | (11.89, -26.59) | 2.001 / 1.998 | 70.637 / 68.290 | 14.136 | -2.02 |
| 5 | 16:20:38 | 78.66 | 2.5×2.5 µm | (11.89, -26.59) | 2.001 / 1.998 | 67.695 / 68.290 | 14.135 | -1.33 |

map1与后四图位置、field和速度不同，只保留为背景；时间结论只来自map2–5。

## Sensitivity QC

| map | 30 nm contact | 40 nm contact | 50 nm contact | retained / 64 |
|---:|---:|---:|---:|---:|
| 1 | 63.448 | 63.686 | 65.161 | 61 / 64 |
| 2 | 70.026 | 71.529 | 73.133 | 59 / 64 |
| 3 | 69.359 | 70.877 | 72.563 | 64 / 64 |
| 4 | 67.917 | 70.176 | 73.199 | 62 / 64 |
| 5 | 66.392 | 67.900 | 70.985 | 63 / 64 |

## Measured force slices

以下为far-linear、64-pixel median，单位pN。

| D (nm) | map1 | map2 | map3 | map4 | map5 |
|---:|---:|---:|---:|---:|---:|
| 20 | 137.44 | 392.35 | 343.32 | 259.07 | 144.44 |
| 50 | 42.45 | 135.68 | 62.34 | 38.02 | 26.04 |
| 100 | 1.24 | -8.66 | -15.23 | -12.74 | -8.62 |
| 200 | -11.07 | -38.50 | 6.65 | -4.51 | -4.28 |

## 四时点pixel-level时间趋势

对每个physical pixel分别用map2–5四个点拟合force对时间的直线，再把slope乘总时间 `78.66 min`。这比直接拟合四个map median更能保留same-pixel结构；但每个pixel只有四个时间点，仍只能视为linear trend diagnostic。

### Primary: far-linear baseline

| D (nm) | map2 → map3 → map4 → map5 medians (pN) | fitted Map5−Map2 median [IQR] (pN) | negative slopes | strict decrease | strict increase | pixel p | row-block p |
|---:|:---|:---|---:|---:|---:|---:|---:|
| 20 | 392.35 → 343.32 → 259.07 → 144.44 | -213.41 [-333.99, -112.69] | 93.8% | 29.7% | 0.0% | 7.49e-12 | 0.00781 |
| 50 | 135.68 → 62.34 → 38.02 → 26.04 | -99.36 [-186.46, -23.83] | 79.7% | 14.1% | 0.0% | 1.48e-08 | 0.00781 |
| 100 | -8.66 → -15.23 → -12.74 → -8.62 | -10.47 [-88.03, +58.70] | 54.7% | 6.2% | 4.7% | 0.288 | 0.641 |
| 200 | -38.50 → 6.65 → -4.51 → -4.28 | +35.52 [-44.64, +89.01] | 37.5% | 6.2% | 6.2% | 0.11 | 0.195 |

### Baseline sensitivity of full-interval trend

| D (nm) | fitted change linear (pN) | fitted change constant (pN) | constant row-block p |
|---:|---:|---:|---:|
| 20 | -213.41 | -233.72 | 0.00781 |
| 50 | -99.36 | -114.45 | 0.00781 |
| 100 | -10.47 | -27.70 | 0.00781 |
| 200 | +35.52 | +13.03 | 0.195 |

## 相邻与全间隔paired changes

| pair | D (nm) | paired Δ median [IQR] (pN) | Wilcoxon p | row-block p | spatial Spearman |
|:---|---:|:---|---:|---:|---:|
| map2→map3 | 20 | -32.36 [-136.23, +73.47] | 0.276 | 0.945 | +0.796 |
| map2→map3 | 50 | -28.58 [-106.27, +66.05] | 0.0525 | 0.945 | +0.527 |
| map2→map3 | 100 | -12.11 [-112.57, +73.01] | 0.407 | 0.844 | +0.407 |
| map2→map3 | 200 | +16.02 [-43.08, +91.66] | 0.0881 | 0.312 | +0.377 |
| map2→map4 | 20 | -107.71 [-235.95, +9.89] | 5.26e-06 | 0.0156 | +0.759 |
| map2→map4 | 50 | -60.61 [-163.08, -1.53] | 3.38e-05 | 0.109 | +0.254 |
| map2→map4 | 100 | -21.50 [-108.67, +59.32] | 0.256 | 0.641 | +0.102 |
| map2→map4 | 200 | +25.26 [-31.84, +94.27] | 0.0203 | 0.148 | +0.006 |
| map2→map5 | 20 | -201.53 [-339.86, -106.56] | 1.31e-11 | 0.00781 | +0.682 |
| map2→map5 | 50 | -97.24 [-177.42, -29.80] | 1.08e-08 | 0.00781 | +0.062 |
| map2→map5 | 100 | -6.79 [-99.30, +46.91] | 0.221 | 0.461 | -0.137 |
| map2→map5 | 200 | +31.14 [-47.63, +82.90] | 0.164 | 0.25 | -0.150 |
| map3→map4 | 20 | -96.35 [-165.18, -28.43] | 1.21e-08 | 0.00781 | +0.898 |
| map3→map4 | 50 | -36.26 [-108.04, +21.26] | 0.000867 | 0.109 | +0.742 |
| map3→map4 | 100 | +10.75 [-58.93, +55.49] | 0.984 | 0.641 | +0.752 |
| map3→map4 | 200 | +10.08 [-49.43, +59.25] | 0.508 | 0.461 | +0.625 |
| map3→map5 | 20 | -202.27 [-298.22, -90.87] | 2.98e-11 | 0.00781 | +0.739 |
| map3→map5 | 50 | -60.63 [-171.17, +11.71] | 9.15e-05 | 0.0391 | +0.289 |
| map3→map5 | 100 | -4.57 [-105.52, +90.42] | 0.669 | 0.844 | +0.112 |
| map3→map5 | 200 | +3.98 [-91.86, +59.53] | 0.579 | 0.844 | +0.054 |
| map4→map5 | 20 | -110.22 [-182.90, -40.59] | 2.76e-09 | 0.00781 | +0.800 |
| map4→map5 | 50 | -22.14 [-95.06, +18.75] | 0.00386 | 0.148 | +0.376 |
| map4→map5 | 100 | +0.82 [-77.47, +51.82] | 0.616 | 1 | +0.254 |
| map4→map5 | 200 | -11.68 [-59.82, +36.37] | 0.199 | 0.844 | +0.186 |

## Interpretation of the four-map sequence

- **20 nm存在最清楚的time/history dependence。** map median为 `392.4 → 343.3 → 259.1 → 144.4 pN`；same-pixel full-interval change为 `-213.4 pN`（linear）或 `-233.7 pN`（constant）。linear和constant的row-block p分别为 `0.00781`、`0.00781`。
- **50 nm的精确量级更受baseline与空间相关影响。** full-interval paired trend为 `-99.4 pN`（linear）和 `-114.4 pN`（constant）；negative slope pixel占 `79.7%`。linear row-block p为 `0.00781`，constant为 `0.00781`。因此不能把任一单一拟合值当作所有pixel的唯一校正。
- **100 nm结论依赖zero-force定义。** linear给出 `-10.5 pN`且row-block不显著；constant给出 `-27.7 pN`且row-block p=`0.00781`。200 nm主要受far-field baseline控制，不宜作为absolute interaction-force time trend。
- **这不是简单统一的linear/exponential relaxation。** 相邻map在20 nm的paired median依次变化 `-32.4`、`-96.3`、`-110.2 pN`；50 nm依次为 `-28.6`、`-36.3`、`-22.1 pN`。四个时点仍不足以唯一确定time law。
- speed、load和sensitivity不是主要解释。map2→5 paired gap-speed只改变 `+0.0068 µm/s`；按纯水no-slip估算，在20/50 nm只对应 `0.117/0.047 pN`。local InvOLS paired median变化 `-1.238 nm/V`，terminal load变化 `+0.0106 nN`。
- 最稳妥的结论是：20–50 nm force随约78.7 min实验history变化，但具有空间异质性，不能用单一offset或单一time slope校正所有pixel。


## Within-map acquisition-order trend

| map | D (nm) | Theil–Sen endpoint change (pN) |
|---:|---:|---:|
| 1 | 20 | -118.14 |
| 1 | 50 | -18.33 |
| 1 | 100 | +25.02 |
| 1 | 200 | +29.63 |
| 2 | 20 | +249.64 |
| 2 | 50 | +54.99 |
| 2 | 100 | +19.77 |
| 2 | 200 | +44.29 |
| 3 | 20 | +217.26 |
| 3 | 50 | +128.85 |
| 3 | 100 | +122.76 |
| 3 | 200 | +84.36 |
| 4 | 20 | +161.46 |
| 4 | 50 | +101.89 |
| 4 | 100 | +93.22 |
| 4 | 200 | +46.87 |
| 5 | 20 | +99.32 |
| 5 | 50 | +15.56 |
| 5 | 100 | -0.67 |
| 5 | 200 | -7.89 |

### Map2–5 at 50 nm: separate 64-point fits

Theil–Sen为主拟合，OLS只作对照。斜率的sample单位是一条force curve；1→64拟合变化为斜率×63。它描述acquisition-order association，不能单独解释为纯时间drift，因为raster位置与时间混杂。

| map | Theil–Sen slope (pN/sample) | 95% slope CI | fitted 1→64 change (pN) | OLS change (pN) | OLS R² |
|---:|---:|:---|---:|---:|---:|
| 2 | +0.873 | [-0.450, +2.290] | +54.99 | +83.48 | 0.049 |
| 3 | +2.045 | [-0.067, +3.885] | +128.85 | +132.40 | 0.083 |
| 4 | +1.617 | [-0.068, +2.991] | +101.89 | +80.01 | 0.065 |
| 5 | +0.247 | [-0.600, +1.077] | +15.56 | +3.05 | 0.000 |

## QC pair comparisons

以下只列map2→map5全间隔QC；完整六对结果在CSV中。

| metric | map2 | map5 | paired Δ median [IQR] | Wilcoxon p |
|:---|---:|---:|:---|---:|
| local_contact_InvOLS (nm/V) | 71.53 | 67.9 | -1.238 [-3.298, +2.024] | 0.0516 |
| far_slope (pN/100 nm) | -2.108 | -1.329 | +0.01354 [-2.97, +3.98] | 0.683 |
| far_noise (pN) | 19.68 | 11.89 | -7.484 [-11.3, -3.309] | 1.09e-11 |
| terminal_load (nN) | 14.13 | 14.13 | +0.01063 [-0.1838, +0.1089] | 0.517 |
| gap_speed (um/s) | 1.992 | 1.998 | +0.006825 [+0.004843, +0.008994] | 3.53e-12 |
| retract_pull_off_force (nN) | -4.226 | -4.061 | +0.08539 [-5.477, +2.342] | 0.598 |
| retract_detachment_travel (nm) | 21.22 | 20.24 | -2.377 [-7.507, +11.26] | 0.913 |

## Calibration fingerprint

| candidate | stored InvOLS reference (nm/V) | stored k reference (N/m) | log-ratio distance |
|:---|---:|---:|---:|
| D4 | 70.112 | 0.177639 | 0.6888 |
| D5 | 75.530 | 0.348252 | 0.0363 |
| D6 | 69.478 | 0.191763 | 0.6151 |

## Claim boundary

- map2–5是same-location、same-speed四时点序列，可判断repeatability/chronology是否存在，但四个时点仍不足以唯一确定linear、exponential或其他time law。
- 64个pixel具有spatial/serial correlation；pixel p值偏乐观，8个physical-row block结果作为保守敏感性检查。
- far-linear与far-constant都保留；两者不一致时不能把变化解释为absolute interaction-force变化。
- Force scale仍以D5 metadata inference为条件；若实际cantilever不同，force需按正确calibrated k整体rescale。
- map1不参与时间拟合，也不能与map2唯一识别velocity effect。

## Outputs

- `map2_map5_pairwise_time_comparison.csv`, `map2_map5_time_trends.csv`, `map2_map5_paired_QC.csv`.
- `map_inventory_QC.csv`, `contact_sensitivity_fits.csv`, `pixel_QC.csv`, `map_force_curves.csv`, `pixel_force_curves.npz`.
- `map2_map5_force_50nm_acquisition_fits.csv` and the four separate `figures/map[2-5]_force_50nm_vs_acquisition_order.png` plots.
- `map2_map5_absolute_force_time_statistics.csv` and `figures/absolute_force_vs_time_mean_std_variance.png`; statistics use per-pixel `|F|`, with sample SD/variance (`ddof=1`).
- `figures/map2_map5_time_force_maps.png`, `figures/map2_map5_time_statistics.png`及五图总览/QC。
- `provenance.json`, `artifact_manifest.sha256`.

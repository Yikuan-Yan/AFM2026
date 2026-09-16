# 16-09-26 low-speed palindrome AFM analysis

## 主要结论

18 张 8x8 force maps 的 1152 条 approach 中，**1137 条可用于 force-distance 重建**；8 条是小于 50 nm 的明显 premature trigger，1 条量程不足，6 条虽有较长记录但没有满足振幅和线性判据的 terminal hard contact。JPK acquisition flags 对这些点全部未报警，因此必须由曲线级 QC 识别并保留为 NaN，不能插值补点。
18 张 map 的 point index 0 全部通过相同判据；因此黑色 height-map 像素、行首身份或 point index 本身都不是坏点标签。

第一张 map 的中心相对其余 17 张偏移约 **0.485 um**。因此 block 1 的 0.1 um/s palindrome mate 只能比较 map median，不能做同像素配对；其余 palindrome pairs 的网格中心一致。

## 曲线处理

- branch loader 以 `minimum_points=1` 解码，确保短段仍带原始 point index 进入 QC。主分析只用 Extend；Retract 固定为 2 um/s，不与 approach 速度混合。
- 1152 条 Retract 全部为 500 点，travel span 为 496.4-497.2 nm。坏点的完整 Retract 说明 ZIP/像素记录仍在，故应解释为 approach-side premature trigger/contact failure；不能用 2 um/s Retract 替代 approach force。
- baseline 是每条曲线在 scanner travel 100-250 nm 的 raw-vDeflection 中位数，只减常数，不减线性斜率。它定义 finite-window force gauge，不等于无限远零力。
- contact 用 terminal 40 nm 物理窗口；force-use 判据同时要求 travel、contact amplitude、InvOLS 范围和 robust-line R2。固定样本数没有被使用。
- batch contact InvOLS = **58.244 nm/V**（18 个 map median 再取 median），embedded InvOLS = 84.506 nm/V。embedded InvOLS 不满足本批 hard-contact slope，因此未用于位移换算。
- contact-window sensitivity：30/40/50 nm 分别得到 58.147/58.244/58.394 nm/V；主值没有依赖单一固定样本数。
- spring constant 暂用所有文件一致的 embedded `k=0.174367775 N/m`，得到条件 force scale **10.156 nN/V**。没有同日独立 thermal calibration，因此绝对力标度仍带 calibration 条件。

## 坏点清单

| map | scan | time | speed (um/s) | point | pixel (r,c) | state | points | travel (nm) | contact amplitude (V) |
|---:|---:|---|---:|---:|---|---|---:|---:|---:|
| 1 | 9 | 10:07:02 | 0.1 | 51 | (6,3) | no_valid_contact | 8821 | 441.0 | 0.064 |
| 2 | 10 | 10:12:05 | 0.2 | 16 | (2,0) | hard_fail | 59 | 5.6 | nan |
| 3 | 11 | 10:14:38 | 0.4 | 17 | (2,1) | hard_fail | 19 | 3.2 | nan |
| 4 | 12 | 10:16:27 | 0.4 | 16 | (2,0) | no_valid_contact | 2227 | 444.8 | 0.293 |
| 5 | 13 | 10:19:20 | 0.2 | 16 | (2,0) | no_valid_contact | 4314 | 431.2 | 0.164 |
| 6 | 14 | 10:24:14 | 0.1 | 1 | (0,1) | hard_fail | 8 | 0.3 | nan |
| 7 | 15 | 10:28:49 | 0.2 | 1 | (0,1) | hard_fail | 23 | 2.0 | nan |
| 7 | 15 | 10:28:49 | 0.2 | 18 | (2,2) | hard_fail | 25 | 2.3 | nan |
| 8 | 16 | 10:31:34 | 0.4 | 18 | (2,2) | no_valid_contact | 2067 | 413.0 | 0.087 |
| 9 | 17 | 10:35:52 | 0.1 | 18 | (2,2) | hard_fail | 14 | 0.6 | nan |
| 11 | 19 | 10:45:58 | 0.4 | 29 | (3,2) | no_valid_contact | 2356 | 470.7 | 0.367 |
| 12 | 20 | 10:48:38 | 0.2 | 17 | (2,1) | hard_fail | 23 | 2.0 | nan |
| 13 | 21 | 10:51:18 | 0.4 | 32 | (4,0) | hard_fail | 8 | 0.9 | nan |
| 14 | 22 | 10:55:38 | 0.1 | 31 | (3,0) | no_valid_contact | 9543 | 477.2 | 0.349 |
| 17 | 25 | 11:09:52 | 0.1 | 30 | (3,1) | insufficient_range | 2201 | 110.0 | -0.000 |

## 低速差异是否小于实验随机性

下表给 palindrome-symmetrized speed contrast，并减去同一 contrast 在 250 nm 的值；这样不受每条曲线常数 baseline gauge 影响。n=3 是三个 block，不是像素数。`mean/SE < 1` 只能说明当前实验未分辨出差异，不能证明等效；正式 equivalence 还需要预先规定可接受 margin。

| D (nm) | interval (um/s) | mean double difference (pN) | block SD (pN) | SE (pN) | |mean|/SE |
|---:|---|---:|---:|---:|---:|
| 20 | 0.1->0.2 | +166.70 | 11.84 | 6.84 | 24.38 |
| 20 | 0.2->0.4 | +263.93 | 16.66 | 9.62 | 27.45 |
| 50 | 0.1->0.2 | +78.87 | 7.36 | 4.25 | 18.56 |
| 50 | 0.2->0.4 | +125.96 | 2.82 | 1.63 | 77.47 |
| 100 | 0.1->0.2 | +41.26 | 4.24 | 2.45 | 16.86 |
| 100 | 0.2->0.4 | +52.57 | 4.47 | 2.58 | 20.36 |
| 200 | 0.1->0.2 | +11.01 | 8.39 | 4.84 | 2.27 |
| 200 | 0.2->0.4 | +6.26 | 7.62 | 4.40 | 1.42 |

在 20/50/100 nm，0.1->0.2 um/s 的三个 block 全部同号，且 |mean|/SE 分别为 24.4/18.6/16.9。因此本批数据明确不支持 0.1 与 0.2 um/s 已进入同一 plateau；到 200 nm 时 95% CI 仍跨零，只能说该远距离差异未分辨。

由重建的 `D(t)` 直接求得：20 nm 处 `U_gap/nominal` 的 map-median 范围为 **0.813-0.886**，200 nm 处接近 1。以下模型因此同时列 nominal speed 与实际 `U_gap(D)` 分支。

## F(v,D) 模型检查

每个 6-map block 先平均同速 palindrome mates，再对 distance x speed 矩阵 double-center。这个操作严格去掉任意 distance-only surface-force 曲线和任意 speed-only 常数 offset；剩余项才用于比较 `v/D` 与经验式。它仍可能包含 speed-dependent instrument response 或非线性 history，不能仅凭拟合命名为 hydrodynamics。

| pooled model, 20-200 nm | K | alpha | D0 (nm) | RMSE (pN) | normalized RMSE | Jacobian condition |
|---|---:|---:|---:|---:|---:|---:|
| additive_null | 0 | 1.000 | 0.00 | 42.54 | 1.000 | nan |
| theoretical_no_slip | 3.15e+05 | 1.000 | 0.00 | 361.38 | 8.494 | nan |
| fitted_v_over_D | 3.27e+04 | 1.000 | 0.00 | 7.66 | 0.180 | 1 |
| shifted_v_over_D | 6.16e+04 | 1.000 | 17.82 | 3.84 | 0.090 | 7.88e+03 |
| velocity_power_over_D | 3.05e+04 | 0.776 | 0.00 | 7.48 | 0.176 | 1.05e+04 |
| empirical_power_shift | 5.75e+04 | 0.784 | 17.76 | 3.48 | 0.082 | 6.95e+04 |
| actual_gap_theoretical_no_slip | 3.15e+05 | 1.000 | 0.00 | 285.90 | 6.720 | nan |
| actual_gap_fitted | 4.05e+04 | 1.000 | 0.00 | 4.99 | 0.117 | 1 |
| actual_gap_shifted | 5.28e+04 | 1.000 | 6.39 | 3.70 | 0.087 | 7.15e+03 |

在 pooled in-sample 比较中，`v/(D+D0)` 相对 fitted `v/D` 把 RMSE 降低 **49.9%**；完整经验式相对 fitted `v/D` 降低 **54.5%**，但相对已经带 D0 的线性速度式只再降低 **9.3%**。三个 block 的 double-centered shape 两两相关范围为 **0.993 到 0.995**。block-wise alpha 为 0.647-0.938，而 D0 为 17.32-18.49 nm。经验式多两个自由参数且仅有三个速度，Jacobian condition 也很大；因此较低 RMSE 不能把 alpha 解释成已识别的 shear-thinning exponent，也不能把 D0 直接解释成 slip length。

用实际 `U_gap(D)` 后，最佳 shifted model 的 normalized RMSE 为 0.087。固定 bulk no-slip coefficient 仍显著过预测，说明这批数据支持可重复的 velocity-distance coupling，但不验证以 bulk viscosity、prior radius、nominal no-slip geometry 组成的绝对 prefactor。

固定 no-slip 系数仅作为 campaign-context comparator：假定 99.7 wt% glycerol、T=25.6 C、R=4.546849 um。这些量不在 JPK header 中，主 empirical/fitted comparisons 不依赖固定理论系数。

## 产物

- `curve_qc.csv`: 1152 条 approach 的长度、量程、baseline、contact、retract inventory 与最终状态。
- `map_inventory_and_summary.csv`: acquisition order、grid shift、校准、QC 数量及 20/50/100/200 nm force。
- `contact_calibration_maps.csv`: 等 map 权重的 contact InvOLS。
- `curve_force_bins.npz` / `map_force_by_separation.csv`: NaN-preserving force、actual gap-speed pixel matrices 与 map summaries。
- `palindrome_pair_force.csv`: 同速 early/late、pair mean 与 history half-difference。
- `speed_interval_detail.csv` / `speed_interval_summary.csv`: 低速 interval 的 block-level 结果。
- `hydrodynamic_double_centered.csv`, `hydrodynamic_model_fits.csv`, `block_reproducibility.csv`: `F(v,D)` 检查。
- `figures/`: QC、calibration/baseline、force-time、speed intervals 与 model comparison。
- `provenance.json`, `input_manifest.sha256`, `artifact_manifest.sha256`: 输入 hash、定义、软件和输出 hash。
- `NUMERICS_AUDIT.md`: units、shape、NaN、window/range sensitivity、optimizer 与 synthetic closure 检查。

本脚本没有进行 PB/Debye/surface-potential 拟合。先确认低速速度项是否可重现、校准是否独立成立，再决定能否把某个 subtraction 后的残差解释为 equilibrium surface force。

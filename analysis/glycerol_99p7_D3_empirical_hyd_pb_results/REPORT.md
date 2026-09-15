# 12-09-26 - empirical hydrodynamic subtraction and apparent-PB refit

## 结论

**经验式不能被称为完全、唯一地去除了 hydrodynamics。条件于该经验模型，可以逐点减去完整的 `K v^alpha/(D+D0)`；以下 Debye length 和 surface potential 仍只能称为 empirical-hyd-corrected、model-conditioned apparent parameters。**

主分支把四个 palindrome blocks 共享拟合为 `K v^alpha/(D+D0)`：`K=112407.0 pN nm/(um/s)^alpha`、`alpha=0.6502`、`D0=37.46 nm`。输出表实际保存 `F_full=F-K v^alpha/(D+D0)`。由于输入曲线此前已减过有限窗口 constant baseline，完整减法会留下一个 speed-dependent constant；PB 优化前精确加回该模型在 250 nm 的值，拟合后再从 offset 中减回。这个常数变换不改变 lambda_D、|psi|、残差或 R2。

经验 hyd 项来自同一批 force curves，因此 speed dependence 变弱是部分 in-sample consequence，不能当作对 hydrodynamic mechanism 的独立验证。主分支使用所有 block 共用参数，以免把 block-to-block 时间变化直接吸收到 block-specific hyd amplitude；后者仅作敏感性检查。

## 参数趋势比较

global speed rho 使用 quality_state != unusable 的 map；same-speed time rho 先在四个速度内分别计算再取中位数。每小时变化来自 `log(parameter) ~ elapsed_hours + speed categorical effects`，HC3 interval 只作描述。

| parameter | global rho(speed), original -> corrected | median same-speed rho(time), original -> corrected | speed-adjusted change/h, original -> corrected | corrected palindrome rho(speed), block median |
|---|---:|---:|---:|---:|
| Apparent Debye length | 0.968 -> -0.304 | -0.238 -> -0.890 | -5.0% -> -28.8% | -1.000 [-1.000, 0.500] (n=3) |
| Apparent surface potential | 0.968 -> 0.414 | -0.798 -> -0.502 | -8.0% -> -8.4% | 1.000 [0.500, 1.000] (n=3) |

完整的整体关联和 HC3 描述区间见 `parameter_time_speed_associations.csv`。速度可视化优先使用同 block 的 symmetric pair means，从而压低一阶 acquisition-time drift；校正后并非每个 block 的四速度 pairs 都通过 QC，所以汇总 rho 使用至少三个可用速度的 block。

Hyd 参数作用域与拟合范围的敏感性不是小修正，尤其会改变 lambda_D 的时间方向：

| corrected branch | rho_speed lambda | adjusted lambda change/h | rho_speed potential | adjusted potential change/h |
|---|---:|---:|---:|---:|
| Shared full empirical subtraction | -0.304 | -28.8% | +0.414 | -8.4% |
| Shared empirical 20-250 sensitivity | -0.203 | -39.7% | +0.571 | -8.2% |
| Block-specific sensitivity | +0.103 | +74.3% | +0.187 | -16.2% |

主分支逐速度时间趋势如下；first/last 指该速度首末两个 usable fits，不保证对应整个实验的共同端点：

| speed (um/s) | n lambda | rho_time lambda | lambda first -> last (nm) | n potential | rho_time potential | potential first -> last (mV) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 8 | -0.262 | 14.55 -> 11.31 | 8 | -0.643 | 39.46 -> 33.79 |
| 0.3 | 8 | -0.881 | 12.60 -> 7.62 | 8 | -0.405 | 45.68 -> 38.97 |
| 0.9 | 5 | -0.900 | 14.65 -> 8.44 | 5 | -0.600 | 49.25 -> 42.06 |
| 2.7 | 4 | -1.000 | 18.92 -> 8.03 | 4 | 0.800 | 36.32 -> 46.71 |

## Fit quality

| variant | pass | weak | unusable | median usable R2 | median usable RMSE (pN) |
|---|---:|---:|---:|---:|---:|
| Original | 24 | 7 | 1 | 0.979 | 48.3 |
| Shared full empirical subtraction | 9 | 16 | 7 | 0.937 | 17.4 |
| Shared empirical 20-250 sensitivity | 12 | 16 | 4 | 0.950 | 16.4 |
| Block-specific sensitivity | 5 | 22 | 5 | 0.919 | 18.7 |

所有 corrected branches 都复用原始 nonlinear PB Derjaguin + vdW + constant offset、20-250 nm fit window、spatial-MAD weights、四起点 optimizer 和相同 quality flags。fit quality 只说明该投影的数值状态，不证明 corrected curve 是 equilibrium PB force。

## 可解释范围

1. corrected force curves 可以用来判断：移除经验 interaction 后，原来很强的 apparent-parameter speed ordering 还剩多少，以及同速度随时间的趋势是否保留。
2. 不能把 corrected lambda_D 当作 bulk Debye length、corrected |psi| 当作已识别 equilibrium surface potential；empirical subtraction、contact zero、epsilon_r=42.5、R、Hamaker constant 和 force scale 都是条件。
3. block-specific subtraction 若比 shared subtraction 更大幅消除时间趋势，说明它同时吸收了 history；这不是更完全的 hydrodynamic correction。
4. 真正的 physical subtraction 仍需 raw `D(t)` 和 `U_gap(t)=-dD/dt`，并用独立 viscosity、temperature、blank/control 和 approach-retract sign reversal 检查。当前 raw JPK 不在本地 checkout。

## 文件

- `map_force_empirical_hyd_corrected.csv`: 三个 branches 的完整模型减法曲线、250 nm gauge-centered PB 输入曲线及两种 correction 分量。
- `map_pb_fits_comparison.csv`: 原始与三个 corrected branches 的全部 PB fits、full-curve offset 和 QC。
- `empirical_hyd_parameter_sets.csv`: shared 与 block-specific 经验参数。
- `parameter_time_speed_associations.csv`: time/speed rank association 与 speed-adjusted log-time trend。
- `same_speed_parameter_time_trends.csv`: 每个速度的 usable-count、首末值、rho 和 log-time slope。
- `palindrome_parameter_pair_means.csv`, `palindrome_parameter_speed_summary.csv`: block 内速度对照。
- `low_speed_parameter_pair_comparison.csv`: 0.3-0.1 um/s 的四-block PB 参数 paired SE 与描述区间。
- `low_speed_force_shape_comparison.csv`: 消去 250 nm constant gauge 后的逐距离低速 force-shape 对照。
- `fit_quality_summary.csv`: pass/weak/unusable 与 residual summary。
- `figures/force_curves_before_after_empirical_subtraction.*`: 原始、完整模型减法、PB gauge-centered 三联图。
- `figures/pb_parameters_time_original_vs_empirical_corrected.*`: Debye/potential 时间趋势。
- `figures/pb_parameters_time_empirical_corrected_detail.*`: corrected-only 同速度时间趋势放大图。
- `figures/pb_parameters_speed_palindrome_original_vs_corrected.*`: palindrome-controlled 速度比较。
- `provenance.json`, `artifact_manifest.sha256`: 输入 hashes、模型定义、软件和产物校验。

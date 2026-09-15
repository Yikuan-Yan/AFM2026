# 12-09-26 99.7 wt% glycerol：经验 hydrodynamic correction 与 PB 重拟合

## 核心结论

本次分析确认：经过逐曲线 constant baseline 校正后，force curves 中仍存在很强、跨
block 可复现的 velocity-distance interaction。经验式

```text
F_emp(D,v) = K v^alpha / (D + D0)
```

对该 interaction 的描述明显优于严格的 nominal-speed `v/D`。在这个经验模型成立的
条件下，可以逐点减去其完整拟合值并重新进行 PB 拟合；但这不等于真实物理
hydrodynamic force 已经被唯一识别或完全消除。

使用四个 blocks 共享的经验参数校正后，原来 apparent PB parameters 随速度强烈单调
增大的排序大幅塌缩，但没有完全消失。移除大幅度的 hydrodynamic-like signal 后，剩余
曲线对两个 PB shape parameters 的约束反而变弱。校正后的 apparent surface potential
随时间下降的方向比 apparent Debye length 更稳定；后者在 block-specific hyd correction
下发生符号反转。

当前数据也不能证明 `v=0.1` 和 `v=0.3 um/s` 已进入 equilibrium plateau。四个
palindrome blocks 中，这两个速度的 corrected PB parameters 都保留同方向差值，而且
差值大于 block-level descriptive paired standard error。

## 1. 研究问题与证据范围

本分析依次回答四个问题：

1. 如果 hydrodynamic force 同时依赖速度和间距，逐曲线 constant baseline 是否能将其
   完全消除？
2. 观察到的 velocity-distance interaction 是否符合 classical nominal-speed `v/D`，还是
   需要更灵活的经验式？
3. 减去经验项后，apparent Debye length 和 apparent surface-potential magnitude 还保留
   怎样的时间趋势与速度相关性？
4. `0.1` 与 `0.3 um/s` 是否已经慢到两者差异小于实验重复性的 standard error？

当前本地 checkout 保存了 derived force、chronology、calibration 和 palindrome tables，
但没有原始 JPK archives。因此本包不能从 raw trajectories 重建 instantaneous gap-closing
speed

```text
U_gap(t) = -dD/dt,
```

也不能独立检验 viscosity、approach/retract sign reversal 或 scanner speed 与实际 gap
speed 的对应关系。所有物理结论都受此证据边界约束。

## 2. 实验组织与重复单位

数据包含顺序采集的 32 张 force maps，nominal approach speeds 为 `0.1, 0.3, 0.9,
2.7 um/s`。它们组成四个 8-map palindrome blocks；每个速度在每个 block 内出现两次，
并位于关于 block midpoint 对称的位置。

为了压低一阶 acquisition-time drift，与速度有关的比较先平均同一 block 内两张
same-speed symmetric maps。描述性重复单位是四个 blocks，而不是每张 map 内的 64 个
pixels，也不是同一 force curve 上相互关联的 5 nm distance bins。

Palindrome pairing 可以降低一部分 time-speed confounding，但不能把 sequential experiment
变成四个独立随机重复。因此本文中的 standard error 和 t interval 只用于描述效应尺度，
不作为严格独立重复下的显著性结论。

所有输入路径与 SHA-256 记录在 [provenance.json](provenance.json)。

## 3. 为什么 constant baseline 不能消除全部 drainage

一条已经过 baseline 校正的 approach curve 可以写成

```text
F_obs(D,v,t) = F_surface(D,t)
             + F_hyd(D,U_gap,t)
             - B_curve(v,t)
             + noise.
```

当前 dynamic baseline 对每条 approach curve 只减去一个常数。它能够去掉 detector offset，
也会去掉 finite baseline window 内 hydrodynamic contribution 的平均值；但它不能去掉该
contribution 随 separation 的变化。

如果残余 hydrodynamics 仅为 `F_hyd=F(v)`，speed-conditioned constant subtraction 确实
可以将它消除。如果 `F_hyd=F(v,D)`，baseline 后通常仍保留

```text
F_hyd(D,v) - <F_hyd(Db,v)>baseline_window.
```

因此，baseline 后仍观察到可复现的 velocity-distance interaction，排除了纯 speed-only
residual。这一事实支持 hydrodynamic-like distance dependence，却不自动验证 classical
no-slip drainage。

## 4. Baseline-invariant interaction 检验

令 `F_b(D,v)` 表示 block `b` 内 same-speed palindrome pair 的平均力。对每个
distance-by-speed matrix 做 double centering：

```text
I_b(D,v) = F_b(D,v)
           - <F_b(D,v)>v
           - <F_b(D,v)>D
           + <F_b(D,v)>D,v.
```

该操作精确消去任意 distance-only contribution `A(D)`、任意 speed-only contribution
`B(v)` 以及 block constant，只保留 velocity-distance interaction。

四个 blocks 中，第一 singular component 解释 interaction squared amplitude 的
`99.90%-99.95%`；block interaction shapes 的两两 correlation 为 `0.9972-0.9994`。
因此，该 interaction 高度可复现，并且近似 separable。

在 `20-200 nm` 范围内：

- 自由拟合 amplitude 的严格 `K v/D`，block-median normalized RMSE 为 `0.277`；
- 加入 distance shift、但保持 velocity linearity，降至 `0.124`；
- 同时允许 `D0` 和 nominal-speed exponent `alpha`，降至 `0.036`。

完整 interaction derivation、block results 和 fit-window sensitivity 见
[hydrodynamic-scaling report](../glycerol_99p7_D3_hydrodynamic_scaling_results/REPORT.md)。

### 4.1 为什么旧的 v=0 图中有两个截距

旧 zero-speed extrapolation 图中的两个 marker 不是同一个模型中的两个物理截距，也不
代表两种 surface forces：

- 实心菱形：使用 `0.1, 0.3, 0.9, 2.7 um/s` 四个速度做 linear extrapolation 得到的
  `v=0` intercept；
- 空心方形：只使用 `v<=0.9 um/s` 三个速度得到的 sensitivity intercept。

两者采用不同的 speed fit range。删除 `2.7 um/s` 后截距明显变化，说明 nominal-speed
linearity 不稳定，尤其是高速点不能被视为同一条简单直线上的普通重复。两个截距的差异
正是 model-range sensitivity，而不是发现了第二个 equilibrium force。

这也是后续使用 `K v^alpha/(D+D0)` 描述 interaction、并把线性 `v=0` 截距降级为
protocol-dependent apparent quantity 的原因。旧结果见
[decay/v0 report](../glycerol_99p7_D3_decay_v0_results/REPORT.md)。

## 5. 经验模型与 primary 参数选择

Primary correction 使用四个 palindrome blocks 共享的一组参数：

```text
K     = 112407.0 pN nm / (um/s)^alpha
alpha = 0.6502
D0    = 37.46 nm
```

即

```text
F_emp(D,v) = 112407.0 v^0.6502 / (D + 37.46) pN,
```

其中 `D` 使用 nm，nominal `v` 使用 um/s。

Primary branch 使用 shared parameters，是为了避免每个 sequential block 各自调整
amplitude 和 shape，从而直接把待研究的 block-to-block history 吸收到 hyd correction。
同时保留两条 sensitivity branches：

- shared parameters，fit range 扩展至 `20-250 nm`；
- 每个 block 独立拟合 `20-200 nm` 参数。

`alpha<1` 不能直接解释为 glycerol shear-thinning exponent，`D0` 也不能直接命名为 slip
length。它们可能同时混合 nominal-to-gap speed mapping、cantilever compliance、contact-zero
error、roughness、slip 和 feedback response。

Double centering 会消去任意 `B(v)`，所以数据只直接识别 distance-varying interaction；
将 `F_emp` 定义为 `D` 趋于无穷时归零，是重建完整经验项时采用的 gauge convention，
不是由本数据独立测得的 absolute hydrodynamic zero。

## 6. 完整经验减法与 constant-offset gauge

输出表实际保存完整经验模型减法：

```text
F_full_corrected(D,v) = F_original(D,v) - K v^alpha/(D+D0).
```

由于 `F_original` 已经减过 finite-window constant baseline，这个完整减法会留下一个
speed-dependent constant zero。它不含 distance-shape information，应由 PB fit 中的 nuisance
offset 处理。

为了沿用既有 offset bounds 并保持数值条件，optimizer 接收

```text
F_PB_input(D,v) = F_full_corrected(D,v)
                  + K v^alpha/(250 nm + D0).
```

拟合完成后，再从 fitted offset 中减回同一常数：

```text
b_full_curve = b_PB_input - K v^alpha/(250 nm + D0).
```

这是精确的 reparameterization，不会改变 `lambda_D`、surface-potential magnitude、residual
或 `R2`。数值审计得到的最大 residual-equivalence error 为 `1.82e-12 pN`。

[subtraction and gauge figure](figures/force_curves_before_after_empirical_subtraction.png)
依次显示 original baseline-corrected curves、完整模型减法后的 curves、以及为 PB fit 加回
250 nm constant 的 gauge-centred curves。两种 corrected representations 都保存在
[map_force_empirical_hyd_corrected.csv](map_force_empirical_hyd_corrected.csv)。

因此，本包中的“完整减法”应读作：**完整减去假定经验式的拟合值**。它不能升级为
“完整去除真实 hydrodynamics”。

## 7. PB 重拟合定义

Original 与所有 corrected branches 使用同一模型和 QC implementation：

- equal-potential nonlinear PB force，Derjaguin geometry；
- sphere-plane van der Waals term；
- fitted constant nuisance offset；
- fit interval `20-250 nm`；
- spatial-MAD weights，最低 `5 pN`；
- 四个 deterministic optimizer starts；
- assumed dielectric constant `epsilon_r=42.5`；
- probe radius `R=4.546849 um`；
- Hamaker constant `A_H=2.4e-21 J`；
- temperature `T=25.6 degC`。

自由参数只有 `lambda_D`、potential magnitude 和 constant offset。同表面 force model 不能
识别 potential sign。因此输出必须称为 `apparent lambda_D` 和 `apparent surface-potential
magnitude`，不能称为已经验证的 bulk Debye length、equilibrium surface potential 或 zeta
potential。

## 8. Correction 前后的 fit quality

| branch | pass | weak | unusable | median usable R2 | median usable RMSE (pN) |
|---|---:|---:|---:|---:|---:|
| Original | 24 | 7 | 1 | 0.979 | 48.3 |
| Shared empirical 20-200 nm | 9 | 16 | 7 | 0.937 | 17.4 |
| Shared empirical 20-250 nm | 12 | 16 | 4 | 0.950 | 16.4 |
| Block-specific 20-200 nm | 5 | 22 | 5 | 0.919 | 18.7 |

Corrected RMSE 明显降低，说明大幅度经验 velocity-distance component 已从 force shape 中
移除。但 pass fits 从 24 降至 9 同样重要：大信号被移除后，剩余曲线对两个 PB shape
parameters 的约束变弱。视觉上更平的 curve 不自动意味着 equilibrium PB parameters 更
可信或更精确。

## 9. Correction 后残余速度相关性

使用所有未被标为 unusable 的 fits，global Spearman speed association 为：

| apparent parameter | original rho(speed) | corrected rho(speed) |
|---|---:|---:|
| Debye length | +0.968 | -0.304 |
| Surface-potential magnitude | +0.968 | +0.414 |

Primary corrected branch 中，各 nominal speed 的 usable-map medians 为：

| speed (um/s) | median apparent lambda_D (nm) | median surface-potential magnitude (mV) |
|---:|---:|---:|
| 0.1 | 11.78 | 37.99 |
| 0.3 | 9.65 | 43.40 |
| 0.9 | 9.57 | 46.87 |
| 2.7 | 9.51 | 39.72 |

原来的强单调上升排序大幅塌缩，但没有变成完全共同的 speed-independent parameter：
Debye-length association 变为轻度负值，surface-potential magnitude 仍保留非单调的正
association。

Palindrome-controlled visualization 见
[speed comparison figure](figures/pb_parameters_speed_palindrome_original_vs_corrected.png)。

## 10. Correction 后的时间趋势

Primary shared branch 使用以下描述模型：

```text
log(parameter) ~ elapsed_hours + speed categorical effects.
```

HC3 interval 用于表示趋势尺度，不应读作 independent-replicate kinetic model。

| apparent parameter | usable maps | median same-speed rho(time) | adjusted change/hour | descriptive HC3 95% interval |
|---|---:|---:|---:|---:|
| Debye length | 25 | -0.890 | -28.8% | [-41.0%, -14.0%] |
| Surface-potential magnitude | 25 | -0.502 | -8.4% | [-15.1%, -1.2%] |

[corrected-only time figure](figures/pb_parameters_time_empirical_corrected_detail.png)
只在同一 nominal speed 内连接 usable observations；连接线只是 visual guide，不是 fitted
kinetics。[common-scale time comparison](figures/pb_parameters_time_original_vs_empirical_corrected.png)
则显示 corrected parameter range 相对于 original apparent parameters 的量级变化。

最高速 potential subset 只有 4 个 usable fits，且其 within-speed rank correlation 为正，
不能单独用来否定或建立整体时间趋势。

## 11. 0.1 与 0.3 um/s 是否已经进入 equilibrium plateau

当前证据不能建立该 plateau。

先在每个 block 和 speed 内平均两张 symmetric maps，再计算 block 内 `0.3-0.1 um/s`
parameter difference。四个 blocks 的结果为：

| corrected apparent parameter | mean paired difference | paired SE | descriptive t 95% interval | abs(mean)/SE |
|---|---:|---:|---:|---:|
| Debye length | -2.43 nm | 0.55 nm | [-4.16, -0.69] nm | 4.5 |
| Surface-potential magnitude | +4.88 mV | 0.67 mV | [+2.74, +7.01] mV | 7.3 |

四个 Debye-length differences 全部为负，四个 potential differences 全部为正。残余
parameter difference 并不小于 descriptive paired standard error。

在 force-shape level，进一步消去每个 block 的 250 nm speed-only constant：

```text
Delta_shape(D) = [F0.3(D)-F0.1(D)]
                 - [F0.3(250)-F0.1(250)].
```

选定 distances 的结果为：

| D (nm) | mean Delta_shape (pN) | paired SE (pN) | descriptive t 95% interval (pN) |
|---:|---:|---:|---:|
| 20 | +14.05 | 20.11 | [-49.94, +78.03] |
| 50 | -25.10 | 5.55 | [-42.77, -7.44] |
| 100 | -14.89 | 5.09 | [-31.10, +1.32] |
| 150 | -6.08 | 1.01 | [-9.29, -2.87] |
| 200 | -1.84 | 3.70 | [-13.63, +9.95] |

部分 distances 的差异与 block variability 相当，但另一些 distances 保留系统性 residual
speed-shape difference。Distance bins 相互关联，不能把上表各行当作独立 multiple tests。

完整可复算结果见
[low_speed_parameter_pair_comparison.csv](low_speed_parameter_pair_comparison.csv) 和
[low_speed_force_shape_comparison.csv](low_speed_force_shape_comparison.csv)。

Equilibrium plateau 应使用预先定义的 physical equivalence margin，例如允许的 force 或
parameter difference。Standard error 不是 equivalence margin；“没有拒绝 zero difference”
也不等于已经证明 equivalence。

## 12. Hyd 参数作用域与 fit range sensitivity

| corrected branch | rho(speed), lambda_D | adjusted lambda_D change/hour | rho(speed), potential magnitude | adjusted potential change/hour |
|---|---:|---:|---:|---:|
| Shared 20-200 nm | -0.304 | -28.8% | +0.414 | -8.4% |
| Shared 20-250 nm | -0.203 | -39.7% | +0.571 | -8.2% |
| Block-specific 20-200 nm | +0.103 | +74.3% | +0.187 | -16.2% |

Apparent-potential time direction 在三条 branches 中均为负。Apparent Debye-length time
direction 在 block-specific subtraction 下反转。因此 potential decrease 的方向相对更稳定，
而 Debye-length time trend 对是否允许 block history 进入 hyd parameters 非常敏感。

Block-specific result 不是“更完整 correction”的证据。它说明 correction 本身能够吸收
一部分原本希望研究的 history。

## 13. 当前数据支持与不支持的命题

当前数据支持：

- constant baseline 没有消除观察到的 velocity-distance interaction；
- 在分析区间内，`K v^alpha/(D+D0)` 对 interaction 的描述明显优于 strict nominal-speed
  `v/D`；
- 条件减去该经验式后，原来 apparent PB parameters 的大部分速度排序消失，force residual
  RMSE 降低；
- `0.1` 与 `0.3 um/s` 在 apparent PB parameters 和部分 force-shape distances 上仍保留
  residual difference；
- corrected apparent potential 的时间方向比 corrected apparent Debye length 更稳定。

当前数据不支持：

- 真实 hydrodynamic force 已被完全或唯一去除；
- nominal scanner speed 等于 instantaneous gap-closing speed；
- `alpha` 测得 glycerol shear thinning；
- `D0` 测得 slip length；
- corrected `lambda_D` 就是 bulk Debye length；
- corrected surface-potential magnitude 就是 equilibrium potential 或 zeta potential；
- `0.1` 与 `0.3 um/s` 已被证明处于 quasi-static plateau。

## 14. 决定性的下一步实验

更强的物理判别需要：

1. 为每条 curve 保留 raw `D(t)`、cantilever deflection 和 timing；
2. 重建 `U_gap(t)=-dD/dt`，而不是直接代入 nominal scanner speed；
3. 增加低于 `0.1 um/s` 的速度，例如 `0.03` 和 `0.05 um/s`；
4. 在独立 blocks 中随机化顺序，或继续使用严格平衡的 palindrome design；
5. 预先定义 force 与 parameter equivalence margins；
6. 增加 approach/retract sign reversal、viscosity、temperature 和 blank/control checks；
7. 检验一条共同 zero-speed surface-force curve 能否在不使用 block-specific hyd parameters
   的条件下预测所有速度。

只有当最低两个速度的 contrast interval 完全落入预定义 equivalence margin，而且再降低
速度不再改变 inferred surface curve 时，才应宣称进入 quasi-static plateau。

## 15. 复现方法与文件导航

在 repository root 运行：

```bash
python analysis/analyze_12_09_26_997glycerol_hydrodynamic_scaling.py
python analysis/analyze_12_09_26_997glycerol_empirical_hyd_pb.py --workers 4
sha256sum -c analysis/glycerol_99p7_D3_hydrodynamic_scaling_results/artifact_manifest.sha256
sha256sum -c analysis/glycerol_99p7_D3_empirical_hyd_pb_results/artifact_manifest.sha256
```

只有在 corrected-force 与 PB-fit CSVs 已经验证后，才能使用 `--postprocess-only`。该模式
复用 96 个 saved nonlinear corrected fits，只重建 summaries、figures、documentation、
provenance 和 manifest。

主要文件：

- [REPORT.md](REPORT.md)：自动生成的紧凑结果摘要；
- [map_force_empirical_hyd_corrected.csv](map_force_empirical_hyd_corrected.csv)：完整减法与
  gauge-centred corrected curves；
- [map_pb_fits_comparison.csv](map_pb_fits_comparison.csv)：original 与 corrected PB fits、QC
  和 transformed offsets；
- [parameter_time_speed_associations.csv](parameter_time_speed_associations.csv)：global、
  same-speed associations 与 speed-adjusted time summaries；
- [same_speed_parameter_time_trends.csv](same_speed_parameter_time_trends.csv)：逐速度时间趋势；
- [palindrome_parameter_pair_means.csv](palindrome_parameter_pair_means.csv)：symmetric map-pair
  parameter means；
- [empirical_hyd_parameter_sets.csv](empirical_hyd_parameter_sets.csv)：shared 与 block-specific
  empirical parameters；
- [fit_quality_summary.csv](fit_quality_summary.csv)：pass、weak、unusable、R2 与 RMSE；
- [artifact_manifest.sha256](artifact_manifest.sha256)：本目录全部 published artifacts 的
  SHA-256。

实现脚本：

- [analyze_12_09_26_997glycerol_hydrodynamic_scaling.py](../analyze_12_09_26_997glycerol_hydrodynamic_scaling.py)
- [analyze_12_09_26_997glycerol_empirical_hyd_pb.py](../analyze_12_09_26_997glycerol_empirical_hyd_pb.py)

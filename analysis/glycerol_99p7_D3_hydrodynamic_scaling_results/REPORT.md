# 12-09-26 - 99.7 wt% glycerol - hydrodynamic scaling test

## 结论

**数据验证了 baseline 校正后仍存在很强的 velocity-distance coupling，因此纯粹的 `F_hyd=F(v)` 不能描述残余信号；但数据不验证把 nominal scanner speed 直接代入的严格 no-slip `F_hyd=6 pi eta R^2 v/D`。**

在每个 8-map block 内，先平均同速度的 palindrome 对称 map，再把 distance-by-speed force matrix 双重中心化。这样会精确消去任意的 distance-only 项 `A(D)` 和任意的 speed-only 常数 `B(v)`；如果 hydrodynamics 只剩 `B(v)`，中心化后的 interaction 应为零。实际 interaction 的第一 singular component 在四个 block 中解释 99.90%-99.95% 的平方幅值，说明残余不是随机杂乱项，而是高度可复现、近似 separable 的 `v x D` 耦合。
四个 block 的 interaction 图样两两 Pearson correlation 为 0.9972-0.9994；这描述跨 block 的形状复现度，不把相关的 distance bins 当作独立样本做显著性检验。

但该耦合的形状不是严格 `v/D`。20-200 nm 内，自由拟合振幅的 `K v/D` 仍留下 interaction RMS 的 25.6%-28.7%；经验式 `K v^alpha/(D+D0)` 只留下 3.3%-3.8%，四个 block 给出 `alpha=0.588-0.719`、`D0=33.5-40.6 nm`。后者只是 compact empirical descriptor，不是 shear-thinning exponent 或 slip length 的识别。

有限距离 baseline 的物理含义也需要写清楚。若 baseline 对应间距窗口 `D_b`，理想 lubrication force 经过 constant subtraction 后成为

```text
F_corr(D,v) = F_surface(D) - <F_surface(D_b)> + K v [1/D - <1/D_b>]
```

所以 residual 本来就不应消失；baseline 只把 hydrodynamic zero 移到有限距离。四个速度使用不同 `D_b` 时，两个被减去的窗口平均都可产生 speed-only offset。下面的双中心化正是为了把该 offset 消掉，只检验剩余的 distance-dependent interaction。

## 检验定义

设同一 block 内、同速度对称 pair 的平均力为 `F_b(D,v)`。分析量为

```text
I_b(D,v) = F_b(D,v) - <F_b>_v - <F_b>_D + <F_b>_{D,v}
```

因此 `F_surface(D)`、constant baseline subtraction 留下的 speed-only offset，以及 block 常数都不进入 `I_b`。若 `F_hyd=K v/D`，则必须有

```text
I_b(D,v) = K [1/D - <1/D>_D] [v - <v>_v].
```

5 nm bins 来自同一批曲线，彼此相关；它们只定义曲线形状，不被当作独立重复。重复单位是四个 acquisition blocks，表中范围是 block range，不是 confidence interval。

## 20-200 nm block 结果

bulk no-slip 参照使用仓库既有 Cheng correlation：99.7 wt%、25.6 C 时 `eta=807.771 mPa s`，`R=4.546849 um`，因此 `K_bulk=6 pi eta R^2=314782 pN nm/(um/s)`。

| block | first SVD fraction | fitted K for v/D | K/K_bulk | v/D NRMSE | alpha | D0 (nm) | flexible NRMSE |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.99900 | 24950 | 0.079 | 0.284 | 0.588 | 36.3 | 0.038 |
| 2 | 0.99915 | 23111 | 0.073 | 0.287 | 0.632 | 40.6 | 0.033 |
| 3 | 0.99923 | 22608 | 0.072 | 0.269 | 0.719 | 39.7 | 0.038 |
| 4 | 0.99949 | 22518 | 0.072 | 0.256 | 0.675 | 33.5 | 0.035 |

`v/D` 的 block-median NRMSE 为 0.277；允许一个 distance shift、但保持 velocity linear 时降至 0.124；同时允许 sublinear nominal-speed dependence 后降至 0.036。固定 bulk viscosity 的 no-slip 预测 NRMSE 为 12.27，量级也不闭合。

## 速度线性检验

为排除 constant baseline 的影响，再取 20 nm 与 200 nm 的 distance double difference，并除以相邻速度差。若 force 对 nominal speed 线性，三个速度区间反推出的 `K` 应相同；实际随速度升高系统性下降：

| speed interval (um/s) | block-median K (pN nm/(um/s)) | block range | median / K_bulk |
|---:|---:|---:|---:|
| 0.1-0.3 | 41470 | 33594-44301 | 0.132 |
| 0.3-0.9 | 27156 | 20879-32174 | 0.086 |
| 0.9-2.7 | 17786 | 16872-18758 | 0.057 |

这不是由一个 speed-only baseline 常数造成的，因为 distance double difference 会把该常数精确消掉。它说明至少以 nominal scanner speed 为横轴时，`F proportional to v` 不成立。

## 物理解读边界

1. **可以说**：残余信号包含强而可复现的 velocity-distance interaction；纯 `F(v)` 被数据排除，hydrodynamic-like distance dependence 是合理解释。
2. **不能说**：本数据已验证 classical no-slip `v/D` 或测得 bulk viscosity。自由 `v/D` 振幅只有 bulk no-slip 系数的约 7%-8%，而且形状残差有系统性。
3. 公式中的速度应是 instantaneous gap-closing speed `U_gap=-dD/dt`，不是自动等于 nominal scanner speed。高黏度下 cantilever deflection/relaxation 会使两者显著不同；当前本地只保留 derived CSV，raw JPK 不在 checkout，因而本包不能重建 `U_gap(D,t)`。
4. `D0` 也可能混合 contact-zero bias、hydrodynamic compliance、slip、粗糙度和控制回路响应，不能直接命名为 slip length。`alpha<1` 可能主要反映 nominal-to-gap velocity mapping，不足以证明 glycerol shear thinning。
5. pair averaging 只压低 block 内一阶时间漂移；四个 block 的幅值仍随 history 改变。距离 bins 的高相关性与仅四个 speeds 也限制了参数 identifiability。

## 下一步判别实验

决定性检验应从每条 raw curve 同时重建 `D(t)` 与 `U_gap(t)=-dD/dt`，再拟合 `F(D,t)=F_surface(D)+6 pi eta R^2 U_gap/D+C_curve`；同时保留 block/palindrome 结构，并检查 approach/retract 的 hydrodynamic sign reversal。若使用 bulk `eta` 后系数、速度线性和 `1/D` 三项同时闭合，才可称为验证 classical drainage law。

## 数值与可复现性

- 重建的 20/50/100/200 nm pair means 与既有 v0 包最大差为 `0.000e+00 pN`。
- synthetic additive-plus-`v/D` self-check 的 interaction 最大误差为 `3.695e-13 pN`，K 相对恢复误差为 `0.000e+00`。
- 将上限扩到 250/300 nm 后，经验式 alpha 的全部 block 范围为 `0.615-0.750`，D0 为 `35.6-47.5 nm`；定性结论不变。
- 所有 nonlinear fits 均使用多起点 bounded least-squares；CSV 保存 solver status。没有从相关 distance bins 构造 p-value 或 nominal CI。

## 文件

- `palindrome_pair_force_matrix.csv`: 20-300 nm block x speed pair means。
- `double_centered_interaction_20_200nm.csv`: 主检验矩阵。
- `block_interaction_reproducibility.csv`: 四个 block 的 pairwise shape correlation 与 scale。
- `interaction_model_fits.csv`: 主范围与四个 fit-window sensitivities。
- `distancewise_speed_slopes.csv`: 固定距离的 all-four-speed descriptive slopes。
- `velocity_interval_double_differences.csv`: 以 200 nm 为 reference 的 baseline-invariant speed contrasts。
- `figures/hydrodynamic_scaling_diagnostics.*`: 速度区间 collapse 与模型误差。
- `provenance.json`, `artifact_manifest.sha256`: 输入 hash、定义、软件与产物校验。

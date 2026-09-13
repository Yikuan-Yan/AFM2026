# 12-09-26 · 99.7 wt% glycerol · D3 · force/PB time figures

## 物理图像

动态 baseline 去掉的是每条 approach 的速度相关 detector offset；它没有去掉随距离增长的 viscous drainage。重建后的 20–200 nm 绝对力随 approach speed 显著增大，并在同速度下随 acquisition history 改变。因此 PB 图中的两个参数只能读作把有限速度曲线投影到 equilibrium PB 形状后的 apparent parameters，不能读作溶液真实 Debye length 或 equilibrium surface potential。

## 对应图

- `figures/map_median_force_slices.png`：与旧图相同的 20/50/100/200 nm map-median absolute-force 时间序列。
- `figures/debye_length_surface_potential_vs_time.png`：与旧图相同的 apparent Debye length / surface-potential 时间序列。

## 本批数据专用重建

- baseline：直接读取 held-out-map 验证过的四个 speed-conditioned 120 nm 窗口；每条曲线只减一个常数。
- contact：每条曲线只在 terminal 10 nm scanner-travel 物理窗口内检测，不使用固定样本数。
- 本批 contact-derived InvOLS = **49.9844 nm/V**（32 个 map contact 中位数的中位数）；旧 D3 值 91.5630 nm/V 不满足本批硬接触单位斜率。
- force scale = **11.9146 nN/V**，使用独立 D3 `k=0.238366891 N/m`；它与 JPK embedded product 11.8360 nN/V 很接近，但这种乘积一致不等于两项独立标定都正确。
- 每个 map 可接受 terminal contacts 为 37–64/64；重建 endpoint separation 的 map median 绝对值最大为 0.387 nm。

## 绝对力范围

| separation (nm) | 32-map median-force range (pN) |
|---:|---:|
| 20 | 582.4–3572.1 |
| 50 | 217.2–1845.3 |
| 100 | 57.1–858.7 |
| 200 | 4.2–156.2 |

## PB 映射质量与边界

- 24/32 pass，7/32 weak，1/32 unusable。判据只检查 optimizer、boundary、Jacobian、R²、信号和局部参数尺度。
- 模型：equal-potential nonlinear PB Derjaguin + van der Waals + constant offset，20–250 nm；假定 `R=4.546849 µm`、`T=25.6 °C`、`εr=42.5`。
- `εr` 外部参考：https://doi.org/10.3390/ma11040650；它仍是当前样品在假定温度下的模型输入，不是从本批 AFM 数据识别出的参数。
- 未做 hydrodynamic subtraction；99.7% glycerol 的速度相关瞬态和 drainage 是 load-bearing model mismatch。surface potential 还条件依赖 assumed dielectric constant、probe radius、contact zero 和 force calibration。

## 描述性时间/速度相关

| observable | n maps | rho(speed) | rho(time) | median same-speed rho(time) |
|---|---:|---:|---:|---:|
| absolute_force_20nm | 32 | 0.969 | -0.194 | -0.798 |
| absolute_force_50nm | 32 | 0.969 | -0.181 | -0.798 |
| absolute_force_100nm | 32 | 0.969 | -0.169 | -0.976 |
| absolute_force_200nm | 32 | 0.920 | -0.213 | -0.690 |
| apparent_debye_length | 31 | 0.968 | -0.008 | -0.238 |
| apparent_surface_potential | 31 | 0.968 | -0.142 | -0.798 |

这些 rho 是 sequential、non-randomized acquisition 的描述统计；speed 与 time/history 仍不能因一张相关图而因果分离。

## 文件

- `contact_calibration_curves.csv` / `contact_calibration_maps.csv`：物理窗口接触标定和 QC。
- `map_force_by_separation.csv` / `map_reconstruction_diagnostics.csv`：绝对力时间图的数值来源。
- `map_pb_apparent_fits.csv`：PB optimum、残差和 identifiability flags。
- `rank_associations.csv`：map-level 速度/时间秩相关。
- `provenance.json` / `artifact_manifest.sha256`：假设、版本与输出哈希。

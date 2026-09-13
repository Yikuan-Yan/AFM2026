# 12-09-26 · 99.7 wt% glycerol · linear decay and v→0 extrapolation

## 物理图像

同一 approach speed 的 map-median force 随实验进行通常下降；这可用一条描述性直线量化。每个 8-map block 内，同速度的两张 map 位于 palindrome 对称位置；先取 pair mean，能在时间趋势近似线性时消去一阶漂移，再把四个速度的 pair means 线性外推至 v=0。

## 一小时衰减

定义：`100 × [F_fit(t=0) − F_fit(t=60 min)] / F_fit(t=0)`；正值表示衰减。t=0 是第一张 map 的 midpoint。回归单位是 map，每个速度 n=8；区间为 HC3/delta-method 95% interval。

| D (nm) | speed (µm/s) | fitted equation, t_h from first map midpoint (pN) | expected decay in first 1 h | 95% interval | R² |
|---:|---:|---:|---:|---:|---:|
| 20 | 0.1 | `F=720.7-86.3 t_h` | 12.0% | [1.2, 22.8]% | 0.357 |
| 20 | 0.3 | `F=1170.0-164.5 t_h` | 14.1% | [5.4, 22.7]% | 0.596 |
| 20 | 0.9 | `F=2075.8-337.8 t_h` | 16.3% | [9.2, 23.3]% | 0.786 |
| 20 | 2.7 | `F=3554.2-362.1 t_h` | 10.2% | [5.9, 14.5]% | 0.698 |
| 50 | 0.1 | `F=250.7-18.8 t_h` | 7.5% | [-12.5, 27.5]% | 0.189 |
| 50 | 0.3 | `F=464.7-40.8 t_h` | 8.8% | [-0.8, 18.4]% | 0.506 |
| 50 | 0.9 | `F=951.7-106.2 t_h` | 11.2% | [5.7, 16.7]% | 0.838 |
| 50 | 2.7 | `F=1845.8-189.5 t_h` | 10.3% | [7.6, 13.0]% | 0.928 |
| 100 | 0.1 | `F=67.5-2.3 t_h` | 3.4% | [-50.9, 57.8]% | 0.015 |
| 100 | 0.3 | `F=183.3-25.2 t_h` | 13.7% | [6.8, 20.7]% | 0.871 |
| 100 | 0.9 | `F=429.9-61.0 t_h` | 14.2% | [8.7, 19.6]% | 0.899 |
| 100 | 2.7 | `F=852.9-117.8 t_h` | 13.8% | [8.4, 19.2]% | 0.940 |
| 200 | 0.1 | `F=13.4+1.3 t_h` | -9.5% | [-209.0, 190.0]% | 0.007 |
| 200 | 0.3 | `F=51.0-6.1 t_h` | 12.0% | [-11.1, 35.1]% | 0.308 |
| 200 | 0.9 | `F=116.9-15.2 t_h` | 13.0% | [8.3, 17.6]% | 0.884 |
| 200 | 2.7 | `F=156.7-44.0 t_h` | 28.1% | [4.3, 51.8]% | 0.789 |

OLS 与 Theil–Sen slope 在 15/16 个 distance×speed 组合中同号。唯一异号组合是 0.1 µm/s、100 nm；在 200 nm 两种方法都给出很小的正 slope。二者的区间都跨过零，不能解释为已分辨的衰减或增长率。

## block 内 v=0 外推

每个 block 的四个 pair centres 偏离该 block midpoint 最多 0.308 min，因此 pair averaging 对一阶时间项的抵消误差很小。主结果使用全部四个速度；`≤0.9` 三点截距仅作为 2.7 µm/s 瞬态敏感性检查。

| D (nm) | all-4-speed v=0 median [block range] (pN) | ≤0.9-speed median [range] (pN) | median difference (pN) | all-4 R² median [range] |
|---:|---:|---:|---:|---:|
| 20 | 739.9 [641.6, 831.0] | 556.0 [491.1, 645.5] | 177.8 | 0.978 [0.968, 0.988] |
| 50 | 264.4 [248.8, 290.3] | 176.7 [157.6, 195.9] | 87.7 | 0.983 [0.980, 0.984] |
| 100 | 83.7 [81.4, 96.5] | 35.6 [30.3, 40.2] | 49.7 | 0.977 [0.974, 0.977] |
| 200 | 36.6 [31.7, 37.8] | 9.2 [-0.8, 12.1] | 27.5 | 0.705 [0.628, 0.846] |

虽然 20–100 nm 的 all-4 speed 直线通常有高 R²，但把 2.7 µm/s 去掉后 v=0 截距明显下降；这种 fit-range dependence 在 100–200 nm 尤其大。因此这里得到的是 protocol-dependent apparent zero-speed intercept，不是已识别的 equilibrium surface force。四点 fit 的单-block 截距 95% intervals 也很宽，详见 CSV。
本次 16/16 个单-block all-4-speed 截距的 HC3 95% interval 都包含零；点估计可用于比较 block 和 fit range，但不能单独证明正的零速表面力。

## 图与数值文件

- `figures/same_speed_linear_decay_vs_time.png`：同速度时间直线、HC3 mean-response bands 与一小时衰减。
- `figures/blockwise_speed_extrapolation_v0.png`：四个 palindrome blocks 的 pair means 和 v=0 外推。
- `same_speed_time_linear_fits.csv`：16 个时间回归及 HC3/delta-method interval。
- `blockwise_speed_pair_means.csv`：每个对称 pair 的两张 map、均值、半差和时间中心。
- `blockwise_speed_linear_fits.csv`：每个 block×distance 的 all-4 与 low-3 sensitivity fits。
- `v0_summary.csv`：四个 block 的 v=0 汇总。

## 解释边界

- 线性时间模型描述约 1.5 h 的本次序列，不证明长期指数衰减机制。
- pair averaging 只压低一阶时间漂移；surface/probe conditioning、raster history 与速度相关 baseline/transient 仍可能存在。
- v=0 的数学截距不等于 thermodynamic equilibrium，除非 residual velocity dependence、contact/baseline 系统误差及等待时间效应另有控制。

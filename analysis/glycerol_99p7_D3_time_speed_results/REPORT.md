# 12-09-26 · 99.7 wt% glycerol · D3：时间/速度相关性

## 物理图像

动态 baseline 校正后，剩余信号仍可能同时随 approach speed、实验经过时间、surface/probe conditioning 和 raster history 改变。四个 rotated palindrome blocks 让每个速度在每个 block 内出现两次，并位于对称 acquisition position；因此可以分别观察同速度早晚变化，以及先对称平均、再比较速度排序。这个设计能压低一阶时间漂移的混淆，但不能把顺序实验变成随机化速度因果实验。

## 分析定义

- 输入是 speed-conditioned dynamic-baseline 包；approach 每条曲线只做已验证的 constant-median baseline subtraction。
- 统计单位是 8×8 map median，不把 64 pixels 当成独立重复。
- 只计算 Spearman rank correlation、same-speed chronology、palindrome later−earlier exact difference 和 early/late arithmetic pair mean；不做 regression、partial correlation、time model 或 surface-force fit。
- 不报告常规 p-value，因为 32 张 map 是串行、相关、非随机化观测；ρ 只描述当前序列的 rank association。
- scanner travel 不是 tip–surface separation；320–360/380–420 nm 是固定 travel window，而不是绝对表面距离。

## Protocol 的 speed/time 平衡

全部 32 maps 中 speed 与实际 map-midpoint elapsed time 的 Spearman ρ = **0.000**。每个 palindrome block 内使用原始 8-map 顺序计算的 speed/time ρ 为：block 1 0.000, block 2 0.000, block 3 0.000, block 4 0.000。这说明速度标签与线性顺序基本正交，但非线性 drift、block 差异和速度依赖 map duration 仍然存在。

同速度 early/late pair 的实际时间间隔为 1.39–21.64 min；各 pair 的平均时间偏离对应 block midpoint 最多 0.308 min。因此 arithmetic pair mean 很接近共同 block 时间中心，但 later−earlier difference 的幅值不能解释成统一时间间隔下的变化率。

## Map-level rank association

| Observable | global ρ(speed) | global ρ(time) | same-speed ρ(time), median [range] | block pair-mean ρ(speed), median [range] | pair-mean 2.7−0.1, median [range] |
|---|---:|---:|---:|---:|---:|
| Dynamic baseline − initial endpoint (V) | 0.969 | -0.133 | -0.607 [-0.738, -0.167] | 1.000 [1.000, 1.000] | 0.706076 [0.688887, 0.801683] |
| Absolute dynamic-baseline detector level (V) | 0.421 | 0.175 | 0.286 [0.119, 0.429] | 1.000 [0.800, 1.000] | 0.767987 [0.750352, 0.958882] |
| Corrected increment at 320–360 nm travel (V) | 0.969 | -0.173 | -0.940 [-0.976, 0.119] | 1.000 [1.000, 1.000] | 0.0451219 [0.0409588, 0.0475716] |
| Corrected increment at 380–420 nm travel (V) | 0.969 | -0.215 | -0.893 [-0.929, -0.738] | 1.000 [1.000, 1.000] | 0.0757933 [0.0739591, 0.07721] |
| Corrected terminal last 25 nm (V) | -0.969 | 0.231 | 0.905 [0.905, 1.000] | -1.000 [-1.000, -1.000] | -0.705053 [-0.752424, -0.685915] |
| Within-baseline-window drift magnitude (V) | 0.969 | -0.088 | -0.321 [-0.833, 0.071] | 1.000 [1.000, 1.000] | 0.0255019 [0.0246325, 0.0284948] |
| Retract minimum, endpoint referenced (V) | 0.354 | 0.531 | 0.464 [0.333, 0.881] | 0.800 [0.400, 1.000] | 0.149136 [0.0342186, 0.291475] |
| Travel position of retract minimum (nm) | 0.100 | 0.533 | 0.560 [0.357, 0.905] | 0.100 [-0.200, 0.800] | -0.34296 [-2.58169, 11.6148] |

`global ρ(speed)` 利用全部 32 maps；`same-speed ρ(time)` 是四个速度各自八张 map 的结果再取中位数；`block pair-mean ρ(speed)` 在每个 block 内只剩四个速度点，因此 ρ 很离散，只用于方向/排序诊断。

## 主要观测结果

1. **中距离增量保留一致的正速度排序。** 320–360 nm 的 global ρ(speed) = 0.969，四个 block 的 pair-mean ρ(speed) 全为 1；2.7−0.1 µm/s pair-mean 差的 block median 为 0.045122 V。380–420 nm 对应值为 ρ = 0.969、四个 block 均为 1、差值 0.075793 V。
2. **同一速度仍有明显 time/history relaxation。** 320–360 nm 的 same-speed ρ(time) 中位数为 -0.940，范围 -0.976–0.119；其中 0.1 µm/s 是弱瞬态/噪声主导的例外。380–420 nm 的四个速度均为负，范围 -0.929–-0.738。
3. **速度排序只在一段 travel 区域内稳定。** global 和 block-pair-median speed ρ 同时不低于 0.75 的连续主区间是约 237.5–417.5 nm；在早期区间看到的是相对动态 baseline 的 startup tail，约 427.5 nm 后则进入 terminal/contact-like 排序翻转。
4. **terminal response 与中距离相反。** terminal global ρ(speed) = -0.969，四个 block 的 pair-mean ρ 均为 −1；same-speed ρ(time) 中位数却为 0.905。这表明 terminal load/contact alignment 与中距离增量不是同一个可直接互换的 observable。
5. **retract 更突出 history。** retract minimum 的 global ρ(time) = 0.531，same-speed ρ(time) 中位数 = 0.464；later−earlier 多数为正，即负向 minimum 随实验进行通常变得较浅。

## Same-speed palindrome history

下表给出每个 block 中同速度 later−earlier 的四个 exact differences 的中位数。若它明显不为零，说明即使速度固定，观测量仍随该 block 内的时间/history 改变。

| Observable | 0.1 µm/s | 0.3 µm/s | 0.9 µm/s | 2.7 µm/s |
|---|---:|---:|---:|---:|
| Dynamic baseline − initial endpoint (V) | 0.000480039 | -0.00514092 | -0.00530835 | -0.00985456 |
| Corrected increment at 320–360 nm travel (V) | -0.000869801 | -0.000200315 | -0.00137633 | -0.0012206 |
| Corrected increment at 380–420 nm travel (V) | -0.0030507 | -0.000496599 | -0.00242424 | -0.00301255 |
| Corrected terminal last 25 nm (V) | 0.00796072 | 0.00893864 | 0.0100861 | 0.0144524 |
| Retract minimum, endpoint referenced (V) | 0.0993448 | 0.075802 | 0.0498176 | 0.0376992 |

## 解释边界

- 正的 speed rank association 只说明在本 protocol 的 pair-mean map signal 中较快 approach 通常对应较大增量；它不等于已识别的 hydrodynamic coefficient。
- same-speed time association 与 later−earlier difference 是 history/systematic 的直接证据；它们限制任何把跨速度差异全部归因于速度的解释。
- dynamic baseline 本身包含速度相关 motion-start response；绝对 baseline detector level 还会漂移和 reset，因此二者分别报告。
- retract speed 始终为 1 µm/s。retract minimum 对 preceding approach speed 的相关性是接触历史关联，不是 retract velocity dependence。
- terminal last-25-nm response 属于 contact-like 区域且随速度排序可与中距离相反，不能作为相同 separation 的 equilibrium force 比较。

## 文件

- `map_time_speed_observables.csv`：32-map chronology 与八个 map-level observables。
- `rank_associations.csv`、`rank_association_summary.csv`：global、same-speed 与 block pair-mean rank associations。
- `palindrome_history_pairs.csv`、`palindrome_history_summary.csv`：同速度 exact later−earlier differences 和 pair means。
- `map_binned_corrected_curves.csv`：每张 map 的 5 nm corrected median curve。
- `distancewise_rank_associations.csv`、`distancewise_palindrome_pairs.csv`：逐 scanner-travel bin 的时间/速度诊断。
- `figures/`：chronology、distance-resolved correlation、scalar summary 与 palindrome speed profiles。
- `provenance.json`、`artifact_manifest.sha256`：输入身份、算法边界和产物哈希。

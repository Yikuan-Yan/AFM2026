# 27-08-26 D4–D6 标定与 D4 纯水回文测试：两组 pilot

> 当前只有计划中的前两组回文 block。本报告完成基础标定、raw map 重建、同点位早/晚配对与描述性速度拟合；不把两组数据称为完整的 zero-speed extrapolation validation。

## 直接结果

用户给定实际温度为 **25.6 °C**，标定环境为空气；TND header 中记录的是 25.0 °C，本计算明确以 25.6 °C 覆盖 header。热谱全部从 measured PSD 重搜 fundamental 并重拟合，未使用导出的 fit-data。

| Cantilever | Air InvOLS (nm/V) | k (N/m) | f0 (kHz) | Q |
|---|---:|---:|---:|---:|
| D4 | 72.76 ± 3.76 | 0.1648 ± 0.0173 | 13.1394 ± 0.0023 | 49.61 ± 0.95 |
| D5 | 79.90 ± 2.69 | 0.3113 ± 0.0217 | 16.5647 ± 0.0014 | 69.81 ± 1.33 |
| D6 | 70.90 ± 1.06 | 0.1817 ± 0.0066 | 13.2811 ± 0.0018 | 52.48 ± 0.82 |

± 为 5 次重复的 sample SD，不是 traceable absolute uncertainty。由于 D5 的硬接触段只有约 40 nm，三支 cantilever 统一使用末端 35 nm 作为 primary contact fit，并以 30/40 nm 作为 window sensitivity；直接固定 50 nm 会把 D5 的接触转折混入拟合。

纯水 map 从全部有效硬接触重新得到 global InvOLS = **54.214 nm/V**；力使用空气热标定得到的 D4 **k = 0.164829 N/m**。JPK 文件中写入的 sensitivity 和 force conversion 均未作为最终力标尺。

## 回文采集与 map-level QC

- Block 1 时间顺序：`0.2 → 0.1 → 0.05 → 0.05 → 0.1 → 0.2 µm/s`。
- Block 2 时间顺序：`0.1 → 0.05 → 0.2 → 0.2 → 0.05 → 0.1 µm/s`。
- 12 张 map 均为 8×8 pixels、10×10 µm，同一物理区域；每张 map 是一个速度实验单位，pixels 只作为同点位 paired spatial observations。

| Block | Speed (µm/s) | Map InvOLS early → late (nm/V) | Gap speed early → late (µm/s) | Far slope early → late (pN/100 nm) |
|---:|---:|---:|---:|---:|
| 1 | 0.05 | 53.794 → 53.886 | 0.0491 → 0.0495 | -1.16 → -2.88 |
| 1 | 0.1 | 53.703 → 54.007 | 0.0972 → 0.0994 | -1.17 → -5.88 |
| 1 | 0.2 | 53.591 → 54.009 | 0.1931 → 0.1991 | 0.51 → -6.54 |
| 2 | 0.05 | 54.367 → 54.618 | 0.0474 → 0.0478 | -6.76 → -7.76 |
| 2 | 0.1 | 54.443 → 54.687 | 0.0951 → 0.0964 | -6.25 → -11.04 |
| 2 | 0.2 | 54.579 → 54.613 | 0.1899 → 0.1903 | -6.65 → -7.65 |

Retract segment 的 protocol 在所有 map 中相同，因此下列变化不能简单归因于 approach nominal speed；它们主要作为随时间/接触历史的 systematics probe。

| Order | Block | Approach U | Approach snap fraction | Retract snap-off fraction | Detachment travel (nm) | Pull-off (nN) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 0.2 | 0.69 | 0.98 | 99.7 | -1.789 |
| 2 | 1 | 0.1 | 0.73 | 0.98 | 96.7 | -1.752 |
| 3 | 1 | 0.05 | 0.52 | 0.98 | 105.2 | -1.776 |
| 4 | 1 | 0.05 | 0.53 | 0.97 | 105.9 | -1.682 |
| 5 | 1 | 0.1 | 0.31 | 1.00 | 84.7 | -1.495 |
| 6 | 1 | 0.2 | 0.14 | 0.97 | 80.6 | -1.150 |
| 7 | 2 | 0.1 | 0.42 | 0.89 | 93.1 | -0.902 |
| 8 | 2 | 0.05 | 0.33 | 1.00 | 101.4 | -1.055 |
| 9 | 2 | 0.2 | 0.11 | 0.92 | 102.5 | -0.557 |
| 10 | 2 | 0.2 | 0.08 | 0.94 | 100.6 | -0.403 |
| 11 | 2 | 0.05 | 0.20 | 0.95 | 98.8 | -0.736 |
| 12 | 2 | 0.1 | 0.17 | 0.95 | 88.4 | -0.582 |

Map-order trend 是描述性统计（连续顺序采集、并非 randomized independent maps），但可用于识别主要混杂量：

| Metric | first → last | slope/map | Spearman ρ | two-sided p |
|---|---:|---:|---:|---:|
| water_InvOLS_nm_per_V | 53.6 → 54.7 | 0.108 | 0.993 | 1.3e-10 |
| far_slope_pN_per_100nm | 0.506 → -11 | -0.884 | -0.979 | 3.09e-08 |
| force_50nm_pN | 427 → 441 | 31.6 | 0.517 | 0.0849 |
| force_100nm_pN | 34.2 → 27.1 | 3.84 | 0.434 | 0.159 |
| approach_snap_detected_fraction | 0.688 → 0.172 | -0.0535 | -0.797 | 0.0019 |
| retract_pull_off_force_nN | -1.79 → -0.582 | 0.137 | 0.930 | 1.17e-05 |

Block 1 的 map-median InvOLS 为 53.832 ± 0.168 nm/V，Block 2 为 54.551 ± 0.121 nm/V，均值差 +0.720 nm/V（+1.34%）。这是时间相关 optical/contact response 变化的直接证据；primary force 仍使用预先声明的全数据 global InvOLS，未逐图重新缩放。

## 同点位回文分解与速度 pilot

对每个 block、每个速度，将相同 pixel 的早/晚曲线定义为：

`F_sym = (F_early + F_late)/2`，`F_history = (F_late - F_early)/2`。

下表均为逐 curve far-field 线性 drift 修正后的 map-pixel median。每个 block 内仅有三个 F_sym map-level 速度点，因此 slope 与 U→0 intercept 是描述性 pilot。

| Block | D (nm) | Fsym 0.05/0.1/0.2 (pN) | Fhistory 0.05/0.1/0.2 (pN) | observed dF/dU | hyd theory dF/dU | ratio | R² |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 20 | 1870.2 / 1892.6 / 1920.8 | -114.7 / -405.3 / -552.5 | 336.9 | 17.16 | 19.6 | 0.985 |
| 1 | 50 | 197.9 / 230.1 / 254.8 | -40.1 / -120.0 / -175.3 | 368.7 | 6.86 | 53.7 | 0.931 |
| 1 | 100 | 11.9 / 14.4 / 17.5 | -4.4 / -10.2 / -16.2 | 37.8 | 3.43 | 11.0 | 0.983 |
| 1 | 200 | 2.5 / 4.2 / 3.6 | -0.5 / -0.0 / -0.3 | 5.6 | 1.72 | 3.3 | 0.253 |
| 2 | 20 | 2652.0 / 2481.2 / 2691.3 | -60.1 / -102.5 / -20.3 | 545.5 | 17.16 | 31.8 | 0.125 |
| 2 | 50 | 601.1 / 557.4 / 623.6 | -74.4 / -112.8 / -20.3 | 233.3 | 6.86 | 34.0 | 0.252 |
| 2 | 100 | 63.8 / 62.1 / 67.1 | -22.1 / -34.3 / -5.3 | 27.1 | 3.43 | 7.9 | 0.593 |
| 2 | 200 | 3.4 / 5.0 / 5.6 | -0.6 / -2.7 / -1.1 | 13.9 | 1.72 | 8.1 | 0.793 |

水在 25.6 °C 的 Cheng correlation viscosity 为 0.8806 mPa·s；球半径沿用 R=4.546849 µm。no-slip lubrication prediction 为 `Fhyd=6πηR²U/D`：在 100 nm、0.05/0.1/0.2 µm/s 仅为 0.172/0.343/0.686 pN。当前两组在 100 nm 的 |F_history| 中位量级为 13.19 pN，因而物理 hyd signal 远小于时间/history residual；若 observed slope 很大、变号或不呈 1/D，不能解释成 hydrodynamics。

### Far-field 零点定义的影响

- D=20 nm：linear drift corrected − far constant referenced 的 map-median 为 -15.31 pN，范围 -36.90 到 +5.73 pN。
- D=50 nm：linear drift corrected − far constant referenced 的 map-median 为 -16.66 pN，范围 -32.70 到 +0.08 pN。
- D=100 nm：linear drift corrected − far constant referenced 的 map-median 为 -12.98 pN，范围 -21.16 到 +1.53 pN。
- D=200 nm：linear drift corrected − far constant referenced 的 map-median 为 -6.03 pN，范围 -11.22 到 +1.88 pN。

在 100–200 nm，这个 baseline-definition shift 已远大于理论 0.05–0.2 µm/s hyd force。因此第三 block 之后仍应同时报告两种零点定义，并用跨 block 可重复性决定能否分离真正的 1/D hyd term；不能把 far-field straight-line subtraction 当作无物理代价的预处理。

## 基础判断

1. 这套 8×8 回文设计能清楚测出 map-level 的时间/history 漂移，并避免把 64 pixels 当作 64 个速度重复。
2. 水中 InvOLS 与空气 InvOLS 的差异很大，验证了必须从水中 hard contact 重算 sensitivity；同时两 block 间仍有约百分之一量级 response 漂移，需要第三 block 判断它是单调时间漂移、block 跳变还是可重复的顺序效应。
3. Approach snap fraction、retract detachment travel 与 pull-off 随 map order 的变化是独立于 far-field force 的接触历史指标。若它们在同一时段共同跳变，说明 velocity label 之外还存在 surface/contact-state 漂移，不能用单一 hyd term 修正。
4. 0.05–0.2 µm/s 下理论 hyd force 在 20–200 nm 只有亚 pN 到数 pN。第三 block 的主要价值不是再提高 pixel 数，而是检验回文对称化后的 F_sym(U) slope 是否跨 block 同号、同量级、并近似 1/D。
5. 当前力是 apparent finite-speed force：包含 equilibrium surface force、极小的 hyd contribution 和未完全消除的 history/systematic residual。不能据这两组 alone 报告 zeta potential、Debye length 或 validated U→0 force。

## 数值与数据 QC

- 共 12 张 map、768 对 approach/retract curve；所有 map 的两条 branch 均无 parser-skipped curve。
- 输入 JPK ZIP 已逐文件 CRC 检查；TND frequency 严格递增、PSD 为正且无 NaN/Inf；所有 SHO optimizer 报告 success。
- 接触、far field、scanner speed 与 actual gap speed 都在 SI 单位中计算，输出时才转换为 nm、µm/s、pN、nN。
- 同时保留 `far_constant_referenced` 与 `linear_drift_corrected` map force；primary pair analysis 使用后者，far slope 的 8×8 空间图单独保留，避免把 slope 悄悄当作纯仪器项。

## 输出文件

- `calibration_summary.csv`、`thermal_refits.csv`、`calibration_force_contact_fits.csv`：D4–D6 原始标定与窗口敏感性。
- `water_contact_sensitivity_curves.csv`、`map_inventory_QC.csv`、`pixel_QC.csv`：水中 InvOLS、map acquisition/QC、逐 pixel approach/retract 指标。
- `map_force_curves.csv`：两种 baseline 定义下每张 map 的 median/IQR F–D。
- `palindrome_pair_curves.csv`：同点位 F_sym 与 F_history。
- `pilot_velocity_fits.csv`、`force_slices_20_50_100_200nm.csv`：两 block 的描述性三速度拟合。
- `baseline_sensitivity_slices.csv`、`chronological_trends.csv`：零点定义变化和 map-order 描述性趋势。
- `figures/`：标定、时间 QC、8×8 far-slope、回文分解和 hyd slope 对照图。
- `provenance.json`、`artifact_manifest.sha256`：参数、输入 hashes、软件版本与输出身份。

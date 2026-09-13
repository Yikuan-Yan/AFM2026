# 12-09-26 · 99.7 wt% glycerol · D3：速度条件化动态 baseline

## 物理图像

approach 开始后，cantilever/液体/控制回路先经历速度相关瞬态，随后进入一段近似稳定的恒速运动平台，最后才出现靠近表面的快速变化。因此这里的 baseline 不是开头的静止端点，也不是绝对零力；它被定义为**启动瞬态结束后、近表面上升前的 in-motion 动态平台**。减去它以后得到的是相对该运动平台的增量信号。

## 选中的 baseline 模式

- 所有窗口都用 scanner travel 的 nm 定义，绝不使用固定点数；这避免 0.1–2.7 µm/s 不同采样长度造成不同物理范围。
- 每个候选窗口宽 120 nm，并带前后各 20 nm guard。检测分数由窗口首末四分位差、两侧 guard 与窗口中位数的不连续量、5 nm bin 相邻粗糙度组成。
- blocks 1+3 只负责选择，blocks 2+4 完全留出用于评价；分数先在每张 8×8 map 内取中位数，再跨 map 取中位数，避免把 64 pixels 当成独立重复。为避免追逐 5 nm 网格上的单点噪声，在包含全局最小值、且不高于最小值 1% 的连续平底区间中选择中点。
- 窗口一旦按速度选定，每条曲线只减去该窗口 raw vDeflection 的 median；不做线性斜率扣除，避免把可能的距离依赖物理信号拟合掉。

| approach speed (µm/s) | selected window (nm) | time after motion start (s) | train score (mV) | held-out score (mV) | held-out vs fixed 100–220 nm | LOBO start range (nm) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 35–155 | 0.3500–1.5500 | 2.637 | 2.725 | +4.6% | 25–65 |
| 0.3 | 40–160 | 0.1333–0.5333 | 5.186 | 5.044 | +16.5% | 30–40 |
| 0.9 | 70–190 | 0.0778–0.2111 | 13.694 | 12.549 | +15.1% | 65–70 |
| 2.7 | 155–275 | 0.0574–0.1019 | 51.760 | 48.934 | +25.1% | 135–155 |

这里的 detector score 只用于比较窗口稳定性，单位虽为 V，但不是误差条、likelihood 或置信区间。LOBO 是每次留出一个完整 acquisition block 后重新选择窗口；它检查窗口选择是否依赖某一 block。

## 应用后的观测量

下表以八张 map 的 map median 为统计单位。nN 只采用仓库中既有 D3 比例尺作辅助读数，raw V 仍是主结果。

| speed (µm/s) | dynamic baseline − initial endpoint (V) | 320–360 nm corrected (V) | 380–420 nm corrected (V) | terminal 25 nm corrected (V) | within-window drift (V) |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 0.001664 | 0.010000 | 0.051833 | 1.348937 | 0.001160 |
| 0.3 | 0.018108 | 0.019894 | 0.077938 | 1.300846 | 0.002256 |
| 0.9 | 0.168028 | 0.038834 | 0.110846 | 1.137568 | 0.006576 |
| 2.7 | 0.703521 | 0.054851 | 0.128071 | 0.650320 | 0.026846 |

baseline window 自身的 corrected median 按定义为零，因此没有把它当作验证证据。真正的检查是：留出 blocks 的 detector score、窗口两侧 guard、完整 block 的 LOBO 选择，以及校正后曲线是否仍保留 300 nm 之后的变化。

## 结论与边界

1. **必须按速度延迟 baseline 检测。** 最高速 2.7 µm/s 的稳定窗口显著晚于低速；统一使用开头或固定 100–220 nm 会把启动瞬态混进 baseline。
2. **0.1 µm/s 的精确起点约束较弱。** 该速度的瞬态幅度已接近 detector score 的噪声尺度，LOBO 给出 25–65 nm；35–155 nm 是可复现的操作选择，不应解释为精确的 relaxation length。
3. **采用 constant median 而不是 linear detrend。** 这只移除速度/曲线特异的动态平台常数，不主动消去平台内残余斜率；表中的 within-window drift 就是仍保留的非平坦程度。
4. **校正结果是 incremental response，不是 equilibrium force。** 速度相关平台可能含 bulk drag、cantilever relaxation、photodiode/control transient 及真正的远场物理响应。减去它是一种操作性比较基准，不证明这些成分是纯仪器背景。
5. **retract 未使用此 baseline。** 所有 retract 均为 1 µm/s，且已有强负向、历史相关响应；把 approach 的平台规则套给 retract 会改变目标观测量。
6. **没有做 surface-force 拟合。** 本包不做 contact/separation、hydrodynamic、PB、速度外推、回归或机制识别。

## 文件

- `candidate_window_scores.csv`：每个速度、每个 120 nm 候选窗口的 train/held-out map-median 分数及组成。
- `window_selection_stability.csv`：完整 acquisition block 的 leave-one-block-out 选择。
- `curve_dynamic_baselines.csv`：2048 条 approach 的选定窗口、baseline、guard、近表面增量和 detector diagnostics。
- `corrected_curve_bins.npz`：逐曲线 5 nm binned raw-endpoint 与 dynamic-baseline corrected arrays；行号与 curve CSV 对齐。
- `map_dynamic_baseline_summary.csv`、`speed_dynamic_baseline_summary.csv`：map 与 speed 层级汇总。
- `binned_corrected_speed_summary.csv`、`line_position_residuals.csv`：曲线形状与 line-position diagnostics。
- `figures/`：窗口选择、校正效果与留出验证。
- `provenance.json`、`artifact_manifest.sha256`：参数、输入哈希和产物校验。

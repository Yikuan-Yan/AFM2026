# 12-09-26 · 99.7 wt% glycerol · cantilever D3 · observational analysis

## 物理图像与结论

这批数据的主导特征不是一条可直接放入 PB 或 lubrication 拟合的平衡 force–distance 曲线，而是 **高黏度液体中、每次 Z motion 启动后出现的速度相关瞬态和准平台，加上逐行 baseline adjustment 与采集历史留下的行内位置结构**。在 2.7 µm/s 时，approach 从最远端开始后的几十纳米内就上升到很大的正平台；这一变化远早于末端硬接触样的陡升，因此前段不能作为平坦、无力的 far field。固定像素也保留四速度分离，说明它不是单纯由空间平均造成；但全图中行首/行末与内部点明显不同，说明它也不能被解释成样品表面的静态空间力图。

本分析因此只报告 raw `vDeflection`、从最远端起算的 scanner travel、固定物理窗口中位数/IQR、极值以及四个 8-map palindrome block 的同速 later-minus-earlier 差。**没有进行 baseline、contact、hydrodynamic、PB、回归或其他模型拟合，也没有输出 surface separation、Debye length、surface potential 或 `U→0` 外推。**

## 直接观测

### 1. Approach 的启动瞬态与速度平台

每条曲线以最初 0–10 nm 的 raw detector 中位数作一个纯操作性常数参考；100–300 nm 窗口用于描述 motion 启动后已经较稳定的电压增量。它不是绝对零力。下表先给 raw V，再把同一电压差按仓库内既有 D3 独立标定换成仅供量级阅读的 nN：

| approach speed (µm/s) | maps | settled ΔV map median (V) | map range (V) | prior-D3-scale ΔF (nN) | map range (nN) | ΔV / speed (V per µm/s) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 8 | 0.00285 | 0.00181–0.00478 | 0.062 | 0.040–0.104 | 0.02849 |
| 0.3 | 8 | 0.02187 | 0.01464–0.03598 | 0.477 | 0.320–0.785 | 0.07290 |
| 0.9 | 8 | 0.17446 | 0.14686–0.22418 | 3.808 | 3.205–4.893 | 0.19385 |
| 2.7 | 8 | 0.69914 | 0.67292–0.80215 | 15.259 | 14.687–17.507 | 0.25894 |

四个幅值不仅相差很大，`ΔV/speed` 也不恒定；因此在未解决 startup response、baseline semantics 和绝对间距前，不能把这些平台直接当作简单线性 viscous coefficient。图中 2.7 µm/s 平台在约前 100 nm 内形成，而末端约 400 nm 以后才出现所有速度共有的硬接触样上升。最小可信结论是：**运动启动/体相拖曳/仪器响应在 nominal far side 已经是主量级之一**；仅凭本数据不能把三者再分开。

### 2. 逐行采集位置结构

Header 明确记录 baseline adjustment `enabled=true`, `beginOfLine=true`, `interval=1`, `deadtimeBeforeSamples=100`, `averageSamples=100`。按 serpentine raster 的采集位置统计，100–300 nm 的 endpoint-referenced 平台在行首最低、第二点最高，随后向行末下降。以下为各速度 position 1 / 2 / 8 的中位数（nN，仅 prior D3 scale）：

| speed (µm/s) | line start pos 1 | pos 2 | line end pos 8 |
|---:|---:|---:|---:|
| 0.1 | 0.019 | 0.104 | 0.043 |
| 0.3 | -0.041 | 0.828 | 0.128 |
| 0.9 | 0.108 | 4.540 | 1.478 |
| 2.7 | 9.144 | 16.197 | 11.066 |

把每张 map 的 raw V 先减去自身空间中位数后可见：结构主要集中在最初 0–10 nm detector level；100–300 nm 的 settled raw level 平坦得多。因此彩色 heatmap 的边缘/行内图样主要是 motion-start reference 与逐行采集状态，而不是证据充分的表面相互作用空间图。Pixel 仍是相关的 raster 采样，不能作为 64 个独立实验。

### 3. 旧分析方法在这里失效的位置

- 旧流程把 approach 的前 20% 当作 far-field baseline；本实验的前 20% 正好是最强的速度相关弯曲启动段。线性 subtraction 会删掉或重塑待观察的动态信号。
- 四速度保持约 2 kHz 采样率，而不是保持相同点间距：每条 approach 分别约 10000、3330、1110、370 点。旧 contact search 的固定 150 点分别覆盖约 **7.5、22.5、67.5、202.7 nm** scanner travel，物理窗口不等价。
- 高速平台按 prior D3 scale 可达十余 nN，对应几十纳米 cantilever deflection；若把它先作为基线删去或保留在 `h+δ` 中，都会显著改变 inferred separation。没有独立的 absolute force-zero/contact plane 时，本次不构造 surface separation。
- 文件自身的 embedded calibration 与既有 D3 独立标定冲突，进一步限制了绝对力幅值；本报告的科学主结果因此是 V，nN 只是指定标度下的换算。

### 4. Retract 是固定 1 µm/s，但存在大而历史相关的负向最低点

所有 retract 都是 1 µm/s。下表按其前一条 approach speed 分组；颜色分组不是 retract speed 扫描，也不能据此归因于 approach speed：

| preceding approach speed (µm/s) | retract minimum map median (V) | map range (V) | prior-D3-scale median (nN) | map range (nN) | median minimum position (nm from far endpoint) |
|---:|---:|---:|---:|---:|---:|
| 0.1 | -0.87245 | -1.21381–-0.63692 | -19.04 | -26.49–-13.90 | 346.7 |
| 0.3 | -0.93727 | -1.16948–-0.59894 | -20.46 | -25.52–-13.07 | 341.9 |
| 0.9 | -0.79100 | -0.96033–-0.62214 | -17.26 | -20.96–-13.58 | 349.5 |
| 2.7 | -0.81287 | -0.91375–-0.50664 | -17.74 | -19.94–-11.06 | 344.8 |

Retract minimum 的 map spatial IQR 很大，并在时间上先增强、后减弱；其位置也移动。它可以描述为高黏度条件下的 suction/adhesion-like hysteresis，但仅凭曲线不能区分 viscous drainage、cavitation/bridging、surface adhesion、contact history 或 instrument baseline。全部 2048 条 approach 和 2048 条 retract 都未触及 raw encoder saturation，因此这些极值不是 ADC 截断造成的。

### 5. Palindrome 配对只支持“存在 history”，不支持一次速度拟合

四个 block 各自按位置 1↔8、2↔7、3↔6、4↔5 形成 16 个同速 early/late map pair。Approach settled 平台的 later-minus-earlier 中位数为 **-0.063 nN**，范围 **-1.349 到 +0.538 nN**；retract minimum 的对应中位数为 **+1.553 nN**，范围 **-1.501 到 +4.392 nN**。符号并不一致，说明时间/接触历史对两个 branch 的影响不可忽略。这里没有把 16 对做成显著性检验，也没有拟合速度斜率。

Absolute raw detector level 从早期逐步漂移；最大相邻跳变发生在 acquisition order 21→22，initial level 改变 **+2.653 V**。这类 reset/discontinuity 不出现在 encoder conversion 常数中，故跨越该点的 absolute-V 比较不能直接解释成力变化。Endpoint-referenced 增量仍保留，但也不能恢复绝对零力。

## 数据完整性与标度边界

- Keeper share：32 个文件，32 张完整 8×8 map；共 2048 条 approach 和 2048 条 retract。下载清单逐文件记录 SHA-256，JPK ZIP CRC 全部通过；parser skip、instrument failure flag 与 raw saturation 均为 0。
- 每张图覆盖同一 2×2 µm 区域，500 nm nominal Z span；approach 为 0.1/0.3/0.9/2.7 µm/s，各 8 张；retract 固定 1 µm/s。四个 8-map block 都是同速对称 palindrome，且每个 block 的起始速度轮换。
- 既有 D3 独立标定：InvOLS = 91.562952 nm/V，k = 0.238366891 N/m，乘积 = 21.825576 nN/V。
- 这些 JPK 文件内嵌：InvOLS = 62.670658 nm/V，k = 0.188859650 N/m，乘积 = 11.835958 nN/V。Embedded/prior 幅值比 = 0.542。本分析不拟合新 InvOLS，也不裁决哪一套是当前绝对力标尺。

## 可复现产物

- `map_inventory.csv`：文件、时间、四速度设计、采样点、baseline-adjust settings、校准字段与完整性。
- `pixel_observables.csv`：2048 个物理像素的 approach/retract 原始窗口描述量。
- `map_summary.csv`、`speed_summary.csv`、`palindrome_pair_differences.csv`：等 map 权重的描述性汇总。
- `binned_branch_summary.csv`：5 nm scanner-travel bins 的全曲线中位数/IQR；binning 只用于显示，不提高物理分辨率。
- `fixed_pixel_binned_curves.csv`：row 3, column 3 的逐 map 原始曲线，控制空间与行内采集位置。
- `line_position_profile.csv`：行内 position 1–8 的 startup/settled 结构。
- `figures/approach_retract_speed_overview.png`：四速度 approach 与固定速度 retract。
- `figures/startup_and_line_position.png`：startup 时间/位移图与行内位置效应。
- `figures/chronology_and_history.png`：palindrome block、raw-level reset 与 retract history。
- `figures/fixed_pixel_raw_chronology.png`：固定像素的 32-map 对照。

## 当前可说与不可说

可以说：99.7 wt% glycerol 中，motion-start/体相拖曳/仪器响应在 nominal far side 已经达到与表面近场同等级的重要性；逐行采集状态和长时间 history 都清晰存在；retract 有强负向 hysteresis。

不能说：平台已经给出 hydrodynamic coefficient；四组就是平衡 surface force；负向 retract 极值已经识别为某一种 adhesion/cavitation 机制；既有 D3 或 embedded 标定已被本数据重新验证。下一步若要做定量模型，首先需要独立 blank/large-gap motion transient、current D3 calibration，以及不依赖前 20% 的 absolute force-zero/contact strategy。

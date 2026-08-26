# Silica–silica 球–面 AFM 后续实验方案

- **版本：** 1.1
- **日期：** 2026-08-24
- **目标温度：** 25.6 °C
- **体系：** silica colloidal probe–silica plane；水–甘油混合物；不外加盐
- **现阶段主浓度：** 0、20、30、40 wt% glycerol；10 wt% 不进入主实验和主拟合
- **后续扩展范围：** 0–99.5 wt% glycerol；目标 equilibrium force 距离范围 20–200 nm
- **主 cantilever：** cantilever 1，固定使用独立标定的 `k = 0.2969899087 N/m`
**主球半径：** `R = 4.546849 µm`（局部球冠拟合值，不是全颗粒无条件真值）

## 1. 先给出实验决策

推荐采用下面的正式设计，而不是只做一个六-map 回文序列：

- 每张 map 为 `8 × 8 = 64` 条 force curves，保持当前 `10 × 10 µm²` 物理视野不变。
- 每个回文 block 在同一物理区域连续测六张 map，速度在 map 之间改变、map 内固定。
- 每个浓度、每次独立液体制备测三个旋转回文 block：

| 回文类型 | 六张 map 的速度顺序（µm/s） | 作用 |
|---|---|---|
| A | `1 → 2 → 4 → 4 → 2 → 1` | 1 µm/s 位于最外层，4 µm/s 位于中心 |
| B | `2 → 4 → 1 → 1 → 4 → 2` | 2 µm/s 位于最外层，1 µm/s 位于中心 |
| C | `4 → 1 → 2 → 2 → 1 → 4` | 4 µm/s 位于最外层，2 µm/s 位于中心 |

三个 block 分别使用三个新的、相互可比的 `10 × 10 µm²` 区域；A/B/C 与区域的对应关系、三个 block 的实际先后顺序在实验开始前随机化并记录。每种速度在一次液体制备中因此有 6 张 map，即 384 条曲线；每个浓度共有 18 张 map，即 1152 条曲线。

区域之间不重叠，建议中心至少相隔约 15 µm，并预留明确的 sacrificial/anchor 区。block 内六张 map 必须回到同一 grid coordinates；block 间则故意换区域，以免一个 pixel 累积 18 次正式接触。

`8 × 8` 的曲线数足够做 map 内稳健空间统计。它的真正优势不是 64 条曲线本身，而是允许把时间用于 18 张平衡 map。速度处理的实验单位是整张 map；64 个 pixel 是空间配对样本，不能被当作 64 个独立 velocity replicates。

一个 A block 只有 384 条曲线，可以作为 pilot；它不能独立识别非线性 drift、区域差异和 carryover，不能支撑正式的零速外推。正式实验最低要求为 A/B/C 三个 block。若目标是对 glycerol 浓度趋势作可推广的结论，还需要独立制液/独立实验日重复，见第 4 节。

## 2. 为什么必须重做实验设计

当前 20-08-26 数据不是简单的“hydrodynamic force 太大”，而是 velocity 与多个系统效应叠加：

1. 原三速 map 基本按 `2 → 1 → 4 µm/s` 连续获得，速度与采集时间、表面接触历史完全混杂。
2. 所有 5 个“后测速度更慢”的比较都给出正 hydrodynamic amplitude，而所有 9 个最后测 4 µm/s 的比较都给出负 amplitude。物理 no-slip force 不应仅因顺序改变而反号，因此存在强 map-history 分量。
3. 30 wt% 的同速对照仍有约 `27.2 pN` 的线性 nuisance residual；这给出了当前系统在一张 map 到下一张 map 之间的实际误差尺度。
4. 30 wt% 曾同时出现约 6.7 nN 与 18.5 nN 两个 load regime；它们不能作为 velocity 对照混合。
5. contact InvOLS、terminal load、snap-in、contact height 与 far-field slope 都会随 map 或 map 内采集顺序变化。
6. 当前 retract 在主要拟合区间没有足够的 free non-contact support，且部分曲线发生 detector clipping；因此不能把已黏附 retract 强行用作 hydrodynamic 反号检验。
7. 10 wt% 的 load、硬接触有效率和信号强度属于另一测量状态，所以它只保留为历史 QC，不进入主趋势。

这些判断的本地数值依据分别保存在 [velocity joint-fit report](velocity_joint_fit_results/REPORT.md)、[velocity/systematics report](velocity_systematics_results/REPORT.md)、[surface-force report](surface_force_results/REPORT.md) 和 [cantilever calibration report](results/REPORT.md)。

回文设计利用时间对称性解决第一阶 drift。对同一 block、同一速度 `U`、同一物理 pixel `p`，定义

```text
F_sym,b,U,p(D) = [F_early,b,U,p(D) + F_late,b,U,p(D)] / 2
F_hist,b,U,p(D) = [F_late,b,U,p(D) - F_early,b,U,p(D)] / 2
```

在一个真正对称的速度序列中，两个同速 map 的时间中点关于 block 中心对称；`F_sym` 近似消除线性随时间变化的项，`F_hist` 直接测量该速度下的时间/接触历史。三个旋转序列又让 1、2、4 µm/s 各自在外层、中层和中心出现一次，从而避免“某个速度永远隔得最久”这一新混杂。

回文只能消除近似线性 drift，不能神奇地消除突跳、不可逆污染、任意非线性 aging 或人为干预。因此必须同时保留同速差、实际时间、接触次数、raster 方向和所有 QC 指标。

## 3. 物理量、标定与模型口径

### 3.1 原始信号到力和距离

所有正式分析都从 JPK archive 中的 raw int32 channel 及其 conversion metadata 重建，不使用文件内写入的 sensitivity 或导出的 force channel：

```text
delta = S_c [V_raw - V_baseline]
F = k delta
D = measuredHeight + delta - contactHeight_p
```

其中：

- `S_c` 是对当前浓度、当前独立实验日的全部合格硬接触曲线重新估计的共同 InvOLS/sensitivity；
- `k = 0.2969899087 N/m` 来自 cantilever 1 已完成的 thermal/contact calibration；
- `contactHeight_p` 对每个 pixel 独立确定，不能用一张 map 的统一 contact zero；
- 距离符号按当前 raw pipeline 的 scanner convention 固定，并用硬接触斜率和远场方向检查。

当前数据已显示 solution-dependent optical response：0、20、30、40 wt% 的共同 sensitivity 分别约为 66.609、60.718、58.673、57.903 nm/V。它们与最初独立硬基底标定的 85.30 nm/V 不同，说明不能把一个空气/另一液体中的 InvOLS 复制到所有溶液中。

### 3.2 实际速度

hydrodynamic force 应使用 gap-closing speed，而不是只使用软件输入的 nominal scanner speed：

```text
U_gap(D) = -dD/dt = -d[measuredHeight + delta]/dt
```

接近表面时 cantilever 会偏转，球的实际速度可低于 cantilever base/scanner 的速度。每条曲线应由 `measuredHeight`、`delta` 与原始 time channel 计算 `U_gap(D)`；加速、减速和 reversal 区域不进入拟合。

### 3.3 hydrodynamic 数量级

no-slip、小间隙球–面近似为

```text
F_hyd(D,U) = 6 pi eta R^2 U / D.
```

在 25.6 °C、`R = 4.546849 µm` 下，每 `1 µm/s` 的 nominal speed 对应的理论力如下；4 µm/s 时将表中数值乘 4。近接触处仍应改用实际 `U_gap(D)`。

| glycerol wt% | η (mPa·s) | 25 nm | 50 nm | 100 nm | 150 nm | 200 nm |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.8806 | 13.73 pN | 6.86 pN | 3.43 pN | 2.29 pN | 1.72 pN |
| 20 | 1.5013 | 23.40 pN | 11.70 pN | 5.85 pN | 3.90 pN | 2.93 pN |
| 30 | 2.0886 | 32.56 pN | 16.28 pN | 8.14 pN | 5.43 pN | 4.07 pN |
| 40 | 3.0796 | 48.00 pN | 24.00 pN | 12.00 pN | 8.00 pN | 6.00 pN |

例如在 100 nm，4 与 1 µm/s 的理论差只有约 10 pN（0 wt%）到 36 pN（40 wt%）；它与当前约 27 pN 的同速 history floor 同量级。这正是需要 map-level replication、回文顺序和高黏度 pilot 的原因。

### 3.4 electrostatic 主模型

主模型采用用户指定的“相同表面、共同 potential”模型：球和面均为 silica，在每个浓度中共享同一个 surface-potential magnitude。对无量纲电势 `u=e psi/(k_B T)` 与 `H=kappa D`，平板问题为

```text
u'' = sinh(u),
u(+/-H/2) = u_s,
u'(0) = 0.
```

由中面对称压力积分后，用 sphere–plane Derjaguin prefactor `2 pi R` 得到 `F_EDL`。该几何不是 equal-sphere：不能使用 `R/2` 的两等球 prefactor。观察模型至少包括

```text
F_obs = g_m [F_EDL(D; |psi_s|, lambda_D)
             + F_vdW(D; A_H)
             + chi F_hyd(D, U_gap)]
        + a_m + b_m(D-D_ref) + time/history terms + error.
```

主分析沿用 `A_H = 2.4e-21 J`，但必须把 Hamaker constant 随混合物组成变化作为系统敏感性测试，而不能把固定值称为已测事实。相同 silica 的 normal-force 数据只确定电势幅值；负号来自 silica chemistry。更准确的报告名称应是“equal-surface PB boundary-potential magnitude”。只有额外假设该 boundary potential 等于 slipping-plane potential 时，才报告 `zeta = -|psi_s|`。

## 4. 独立重复与浓度顺序

### 4.1 推荐正式规模

推荐做 4 个独立实验 session，每个 session 使用独立制备的 0、20、30、40 wt% aliquot，并用四序列 Williams design 平衡浓度在一天中的位置和一阶 carryover：

| session | 浓度顺序（wt%） |
|---:|---|
| 1 | `0 → 20 → 40 → 30` |
| 2 | `20 → 30 → 0 → 40` |
| 3 | `30 → 40 → 20 → 0` |
| 4 | `40 → 0 → 30 → 20` |

这比永远从水测到高 glycerol 更重要。若资源只允许 3 个 session，可以完成前三个作为最低探索性重复，但它不能完全平衡每个浓度的日内位置与 carryover，结论必须保留这一限制。

这个 crossover 设计假定换液后的 carryover 足够小。应先做一次 washout check：测水基线，暴露于 40 wt%，按正式流程冲洗回水，再比较 conductivity/refractive index、sensitivity、adhesion 与低速 force curve。若回水后仍明显不同，说明同一 probe 的化学 history 不可逆；此时 Williams design 只能平衡、不能消除 carryover，主结论应限定为该 probe 的连续暴露实验。若要推广为 silica 材料的一般浓度效应，需用额外独立制备的 silica colloidal probes 重复完整浓度设计，并把 probe identity 作为最高层 random effect。

### 4.2 计数层级

正式 4-session 设计的层级为：

| 层级 | 每个浓度、每个 session | 四个主浓度、四个 session |
|---|---:|---:|
| pixel curves | 1152 | 18,432 |
| maps | 18 | 288 |
| palindrome blocks/regions | 3 | 48 |
| independent solution/session replicates | 1 | 16 concentration-session units |

一张旧 `16 × 16` map 有 256 条曲线；三张旧 map 共 768 条。新的 18 张 `8 × 8` map 共 1152 条，所以曲线总量没有减少，同时 map-level velocity replication 从每速约 1 张提高为每速 6 张。

进一步说，每个浓度、每个 session、每个速度的 6 张 map 组成 3 个 time-symmetric pairs，分别来自 3 个区域；四个 session 后，每个浓度、每个速度共有 24 张 map、12 个 symmetric pairs。速度效应最终依赖这 12 个 pair 和 4 个独立 session，而不是依赖 `24 × 64` 个彼此独立的点。

### 4.3 10 wt% 的处理

- 10 wt% 不进入主设计、主拟合、主浓度趋势或 sample-size 计算。
- 若需要追查旧异常，只在所有主实验结束后，以完全相同 load 和完整 A/B/C 设计单独重测；在看见结果前就标记为 `diagnostic_only`。
- 旧 10 wt% 不因结果异常被删除；它作为 different-measurement-regime 的 QC 证据保留。

### 4.4 时间预算

当前 1 µm Z travel、approach/retract 同速时，每条曲线仅 ramp 的理论时间为：1 µm/s 时 2.0 s、2 µm/s 时 1.0 s、4 µm/s 时 0.5 s。64 个 pixel 对应的无 overhead 下限分别约 2.13、1.07、0.53 min。

由旧 16×16 数据的实际起始时间估算，8×8 map 加入 XY movement、feedback、写盘和切换后约为 2.5–3.5、1–1.5、0.6–1 min。一个六-map 回文 block 约需 9–12 min；三个 block 连同 anchor curves 约 35–50 min。加上换液和不少于 20 min 的 equilibration，一个 concentration-session 计划 60–80 min；四个浓度约 4–6 h。正式四-session 设计约需 16–24 h 上机时间，另加 40 wt% pilot。实际首个 pilot 后应以 JPK 记录的 start/end time 更新排程，而不是继续使用这个估算。

## 5. 样品、液体和环境准备

### 5.1 silica 表面与 colloidal probe

1. 球与平面使用同一类 silica，并沿用既定 cleaning SOP；同一研究内不改变清洗化学品、时间、plasma 条件或清洗到测量的等待时间。
2. 每个 session 使用新鲜、可追踪的 planar silica 区域；记录 substrate ID、批次、清洗时间和安装方向。
3. 主 probe 使用 cantilever 1 上的同一 colloidal sphere。记录 laser spot、photodiode sum、cantilever mounting 和液池重装；一旦重装，视为新的 calibration epoch。
4. 现有 reverse imaging 给出局部半径 `4.546849 µm`，拟合窗口变化给出约 `3.980307–4.546849 µm` 的 form/window spread；主拟合固定前者，系统敏感性覆盖整个范围。因为 `F_EDL ∝ R` 而 `F_hyd ∝ R²`，radius 对 hydrodynamic amplitude 尤其关键。
5. 现有局部球冠残差 RMS 约 1.379 nm、相关长度约 72.9 nm；这些是 apparent local topography 指标，不等于材料本征 roughness。实验前后做相同方法的 reverse image 或至少光学/接触 QC，检查球是否污染、脱胶或换了接触 asperity。
6. level planar sample，使大尺度 tilt 最小；不得通过一个全图统一 contact zero 来“校正”残余 tilt，而应保留每 pixel 的 contact height。

同一 substrate 上即使换到未接触区域，整个平面也已暴露于前一种溶液，所以“新区域”只能减少 mechanical contact history，不能消除 chemical exposure history。每个浓度使用不重叠的 region bank，并由第 5.3 节的换液 QC 与第 4.1 节的 balanced order 处理这一限制；若样品架允许在不重装 cantilever 的情况下使用独立 silica coupons，则优先为每个浓度配置独立 coupon，并在 session 间轮换 coupon–concentration 对应关系。

### 5.2 水–甘油溶液

1. 按质量法配制 0、20、30、40 wt%，使用可追踪天平和密闭容器；记录水、glycerol 的 lot、纯度、称量值、配制时间和 batch ID。
2. glycerol 吸水且开口液滴会改变质量分数。配液、转移和液池暴露时间应固定并尽量缩短；容器和液池在允许条件下密闭。
3. 缓慢混匀并充分 equilibration；避免剧烈振荡引入气泡。上机前检查可见颗粒和气泡，必要时采用与所有浓度一致的过滤/脱气流程。
4. 每个 concentration-session 的新鲜 aliquot 在测量前后各留一份。至少记录 refractive index 或 density 作为质量分数的独立 QC；若有条件，测 bulk viscosity。
5. 不外加盐不等于 ionic strength 为零。痕量离子、silica 释放的 counterions 和空气 CO₂ 会控制 screening。记录每份液体测量前后的 conductivity 与 apparent pH；它们不是 PB 拟合的强制输入，但用于识别 batch drift 和解释 `lambda_D`。混合溶剂中的 pH electrode 与 conductivity-to-ionic-strength conversion 都有介质依赖，因此不得把它们无校准地当作绝对 thermodynamic pH 或精确 ionic strength。
6. literature viscosity、permittivity 与 density 必须按实际 mass fraction 和 25.6 °C 插值。温度误差同时改变 viscosity、dielectric response、surface charging 与 detector drift。

### 5.3 换液、carryover 与 equilibration

1. 不仅按固定体积冲洗，还要按结果确认换液完成：至少冲洗 5 个有效 cell volume，并继续到 effluent refractive index/density 或 conductivity 与 incoming aliquot 在仪器精度内一致。
2. 避免冲洗流直接冲击 cantilever 或引入气泡；换液后重新检查 laser alignment、sum signal 和自由 cantilever baseline。
3. laser 打开、液池装好后固定 equilibration 时间，建议不少于 20 min；正式开始条件为 probe 附近实测温度 `25.6 ± 0.2 °C`，且 5 min 内变化不超过 0.1 °C。
4. 温度传感器应尽可能接近 probe，而不是只记录 room temperature。已有 drainage-force 文献显示 laser/electronics 可使 cell temperature 明显高于室温，且 ±1 °C 足以给 viscosity 带来约 5% 量级不确定度。
5. 每次换液后先在 sacrificial 区做少量低载荷曲线，确认无气泡、无异常黏附、无 detector saturation，再进入正式区域。

## 6. 仪器与 acquisition 参数

### 6.1 在正式实验前一次性冻结的参数

| 参数 | 推荐设置/原则 | 原因 |
|---|---|---|
| grid | `8 × 8` | 把时间转化为 map-level replication |
| field of view | 保持 `10 × 10 µm²` | 与旧数据空间尺度一致；pixel pitch 约 1.25 µm |
| back-and-forth | 保持开启并记录 | 减少横向回程；分析时恢复物理 x |
| slow-scan direction | 同速 early/late 尽量相反 | 分离 y 空间梯度与采集时间梯度 |
| Z start/range | 保持当前约 1.0 µm，并确认最初至少 700–900 nm 为自由远场 | 支撑 baseline、cantilever-drag 和 PB 远场检查 |
| approach speeds | 1、2、4 µm/s | 与旧数据可比；回文内只在 map 间改变 |
| retract speed | 与对应 approach 相同 | 保留 branch-symmetry 诊断的可能性 |
| contact pause | 0 s | 减少 dwell/creep/conditioning；全程一致 |
| force setpoint | pilot 后冻结为一个 target force | 防止再次产生 6.7/18.5 nN strata |
| closed-loop Z | 开启，保存 `measuredHeight` | 使用实测而非 command height |
| raw channels | vDeflection、measuredHeight、height、time/seriesTime、hDeflection、error，若可用再存 detector sum | 允许 raw reconstruction 和 cross-talk QC |
| filter/bandwidth | 全速度固定并记录 | 防止速度与 electronics lag 混杂 |
| automatic baseline adjust | pilot 决定后全程冻结；记录是否每 line 调整 | 当前文件显示该功能开启，可能引入 row/map offset |

旧 raw header 中 `force-baseline-adjust` 为 `enabled=true`、`interval=1`、`beginOfLine=true`，并在 deadtime 后平均 100 samples。新实验若继续使用，必须保持完全相同；若 pilot 决定关闭，则所有正式 session 从第一张开始都关闭，不能在正式数据中混用两种状态。

### 6.2 采样率问题

当前设置在 1 µm Z travel 上每个 segment 固定 1000 points，因此 spatial step 约 1 nm，但 acquisition rate 随速度从约 1、2 到 4 kHz 变化。若 JPK 允许，推荐固定时间采样率和 analog bandwidth，例如统一 4 kHz：

- 1 µm/s：4000 points/segment；
- 2 µm/s：2000 points/segment；
- 4 µm/s：1000 points/segment。

8×8 后文件量仍可接受。若软件不能这样设置，则保留 1000 points，但必须在硬表面 pilot 中测定 electronics/filter delay：固定硬接触区域分别用 1、2、4 µm/s，检查 contact response 和 contact position 是否出现与 `U` 成比例的位移。一个固定时间延迟 `tau` 会产生 `U tau` 的假 separation shift，并直接伪造 snap-in 和 `lambda_D` 的速度依赖。

### 6.3 force setpoint 的确定

不能直接沿用旧数据中任一个 load。先在 40 wt% pilot 用低到高的小步调整，选择“能够在至少 90% 曲线中得到 ≥20–30 nm 线性硬接触段，同时不出现明显表面损伤、异常 adhesion 增长或 detector saturation”的最低 load，然后锁定到所有浓度与 session。

对 `k = 0.29699 N/m`，20 与 30 nm cantilever deflection 约对应 5.94 与 8.91 nN，因此约 9–12 nN 是合理的 pilot 起点，而不是预先宣布的最终值。最终 target 必须由 pilot 的硬接触有效率决定。正式 block 中不得调 setpoint；post hoc map median terminal load 必须保持在 target 的 ±5% 内。

JPK 实际控制量若仍是 voltage setpoint，应在每个浓度开始时用当时的 online sensitivity estimate 计算 `V_set ≈ F_target/(k S_online)`。这只用于把仪器工作点放到相同 force 附近；最终 load 仍由全部 raw hard-contact 数据得到的 `S_c` 重建，不能使用软件内写入的 sensitivity 作为最终标定。

### 6.4 raster 方向

若软件支持反转 slow-scan direction 而不改变物理坐标，每个六-map block 使用交替方向：

```text
map position:       1  2  3  4  5  6
slow-scan direction ↑  ↓  ↑  ↓  ↑  ↓
```

这样每一对 time-symmetric 同速 map 的 slow-scan direction 相反。另一些 block 可整体反转为 `↓ ↑ ↓ ↑ ↓ ↑`，并把方向作为模型因子。如果仪器不能可靠反向，则所有 map 保持同一方向，但必须明确：contact/far-slope 的 y-gradient 与 map 内时间仍不可完全区分。

## 7. 每个 concentration-session 的逐步操作

### 7.1 正式开始前

1. 核对 solution batch、substrate、probe、cantilever、温度传感器、laser spot 与数据目录。
2. 换液并满足第 5.3 节的 composition 与 temperature 稳定条件。
3. 在 sacrificial 区完成 force-setpoint、saturation 与 hard-contact 快速检查；正式实验开始后不再调参数。
4. 预先生成三个 block 的区域分配、A/B/C 执行顺序、第一张 map 的 slow-scan direction；将计划写入 metadata CSV，不能根据中途曲线“看起来好不好”改顺序。
5. 在五个新的 sacrificial 坐标采集第一组 5 条独立 force curves，标记为 sensitivity/QC anchor group `S0`。

### 7.2 三个 main blocks

对每个 block：

1. 移动到一个新的 `10 × 10 µm²` 区域，保存区域中心坐标和一张低扰动 reference/contact-height map（若额外 reference map 会引入接触，则使用第一张正式 map 的 contact-height）。
2. 连续执行分配给该区域的六-map 回文序列。只允许在 map 之间改变预定速度与 slow-scan direction，不重新调 laser、不改 setpoint、不加液、不移动区域、不人工重做“难看”的单个 pixel。
3. 每张 map 完成后立即保存，不覆盖旧文件；记录 nominal speed、实际 start/end time、温度、setpoint、scan direction、filter、point count、异常提示和 operator intervention。
4. 只做不改变采集决策的 quick QC：文件可读、64 个 index、无 approach saturation、实际 scanner speed 合理、terminal load 未失控。不得在中途根据 surface-force 方向选择性重复。
5. 六张完成后，在五个新的 sacrificial 坐标采集一组 5 条独立 force curves，作为 `S1`、`S2` 或 `S3`。

完成三个 block 后，固定等待 10 min，不接触正式区域；再在五个新的 sacrificial 坐标采集 `S4`。这样每个 concentration-session 有 25 条独立 force curves，分布在五个时间点。它们首先用于 sensitivity drift；若表面、load、速度和距离范围与正式条件一致，也可作为独立-force secondary cohort，但不能与 map pixel 简单合并计数。

### 7.3 block 被中断时

- laser realignment、solution addition、bubble removal、setpoint change、software restart、区域丢失或超过预定暂停时间都构成 block interruption。
- 已中断 block 不与后续 map 拼成一个“完整”回文。保留原文件并标记原因；修复后在新区域从六张序列的第一张重新开始。
- 不因为物理结果不符合预期而重测；只有预注册的 operational failure 才允许重测。

## 8. 40 wt% pilot 中的额外诊断

40 wt% 给出最大的 hydrodynamic signal，最适合先识别 instrumentation/history。正式大规模采集前完成以下 pilot；这些结果用于冻结设置，不进入之后的 blind primary concentration comparison。

### 8.1 单个 A block

完成 `1 → 2 → 4 → 4 → 2 → 1 µm/s` 的 8×8 block，检查：

- 64 个 pixel 完整度与 hard-contact validity；
- target load 和 saturation；
- 反向 raster 是否能恢复同一物理 grid；
- actual `U_gap(D)` 与 nominal speed；
- same-speed pair 的 contact slope、contact height、snap detection、far slope 与 force residual；
- 4–1 µm/s 的力差是否至少与 same-speed floor 可比。

正、负或零 hydrodynamic amplitude 都不是 pilot 的通过/失败条件；pilot 只判断仪器与设计是否可执行，不能用期望物理结果筛选参数。

### 8.2 time-only 与 contact-dose control

回文仍把 elapsed time 与累计接触次数部分绑定。用两个匹配的新区域做：

1. **time-only：** `2 µm/s map → 等待与四张中间 map 相同的实际时长、期间不接触 → 2 µm/s map`；
2. **contact-dose：** 在另一新区域连续做六张 `2 µm/s` map。

比较前后 same-pixel 的 `F(D)`、InvOLS、load、contact height、snap 和 far slope：

- 两种控制都变化：以 passive time/thermal/chemical drift 为主；
- 只有连续接触明显变化：以 contact conditioning、污染、磨损或 adhesion history 为主；
- 变化与 row/scan direction 相联：以 raster-time、piezo creep 或 spatial gradient 为主。

至少在第一个正式 session 做一次；若 probe 被重装、清洗流程改变或新 session 的同速 drift 明显增大，则重复。

### 8.3 推荐的 free non-contact hydrodynamic control

当前 contact maps 的 retract 在 25–250 nm 内没有足够 free support。若仪器允许在不触发 snap/contact 的安全距离反转，另做一组同样回文的 non-contact 8×8 maps，turnaround 约在 150–200 nm，并覆盖到至少 800 nm：

- approach 与 retract 采用相同速度；
- 不用这组曲线估计 contact zero 或 sensitivity；距离零来自紧邻的 contact reference map；
- 在相同 `D` 上，equilibrium surface force 应同号，而 hydrodynamic contribution 在 approach/retract 之间反号；
- 只有当目标区间 ≥80% pixel 同时具有未黏附、未 clipping 的 approach/retract 数据时，才把 odd/even 分解用于正式 hyd validation。

如果 non-contact turnaround 无法可靠设置，则不以已黏附 retract 代替；主实验仍可完成，但 hydrodynamic identification 的证据等级降低。

### 8.4 区分 hydrodynamics 与真正的 rate-dependent surface force

velocity 除了产生 drainage force，也可能改变 interfacial charge regulation、ion redistribution、preferential glycerol/water structure 或 snap instability。它们是真正的动态表面过程，不能简单归入仪器 drift。建议在 30 或 40 wt% 的 pilot 中另做少量 single-point force spectroscopy：在安全的非接触 separation（例如约 200、100、50 nm；若会 snap 则停止在更远处）加入固定 0.5–2 s hold，记录 deflection 是否随 hold time relax。每个 separation 至少在多个新点重复，并用紧邻的 contact reference 确定 D。

- 若 hold 中没有可分辨 relaxation，而速度差符合 `eta U/D`，支持 hydrodynamic 解释；
- 若固定 D 下仍有明显 relaxation 或其时间常数随浓度强变，主模型需加入 interfacial kinetic term，不能把全部 velocity effect 外推成 `U=0` drainage correction；
- snap-in 本身不用于定义这一 relaxation，因为 snap 是力梯度与 cantilever stability 的动态阈值。

## 9. 预注册 QC 与排除规则

### 9.1 curve 和 map 级规则

| 项目 | 合格/处理规则 |
|---|---|
| raw archive | ZIP CRC 通过；segment/channel byte count 与 metadata point count 一致；保存 SHA-256 |
| approach saturation | 主 approach 任一关键 channel 触及 encoder limit 即排除该 curve；不得截平后继续拟合 |
| hard contact | robust linear contact span ≥20 nm，`R² ≥ 0.995`，斜率符号与 scanner convention 一致 |
| map completeness | ≥58/64 合格为 primary；51–57/64 进入 review；≤50/64 或失败成片聚集则重测/排除该 map |
| terminal load | map median 在冻结 target 的 ±5%；超出则整张 map 标记 load failure，不通过 rescaling 冒充同一 load |
| nominal scanner speed | 远离 reversal 的 measuredHeight median speed 在 nominal ±2%；超出则用 actual speed 并标记；严重非恒速 map 重测 |
| temperature | primary 为 25.4–25.8 °C；block 内范围不超过 0.2 °C。越界数据保留为 secondary，不隐瞒 |
| sensitivity drift | S0–S4 group medians 总漂移 ≤3% 合格；3–5% review 并在模型中放宽 force-scale uncertainty；>5% 时该 concentration-session 不作为 primary |
| XY registration | same-area maps 位移 <0.5 pixel 直接配对；0.5–1 pixel 仅在可靠 cross-correlation 后使用；>1 pixel 重测 block |
| detector sum/laser | 相对 S0 出现突变、明显趋势或重新对光时建立新 calibration epoch，不能跨 epoch 共用一个 sensitivity |
| retract clipping | 保留 censoring 状态；不把 detector floor 当作 pull-off force，不把未观察到 detach 当作零 |

这些阈值在看主 surface-force 拟合前应用。不能因为 `chi` 为负、`lambda_D` 不漂亮或浓度趋势不单调而删除 map。

### 9.2 不作为自动排除、但必须报告的诊断

- snap-in 未检测到；“未检测”不是 snap distance = 0。
- far-field slope 非零；它可能是 optical drift、cantilever drag、长程力或 baseline-control history，不能仅凭斜率大小删除。
- contact-height 平面或 row gradient；先用 raster reversal 判断空间还是时间。
- same-speed force difference；这是 history 的实测结果，不是坏数据的同义词。
- localized adhesion/outlier；除非 raw saturation 或预注册 instrument failure，否则保留并用 robust spatial summaries 表示。

### 9.3 block-level hydrodynamic 可识别性判据

只有同时满足下面条件，才把 `U=0` intercept 称为经过验证的 hydrodynamic extrapolation：

1. `F(D)` 对实际 `U_gap` 的方向在 A/B/C 和独立 session 中一致；
2. `chi` 为正，且与 no-slip order of magnitude 相容；不要求机械地等于 1，但不同 block 不能主要由正负抵消得到；
3. 速度项优于只含 intercept、linear-distance drift 和 time/history 的 nested model；
4. `beta(D)` 具有预期的近似 `eta/D` scaling，而不是简单常数或线性 D slope；
5. speed × acquisition-order、speed × raster-direction interaction 小于主 speed effect；
6. same-speed `F_hist` 小于 4–1 µm/s 的 hydrodynamic contrast，或被显式 history model 解释；
7. leave-one-block、leave-one-speed 和 leave-one-session 结果不反号；
8. 若 free non-contact control 可用，approach/retract 的 odd component 反号且量级相符。

若这些条件失败，正确输出不是强行减掉 `6 pi eta R²U/D`，而是：报告各速度的 model-conditioned PB 参数、慢速结果和 history systematic range，并明确 `U=0` 未被实验识别。

## 10. 数据处理与统计分析计划

### 10.1 raw reconstruction

1. 复制原始 archive 到只读目录；建立 full relative path、文件大小、timestamp、SHA-256、ZIP CRC 和 acquisition metadata manifest。
2. 对每个 channel 使用其自己的 scaling/conversion chain 解码；不要由另一个文件或 JPK export 猜单位。
3. map 的实际 start/end/midpoint 由各 pixel header 与 time/seriesTime 重建；不能盲信 top-level `end-time`。当前 raw map 的 top-level start/end 曾相同，显然不能代表真实十余分钟采集时长。
4. 保留 raw、constant-referenced 和 line-corrected 三种数据视图，不覆盖原数据。
5. 初始 20% samples 作为候选 far field，至少 80 points；用 robust line `V_ff(h)=a h+b` 和 MAD rejection 估计 slope、R²、residual noise。
6. calibration/contact/event branch 可减去完整 far-field line；hyd/PB 主分支只减 far-field constant，并把 `a_m+b_mD` 与 physical `1/D` term 一起拟合。这样不会在拟合前把可能的 hydrodynamic curvature 当 drift 剪掉。

### 10.2 sensitivity

1. 每条合格 approach 的终端硬接触拟合 `V = a + b h`，得到 `S=1/|b|`。
2. 先在每个 raw source/map 内形成 robust anchor，再在同一 concentration-session 内以 source/map 等权聚合；不能让一张 64-pixel map 因曲线多而压倒 5 条独立曲线组成的 anchor group。
3. 最终 force reconstruction 使用该 concentration-session 的共同 `S_c`；map-specific contact response 只作为 QC 或 observation-scale nuisance `g_m`，不逐 map 自由重标定力，否则会把真实漂移和速度效应一起归一化掉。
4. embedded JPK sensitivity 和 embedded force 仅作对照，绝不替代 raw estimate。

### 10.3 map 几何与配对

1. 根据 metadata 读取 grid size、坐标、reflect/back-and-forth，不根据文件 index 猜空间顺序。
2. serpentine 的奇数 acquisition row 恢复到物理 x；反向 slow-scan 的 map 恢复到共同物理 y。
3. 用 contact-height/topography map 做 translation registration；不旋转或非线性扭曲数据来追求更好 force collapse。
4. 同一 block 的 force curves 插值到共同 D grid 后做 same-pixel early/late 配对；插值只在两条曲线共同支持的区间内，不外推。
5. 若保持 10×10 µm²，8×8 数据用 `2 × 2 pixel` spatial blocks，可得到 16 个约 `2.5 × 2.5 µm²` 的 block，与旧 16×16 数据用 4×4 pixel blocks 的物理尺度一致。

### 10.4 先做非参数分解，再做联合拟合

每个浓度、每个 session 的分析顺序固定为：

1. 画每张 map 的 8×8 far-field slope、contact InvOLS、terminal load、contact height、snap detection/distance 与 selected-D force maps；不做视觉 clipping。
2. 对每个同速 pair 计算 `F_sym` 与 `F_hist`，并画随 D、pixel、row/order 的变化。
3. 对 A/B/C 的 `F_sym` 用 actual `U_gap` 拟合每个 D bin 的 slope/intercept，同时含 block、region、map midpoint time、raster direction 与必要的 quadratic-time term。
4. 以 `6 pi eta R² U_gap/D` 为已知形状拟合 `chi`，同时允许每张 map 的 constant/linear-distance nuisance；比较 M0（无 hyd）、M1（`chi=1`）、M2（自由 `chi`）。
5. 只有通过第 9.3 节时才使用 `U=0` force curve 作为 primary equilibrium estimate；否则以三个速度分层结果和慢速结果为主。
6. 对 equilibrium force 拟合 equal-surface nonlinear PB sphere–plane model，主窗口预注册为 25–250 nm。

### 10.5 必做的 sensitivity analyses

- fit windows：25–150、30–250、25–350 nm；
- contact-zero shift：至少 ±2.5 nm；
- sphere radius：3.980307–4.546849 µm；
- force scale：由 concentration-session sensitivity repeatability 与 cantilever-k uncertainty 传播；
- dielectric constant、viscosity：按温度与质量分数的不确定度传播；
- `A_H`：固定主值与合理 composition-dependent range；
- PB boundary：equal constant potential 为主，equal constant charge/charge-regulation 只作模型敏感性；
- hydrodynamics：no-slip 主模型、有效 wall shift/roughness 或 slip 只作诊断；
- leave-one-map、leave-one-block、leave-one-speed、leave-one-session；
- constant-referenced 与 line-corrected preprocessing 双分支；
- 含/不含 independent-force secondary cohort。

### 10.6 统计层级与不确定度

模型与 bootstrap 的 resampling 层级必须是：

```text
session/independent solution preparation
  -> palindrome block / physical region
    -> map
      -> 2×2 spatial block
        -> pixel / separation bins
```

浓度和 actual velocity 是 fixed effects；session、region/block 与 map 是 random/group effects。先 resample session，再 resample block，最后才 resample spatial blocks。不能以 64 个 pixel 或数十个 D bins 的 naive standard error 作为浓度趋势的置信区间。D bins 来自同一条 curve，彼此相关。

每个浓度先独立得到 `lambda_D` 与 `|psi_s|`，再做浓度趋势；不预设必须单调，不用一个高阶 polynomial 强迫平滑。若 session-to-session spread 大于局部 Jacobian CI，以实验层级 spread 为主要不确定度。

## 11. 完整问题清单、诊断和处置

| 类别 | 可能问题 | 在数据中的表现 | 预防/诊断 | 对结论的影响 |
|---|---|---|---|---|
| 设计 | velocity–time confounding | 后测 4 µm/s 系统反号 | A/B/C 回文、actual timestamps | 未解除则不能零速外推 |
| 设计 | 非线性 drift/突跳 | early/late 平均仍不一致 | same-speed `F_hist`、quadratic time、interruption log | 给出 history bound，不强行修正 |
| 设计 | map duration 随速度不同 | 相同 ordinal position 不等于相同时间 | 回文 duration 对称；使用 midpoint time | 不用 map index 代替时间 |
| 设计 | 接触次数与 elapsed time 混杂 | 同速连续 map 漂移 | time-only/contact-dose controls | 区分 passive drift 与 conditioning |
| 设计 | pseudoreplication | pixel CI 极窄但 map 间反号 | map 为 velocity unit，session 为 concentration unit | 防止虚假高显著性 |
| 设计 | 浓度顺序/日内时间混杂 | glycerol trend 等同于一天中的 drift | 4-session Williams design | 未平衡时趋势仅为 conditional |
| 设计 | 换液 carryover | 高到低/低到高不对称 | effluent composition QC、平衡顺序 | 可伪造浓度趋势 |
| 设计 | 单一 colloidal probe | session 很多但 probe chemistry/roughness 不变 | 额外 probe confirmatory replicate | 无额外 probe 时只推广到该 probe/preparation |
| 空间 | 表面 heterogeneity | 同一 map 的 spatial plane/patch | 同区配对、三区域、2×2 blocks | 作为 region variance |
| 空间 | y-position 与 row-time 混杂 | row gradient 随 map 改变 | slow-scan reversal | 不能解除时不解释为空间结构 |
| 空间 | XY drift/misregistration | same-pixel 差呈边缘/feature pattern | closed loop、contact-height registration | >1 pixel 不做配对 |
| 表面 | silica aging/charging kinetics | 随浸泡时间 zeta/adhesion 漂移 | 固定 equilibration、time covariate、balanced order | 可改变 PB 参数本身 |
| 表面 | 接触污染/磨损 | adhesion、InvOLS、snap 随 contact dose 变 | 最低可用 load、新区域、dose control | 不可当 hydrodynamics |
| 表面 | sphere asperity/roughness | apparent wall shift、局部 radius 改变 | reverse image、radius/roughness sensitivity | hyd `R²` 特别敏感 |
| 表面 | sphere/plane 名义同材但状态不同 | 相同 potential 模型残差有系统结构 | 同 cleaning/aging；alternative boundary sensitivity | “共同 zeta”是模型假设 |
| 表面 | plane tilt/local topography | contact-height 平面可达几十 nm | per-pixel contact zero | 统一 zero 会偏 `lambda_D` |
| 液体 | mass-fraction error/吸水/蒸发 | viscosity、RI、force 随时间漂移 | gravimetry、密闭、RI/density before/after | 影响 hyd 和浓度轴 |
| 液体 | 混合未平衡、temperature gradient | map 内慢变或 convection | ≥20 min equilibration、near-probe T | 可产生 baseline/history |
| 液体 | laser heating | cell T 高于 room T | probe-near thermometry | η 对 T 敏感 |
| 液体 | trace ions/CO₂ | `lambda_D` batch-to-batch 大变 | conductivity/apparent pH before/after | no-added-salt 不等于 fixed ionic strength |
| 液体 | 气泡/颗粒 | 局部 jump、超大 adhesion/drag | 统一脱气/过滤、显微检查 | 受影响 block operational failure |
| 液体 | near-surface composition/viscosity | hyd 不按 bulk η scaling | 跨浓度 `chi`、non-contact control | bulk no-slip model 可能失效 |
| 液体 | dynamic charge/solvent relaxation | 固定 D 下 force 随 hold time 变 | non-contact hold spectroscopy | rate effect 不可全归为 hydrodynamics |
| 液体 | electroviscous coupling | D 约为几倍 `lambda_D` 时 drainage 偏离 `eta/D` | 检查 D/λ scaling、PB–flow model sensitivity | PB 与 hyd 未必简单相加 |
| 液体 | buoyancy/gravity 随 density 改变 | 浓度相关但近似 D-independent offset | 每液体独立 far reference | 不应误作长程 EDL |
| 标定 | 使用 embedded sensitivity | 不同浓度力尺度错误 | raw hard-contact consensus | 禁止作为正式 force scale |
| 标定 | sensitivity 随 RI/laser spot 变 | contact slope 随浓度/时间改变 | S0–S4、map contact response | >5% 新 calibration epoch/排除 |
| 标定 | spring constant 错误 | 所有 force 同比例偏移 | 固定 cantilever-1 calibrated k 与不确定度 | 影响 force、主要影响 zeta amplitude |
| 标定 | contact substrate compliance | contact slope 非 1 或 load-dependent | silica hard-contact、load sweep pilot | S 偏差不能由更多 pixel 消除 |
| 仪器 | detector saturation/nonlinearity | encoder floor、flat top | 低 load、保存 raw limit flags | saturated 点不拟合 |
| 仪器 | baseline auto-adjust | line/map 开头 offset 跳变 | pilot 冻结设置、记录每-line behavior | 作为 nuisance，不能暗中变化 |
| 仪器 | laser/photodiode drift | far slope、sum、InvOLS 同时变 | sum/QPD log、S0–S4 | 建立 calibration epochs |
| 仪器 | vertical/lateral cross-talk | hDeflection 与 vDeflection 同步 row pattern | 保存 hDeflection、检查 correlation | 可能伪造 contact response |
| 仪器 | piezo creep/hysteresis | contact height 与 row/time 相关 | measuredHeight、slow-scan reversal、settling | 用实际 height；保留 time term |
| 仪器 | actual speed 非 nominal | hyd amplitude 错误、近接触减速 | `dD/dt` 分离计算 | nominal speed 只作标签 |
| 仪器 | finite filter/electronics delay | snap/contact shift 近似 `U tau` | fixed bandwidth、硬表面多速 pilot | 可伪造 velocity dependence |
| 仪器 | cantilever-body drag | far-field velocity offset/slope | non-contact branch、map baseline nuisance | 不等于 sphere–plane squeeze force |
| 仪器 | vibration/line noise/aliasing | 不同速度的 PSD 或 residual frequency 不同 | 固定 bandwidth、保存 time channel、比较 residual PSD | 可伪造速度相关噪声与 event shift |
| 几何 | probe glue/tilt/off-axis torque | lateral signal、contact slope 和 drag 随方向变 | 固定 mounting、保存 hDeflection、前后成像 | 需要新 calibration epoch 或 geometry sensitivity |
| 数据 | far-field 含长程 EDL/hyd | baseline 拟合删除真实信号 | 1 µm range、constant-only 主分支 | 线性扣除只作 sensitivity branch |
| 数据 | contact zero 不确定 | `lambda_D` 与 zeta 同时偏 | per-pixel zero、±2.5 nm test | 局部 CI 不含此系统误差 |
| 数据 | snap detection selection | 只统计 detected 会改变样本组成 | detection fraction + censored/missing | 不把 missing 设为零 |
| 数据 | retract adhesion/clipping | 假 pull-off/假 sign reversal | explicit censoring、non-contact control | 不用于主 hyd 约束 |
| 数据 | smoothing/interpolation | jump 被抹平、correlation 被放大 | minimal filter、共同支持区插值 | 保存 raw/processed 双份 |
| 数据 | 文件名重复/metadata 错 | 错配浓度或重复计数 | path + SHA-256 + manifest | basename 不是 observation ID |
| 数据 | top-level map time 不可信 | start/end 相同或与 pixel duration 冲突 | 由 pixel header 和 time channel 重建 | 错误时间会破坏回文 drift 分解 |
| 模型 | no-slip 小间隙近似 | `chi` 偏离 1 | D/R 检查、exact/finite-gap sensitivity | `chi` 是 model-conditioned |
| 模型 | slip/effective wall/roughness | hyd 曲线像 shifted D | roughness/radius controls | 不直接宣称 molecular slip |
| 模型 | PB 1:1 point-ion 假设 | 系统残差、window dependence | nonlinear PB + boundary sensitivity | `lambda_D` 可是 effective screening |
| 模型 | charge regulation | constant-potential fit 随 D/window | constant-charge/CR sensitivity | zeta 不是唯一物性常数 |
| 模型 | non-DLVO hydration/solvation | <25 nm 残差和浓度依赖 | 主窗口从 25 nm 开始、改变 Dmin | 近场不用于主 PB 参数 |
| 模型 | Hamaker 随介质变化/retardation | 短程 attraction 拟合偏差 | fixed-primary + composition sensitivity | 不与 zeta/zero 同时自由漂移 |
| 模型 | 参数不可识别 | `lambda_D`、zeta、zero、drift、chi 高相关 | 独立 calibration、回文、profile/leave-out | 报 systematic range 而非窄 local CI |

## 12. 文件命名、日志和 provenance

推荐文件名：

```text
YYYYMMDD_C00_S01_BA_R01_M01_U1p0_DIRup_8x8.jpk-force-map
YYYYMMDD_C40_S01_S0_P01_U2p0.jpk-force
```

其中：

- `C00/C20/C30/C40`：glycerol wt%；
- `S01`：independent session/solution preparation；
- `BA/BB/BC`：palindrome type；
- `R01`：physical region；
- `M01...M06`：block 内 ordinal map；
- `U1p0/U2p0/U4p0`：nominal speed；
- `DIRup/DIRdown`：slow-scan direction。

每个文件在 metadata CSV 中至少有：

```text
relative_path, sha256, session_id, solution_batch_id, glycerol_wt_percent,
water_mass_g, glycerol_mass_g, substrate_id, probe_id, cantilever_id,
block_type, region_id, map_ordinal, nominal_speed, actual_scanner_speed,
start_time, end_time, midpoint_time, cumulative_contact_count,
grid_i, grid_j, ulength, vlength, scan_direction, back_and_forth,
z_range, approach_points, retract_points, filter_bandwidth,
setpoint_setting, reconstructed_terminal_load, temperature_start,
temperature_end, conductivity_before, conductivity_after,
apparent_pH_before, apparent_pH_after, refractive_index_or_density,
laser_sum, baseline_adjust_setting, operator_intervention, qc_status, notes
```

原始数据只读保存；derived CSV、plots、fit config、software version、random seed 与 artifact hashes 单独输出。任何排除都写明 raw identity、规则和原因。

## 13. 现场执行清单

### 开始一个 session

- [ ] 四个浓度 aliquot 独立制备、称量与 batch ID 完整。
- [ ] 当天的 Williams concentration order 已预先确定。
- [ ] probe、substrate、cantilever 1、`k` 和 radius identity 确认。
- [ ] laser、QPD sum、temperature probe、closed-loop scanner 正常。
- [ ] raw channel、bandwidth、point count、baseline-adjust setting 已冻结。

### 开始一个 concentration-session

- [ ] 换液 composition QC 达标，无气泡。
- [ ] probe-near temperature 达到 `25.6 ± 0.2 °C` 并稳定。
- [ ] sacrificial pilot 无 saturation，target load 正确。
- [ ] A/B/C 顺序、region 分配与 raster direction 已写入日志。
- [ ] S0 五条独立曲线完成。

### 每张 map 后

- [ ] 文件立即保存且未覆盖。
- [ ] 64 indices/segments 可读。
- [ ] 无 approach encoder saturation。
- [ ] nominal/actual speed、load、temperature、direction 已记录。
- [ ] 没有未记录的 operator intervention。

### 每个 block 后

- [ ] 六张 map 顺序完整；若中断则整 block 标记，不拼接。
- [ ] 相应 sensitivity anchor group 完成。
- [ ] contact-height、InvOLS、load、far slope 的 8×8 quick plot 生成。
- [ ] 只按 operational QC 决定是否重测，不看 PB/hyd 拟合方向。

### 每个 concentration-session 后

- [ ] 三个区域、A/B/C 共 18 张 map 完整。
- [ ] S0–S4 共 25 条独立 curves 完整。
- [ ] before/after solution QC 和 temperature log 完整。
- [ ] raw files 复制到只读位置并生成 SHA-256/CRC manifest。

## 14. 结果可以和不可以怎样表述

只有在 velocity-identification 判据通过后，才可写：

> 在 25.6 °C、给定 water–glycerol mixture 和 equal-silica sphere–plane 模型下，通过 balanced palindrome force mapping 识别并外推了 velocity-dependent drainage component，随后得到 effective screening length 与 equal-surface boundary-potential magnitude。

在 hydrodynamic 判据未通过但 PB 数值拟合稳定时，应写：

> 得到的是对 acquisition speed/history、contact-zero、calibration 与模型条件敏感的 effective screening length 和 boundary-potential parameter；现有实验未唯一识别零速 equilibrium force。

不得仅由 normal force 宣称电势符号；不得把 boundary potential 无条件称为 zeta；不得把 fitted `lambda_D` 无条件解释为由已知 bulk ionic strength 决定的 textbook Debye length；不得把负或小 `chi` 直接解释为 molecular slip，而忽略 roughness、wall shift、cantilever drag、filter lag 与 surface history。

## 15. 关键文献依据

- Ducker, Senden and Pashley, colloidal-probe AFM 的原始球–面直接力测量方法：[Nature 353, 239–241 (1991)](https://www.nature.com/articles/353239a0)。
- Guriyanova et al., sphere–plane drainage、实际球端速度、cantilever drag、温度/黏度与 roughness 导致 apparent slip 的实验处理：[Microfluidics and Nanofluidics 8, 653–663 (2010)](https://link.springer.com/article/10.1007/s10404-009-0498-2)。
- Honig and Ducker, hydrophilic particles 的 no-slip hydrodynamic boundary condition：[Phys. Rev. Lett. 98, 028305 (2007)](https://doi.org/10.1103/PhysRevLett.98.028305)。
- Maali, Wang and Bhushan, 使用多速度、大球和较硬 cantilever 测 no-slip drainage，并显式处理 electrostatic subtraction：[Langmuir 25, 12002–12005 (2009)](https://pubs.acs.org/doi/10.1021/la902934j)。
- Zhu, Attard and Neto，reliable colloid-probe hydrodynamic force protocol、污染与 boundary-condition 问题：[Langmuir 27, 6712–6719 (2011)](https://pubs.acs.org/doi/10.1021/la104597d)。
- Glycerol–water 无外加盐时 trace ions 仍控制 screening 的实例：[Microelectrophoresis of Silica Rods Using Confocal Microscopy](https://pmc.ncbi.nlm.nih.gov/articles/PMC5348103/)。
- glycerol 对 silica surface charging 的影响不能仅由 bulk dielectric constant 解释：[J. Am. Chem. Soc. 139, 15013–15021 (2017)](https://pubmed.ncbi.nlm.nih.gov/28972749/)。
- glycerol–water viscosity correlation：[Cheng, Industrial & Engineering Chemistry Research (2008)](https://doi.org/10.1021/ie071349z)。
- mixture dielectric properties：[Behrends et al., J. Chem. Eng. Data (2006)](https://pubmed.ncbi.nlm.nih.gov/16626219/)。
- equal-surface electrostatic interaction与 nonlinear PB/Derjaguin：[Hogg–Healy–Fuerstenau](https://pubs.rsc.org/en/content/articlehtml/1966/tf/tf9666201638)、[Stankovich and Carnie](https://pubs.acs.org/doi/10.1021/la950384k)、[Polat and Polat](https://doi.org/10.1016/j.jcis.2009.09.008)。

## 16. 2026-08-24 扩展记录：0–99.5 wt% 与 20–200 nm

本节记录后续研究范围的改变：最终目标从现阶段的 0、20、30、40 wt% 扩展到 `0–99.5 wt% glycerol`，并要求得到 `20–200 nm` 内的 equilibrium force。以下高浓度数值是实验设计用的 no-slip hydrodynamic 估算，不是现有数据已经测得的总表面力，也不是 PB equilibrium force。

### 16.1 固定 0.1/0.2 µm/s 不再是全浓度低-hyd 方案

在 25.6 °C、`R = 4.546849 µm` 下，采用 Cheng mass-fraction viscosity correlation 和球–面 leading lubrication force：

```text
F_hyd(D,U) = 6 pi eta R^2 U_gap / D

F_hyd[pN] = 3.89693 eta[mPa s] U_gap[um/s] (100 nm / D)
F_hyd(20 nm)[pN] = 19.4846 eta[mPa s] U_gap[um/s]
```

表中为 `U_gap = 0.1 µm/s` 的理论力幅值；`0.2 µm/s` 时所有力乘 2。

| glycerol wt% | η at 25.6 °C (mPa·s) | 20 nm | 50 nm | 100 nm | 200 nm |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.8806 | 1.72 pN | 0.686 pN | 0.343 pN | 0.172 pN |
| 20 | 1.5013 | 2.93 pN | 1.17 pN | 0.585 pN | 0.293 pN |
| 40 | 3.0796 | 6.00 pN | 2.40 pN | 1.20 pN | 0.600 pN |
| 60 | 8.6369 | 16.8 pN | 6.73 pN | 3.37 pN | 1.68 pN |
| 70 | 17.5529 | 34.2 pN | 13.7 pN | 6.84 pN | 3.42 pN |
| 80 | 43.8413 | 85.4 pN | 34.2 pN | 17.1 pN | 8.54 pN |
| 90 | 150.027 | 292 pN | 117 pN | 58.5 pN | 29.2 pN |
| 95 | 330.591 | 644 pN | 258 pN | 129 pN | 64.4 pN |
| 97.5 | 520.687 | 1.015 nN | 406 pN | 203 pN | 101 pN |
| 99 | 699.400 | 1.363 nN | 545 pN | 273 pN | 136 pN |
| 99.5 | 774.866 | 1.510 nN | 604 pN | 302 pN | 151 pN |

因此 `0.1/0.2 µm/s` 只能在低黏度区间按具体 force uncertainty 判断是否足够小；它们在接近纯甘油时绝不能称为 hydrodynamics-negligible。尤其在 99.5 wt%，0.1 µm/s 已给出约 `1.51 nN -> 151 pN`（20→200 nm），0.2 µm/s 为其两倍。

### 16.2 “可忽略”的速度阈值及其不可行性

以最严格的 20 nm 端点定义全窗口 hyd force 上限，在 99.5 wt% 时：

| 允许的 `F_hyd(20 nm)` | 最大 `U_gap` | 等价速度 |
|---:|---:|---:|
| 5 pN | 0.000331 µm/s | 0.331 nm/s |
| 50 pN | 0.00331 µm/s | 3.31 nm/s |
| 100 pN | 0.00662 µm/s | 6.62 nm/s |

将全窗口 hyd 压到 5 pN 以下会使完整 8×8 mapping 极慢，不能作为覆盖 0–99.5 wt% 的主策略。正式目标应从“把 hyd 降到零”改为“在可识别的多个 `eta U` 水平测量，并得到 `U -> 0` 截距”。

### 16.3 远场归零不会消除 hydrodynamic force

上述表格是相对于无限远的物理 no-slip 力。如果软件或 preprocessing 在有限参考距离 `D_ref` 设零，constant-referenced hyd 分量变为：

```text
F_hyd,referenced(D) = 6 pi eta R^2 U [1/D - 1/D_ref].
```

例如 `D_ref = 850 nm` 时，在 20、50、100、200 nm 分别仍保留原始 `1/D` 信号的约 `97.6%`、`94.1%`、`88.2%`、`76.5%`。有限距离归零只是减掉一部分 hyd 并改变其距离形状；它不能把该分量定义成零。远场线性 detrending 会进一步吸收真实 `1/D` 信号，因此 line-corrected branch 只用于 QC/sensitivity，不作为 hyd 定量主分支。

### 16.4 扩展实验的速度选择和拟合口径

速度应按 `eta(c) U_gap` 设计，而不是让所有浓度使用相同 nominal velocity。若用 20 nm 处的理论 hyd 幅值作为设计变量，则：

```text
U_target(c) = F_hyd,target(20 nm) / [19.4846 eta(c)].
```

下表给出三个可识别 hyd 水平的设计速度；它们是 acquisition-time/pilot 之前的物理基准，不是未经仪器验证就直接冻结的设置。

| wt% | `F20=50 pN` speed | `F20=100 pN` speed | `F20=200 pN` speed |
|---:|---:|---:|---:|
| 0 | 2.91 µm/s | 5.83 µm/s | 11.7 µm/s |
| 20 | 1.71 µm/s | 3.42 µm/s | 6.84 µm/s |
| 40 | 0.833 µm/s | 1.67 µm/s | 3.33 µm/s |
| 60 | 0.297 µm/s | 0.594 µm/s | 1.19 µm/s |
| 70 | 0.146 µm/s | 0.292 µm/s | 0.585 µm/s |
| 80 | 0.0585 µm/s | 0.117 µm/s | 0.234 µm/s |
| 90 | 0.0171 µm/s | 0.0342 µm/s | 0.0684 µm/s |
| 95 | 0.00776 µm/s | 0.0155 µm/s | 0.0310 µm/s |
| 99.5 | 0.00331 µm/s | 0.00662 µm/s | 0.0132 µm/s |

低浓度端的最高速度可能受 JPK bandwidth、采样率、snap/contact lag 和 surface-history 限制；高浓度端的最低速度可能受 map duration 和 drift 限制。因此正式设置允许分浓度带使用不同的 1:2:4 triplet，但必须在 pilot 前写明选择规则，并覆盖足够的 `eta U` leverage。不能看过 force 结果后再为每个浓度选择有利速度。

扩展数据的主观察模型记录为：

```text
F_app(D,U,c) = F_eq(D,c)
             + A(D,c) U_gap
             + b_map + m_map t
             + error,

A(D,c) approximately proportional to eta(c)/D.
```

主报告中的 20–200 nm force 必须是 map/block/session 层级支持的 `U -> 0` intercept `F_eq(D,c)`。随后才使用 equal-silica sphere–plane PB 模型拟合 effective screening length 和共同 boundary-potential magnitude。任何单一有限速度曲线都只能称为 apparent dynamic force，不能直接称为 equilibrium surface force。

继续沿用以下实验约束：

- 每张 8×8 map 内速度固定，以完整 map 作为 velocity experimental unit；pixel 是空间配对样本，不是独立速度重复。
- 每个浓度采用 time-balanced palindrome map order，并计算实际 `U_gap = -dD/dt`；nominal scanner speed 只作标签。
- 主拟合保留 raw/constant-referenced 数据，将 offset、time drift 和 hyd shape 联合估计；line-corrected 数据只作系统敏感性。
- approach/retract odd-even 可作辅助检验，但只有在 retract 的 20–200 nm 区间有经验证的 free non-contact support 时才能使用。
- 90–99.5 wt% 必须加强密闭、防吸水、前后 refractive-index/density 或实测 viscosity QC。25.6 °C 下 Cheng 模型给出的 99 与 99.5 wt% 黏度已相差约 10.8%，所以微小含水量误差会直接变成相同量级的 hyd force 误差。
- 高黏度下还要用 selected-concentration hold tests 检查 double-layer/charge-regulation relaxation；若 force 存在不可忽略的 dwell-time dependence，不能把全部 rate dependence 强制归入线性 hyd 项。

本扩展尚未改变现有 0、20、30、40 wt% 数据的证据等级；它记录的是后续 acquisition 和分析的设计边界。只有新的高浓度数据满足 map/session/order stability 和 `eta U/D` scaling 后，才能把零速外推称为已识别结果。

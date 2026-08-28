# 27-08-26 D4 纯水回文测试：完整分析

## 实验分组与缺失数据

- **Test A / map1–3**：三个轮换回文block构成完整测试；同浓度液体在map2与map3之间刷新。因此它用于量化换液/时间不连续性，不作为无干预速度因果的primary test。
- **Test B / map3–5**：map3、map4、map5构成无中途换液的primary test；map3同时是两套测试的共享bridge block。
- map5少保存最后一张0.1 µm/s map。分析没有补值：map5的0.1早/晚对标为missing；所有block都完整具有0.05与0.2早/晚对，因此primary速度估计使用block内0.05→0.2差分。map3/map4的0.1对用于检验线性。

回文顺序：

- map1: `0.2 → 0.1 → 0.05 → 0.05 → 0.1 → 0.2 µm/s`。
- map2: `0.1 → 0.05 → 0.2 → 0.2 → 0.05 → 0.1 µm/s`。
- map3: `0.05 → 0.2 → 0.1 → 0.1 → 0.2 → 0.05 µm/s`。
- map4: `0.2 → 0.1 → 0.05 → 0.05 → 0.1 → 0.2 µm/s`。
- map5: `0.1 → 0.05 → 0.2 → 0.2 → 0.05 µm/s`。

实际20–200 nm gap speed和单张map protocol时长为：

- nominal 0.05 µm/s：actual median 0.0491 µm/s，range 0.0476–0.0500 µm/s；protocol约 10.90 min/map。
- nominal 0.1 µm/s：actual median 0.0975 µm/s，range 0.0954–0.0996 µm/s；protocol约 5.58 min/map。
- nominal 0.2 µm/s：actual median 0.1944 µm/s，range 0.1907–0.1998 µm/s；protocol约 2.93 min/map。

## 力标定

温度采用用户给定 **25.6 °C**。D4空气热标定给出 `k=0.164829 N/m`；全部29张水中map的有效硬接触共同给出 global `InvOLS=51.822 nm/V`。文件内写入的sensitivity和force conversion未用于最终力。

| Cantilever | Air InvOLS (nm/V) | k (N/m) | f0 (kHz) | Q | unique TND / files |
|---|---:|---:|---:|---:|---:|
| D4 | 72.76 ± 3.76 | 0.1648 ± 0.0173 | 13.1394 | 49.61 | 4 / 5 |
| D5 | 79.90 ± 2.69 | 0.3113 ± 0.0217 | 16.5647 | 69.81 | 5 / 5 |
| D6 | 70.90 ± 1.06 | 0.1817 ± 0.0066 | 13.2811 | 52.48 | 5 / 5 |

D4有两个byte-identical TND文件；inventory保留两者，但SHA-256 dedup后只用4个unique spectra计算thermal mean和repeatability。该修正使D4 k相对重复计权结果改变约0.12%。
表中k的repeatability SD对D4为 10.5%；它是所有force共享的乘法scale uncertainty，未并入每个map/pixel的IQR。水中local map InvOLS范围为 50.972–54.687 nm/V；按既定规则仍对全部map使用同一个global值，local值只用于QC。

## Primary Test B：无换液block内速度结果

每个map block先对相同速度的早/晚同pixel曲线计算 `F_sym=(F_early+F_late)/2` 和 `F_history=(F_late−F_early)/2`。随后在每个block内用0.05与0.2 µm/s的实际gap speed计算slope和endpoint-linear U→0 intercept。下表的±SD来自map3/map4/map5三个block，不来自64个pixels。

| D (nm) | slopes map3 / map4 / map5 (pN/(µm/s)) | 3-block mean ± SD | hyd theory | ratio | U→0 force mean ± SD (pN) | mid-speed residual* (pN) |
|---:|---:|---:|---:|---:|---:|---:|
| 20 | 316.0 / -491.3 / -49.4 | -74.9 ± 404.3 | 17.16 | -4.4 | 1647.0 ± 814.5 | 51.5 |
| 50 | 112.0 / 67.4 / -40.0 | 46.5 ± 78.1 | 6.86 | 6.8 | 256.3 ± 247.0 | 15.0 |
| 100 | 6.1 / 9.1 / -4.8 | 3.5 ± 7.3 | 3.43 | 1.0 | 25.4 ± 28.6 | -0.2 |
| 200 | 5.5 / -1.3 / 1.2 | 1.8 ± 3.4 | 1.72 | 1.1 | 2.0 ± 0.6 | 0.1 |

注：mid-speed residual仅来自具有完整0.1 pair的map3/map4，定义为实际F_sym(0.1)减去0.05–0.2 endpoint-linear prediction。

如果是真正no-slip drainage，dF/dU应为正、跨block接近一致，并按1/D衰减。三个block的斜率在20 nm出现两种符号，在50–100 nm也有map5负斜率；因此上面的U→0 intercept是provisional calculation，不是validated equilibrium force。

### 时间模型和baseline敏感性

回文平均只严格抵消关于pair中心近似线性的时间项。本数据的force relaxation明显弯曲，因此又在17张原始map层级拟合 `F=平滑时间项+βU`；下表给出二至五次时间多项式以及block-fixed线性/二次时间模型得到的β范围。范围不是置信区间，而是model sensitivity。

| D (nm) | block-pair mean β: far-linear | block-pair mean β: far-constant | time-aware β range: far-linear | time-aware β range: far-constant | no-slip hyd |
|---:|---:|---:|---:|---:|---:|
| 20 | -74.9 | -61.1 | -298.3…-89.0 | -247.5…-55.9 | 17.16 |
| 50 | 46.5 | 40.9 | -6.3…98.4 | -6.9…87.4 | 6.86 |
| 100 | 3.5 | -10.1 | -2.6…16.5 | -13.4…2.0 | 3.43 |
| 200 | 1.8 | -8.8 | 0.9…3.6 | -8.8…-5.7 | 1.72 |

20 nm的time-aware速度项在两种baseline和全部时间模型下均为负，与drainage方向相反；50 nm随时间模型从略负到远高于hyd theory，跨过零。100–200 nm又对far-field零点定义发生量级甚至符号变化。所有target distance、两种baseline和全部time model的HC3 95% CI都包含零。因此当前数据没有任何距离区间满足可靠hyd subtraction和U→0 extrapolation所需的统计、model和baseline稳定性。

## 换液边界

以下为相同速度的F_sym(map3)−F_sym(map2)，包含换液、elapsed time和不同回文block order，不能称为纯换液因果效应：

| D (nm) | 0.05 µm/s | 0.1 µm/s | 0.2 µm/s |
|---:|---:|---:|---:|
| 20 | -48.8 pN | +238.6 pN | -31.5 pN |
| 50 | -38.9 pN | +39.6 pN | -35.8 pN |
| 100 | -1.6 pN | +0.4 pN | -5.9 pN |
| 200 | -0.1 pN | -1.8 pN | -3.3 pN |

## 无换液阶段仍存在的时间漂移

在map3–5的17张连续map中，map-median InvOLS与order的Spearman `ρ=-0.559`；far slope `ρ=-0.436`；pull-off `ρ=0.360`。这些是顺序采集的描述性趋势，p值不等于随机化速度因果检验。

从Test B第一张map3到最后一张map5，far-linear map median变化为：

| D (nm) | first (pN) | last (pN) | relative change |
|---:|---:|---:|---:|
| 20 | 2578.5 | 495.9 | -80.8% |
| 50 | 647.6 | 16.2 | -97.5% |
| 100 | 77.1 | 8.2 | -89.3% |
| 200 | 1.0 | 1.0 | -2.5% |

同期local contact InvOLS只从 51.508 变到 51.131 nm/V（-0.73%），terminal load也维持约21.8 nN；因此20–100 nm的巨大force relaxation不能由sensitivity gain或加载力变化解释。50 nm处完整的8个early/late pair全部为 `F_history<0`，即later map的force更小；在假定pair-sign独立时exact two-sided sign test `p=0.0078`。这些pair共享同一条顺序relaxation，因此p值只作方向性描述；8/8同号本身已经说明无换液并不等于stationary。

## Approach速度还改变了contact/retract状态

这里比较的是每个block内0.2 pair-symmetric值减去0.05 pair-symmetric值；map3/map4/map5分别为：

- InvOLS：+0.217 / +0.167 / +0.064 nm/V；不足以解释force差。
- terminal load：-0.018 / -0.034 / -0.011 nN；相对约21.8 nN很小。
- apparent contact height：+22.2 / -6.3 / +0.1 nm；符号不一致，不支持一个固定的speed-dependent contact-zero shift。
- approach snap-detected fraction：-0.289 / -0.133 / -0.008。高速度在map3/map4显著降低检测率，而且被检测事件的apparent distance可到约90–110 nm；这更像branch-shape/threshold QC异常，不能直接当真实snap-in位置。
- retract pull-off：+0.382 / +0.175 / +0.073 nN，三个block都变得更不负。retract protocol在所有map中相同，所以approach condition与后续接触/脱离状态明确相关；因速度顺序不是随机化的，尚不能区分真实approach-history因果与residual time confounding，但两者都不是可直接相减的approach hydrodynamic force。

## Far-field零点敏感性

- D=20 nm：linear-corrected − constant-referenced的29-map median为 -4.11 pN，范围 -41.48到+12.40 pN。
- D=50 nm：linear-corrected − constant-referenced的29-map median为 -3.68 pN，范围 -28.45到+4.39 pN。
- D=100 nm：linear-corrected − constant-referenced的29-map median为 -3.57 pN，范围 -20.37到+4.08 pN。
- D=200 nm：linear-corrected − constant-referenced的29-map median为 -1.44 pN，范围 -10.72到+2.43 pN。

## 结论边界

- map是速度/时间实验单位；8×8 pixels是同位置paired observations，未被当作独立速度重复。
- map5缺失0.1 late map不会影响三block的0.05→0.2 endpoint contrast；但它破坏了完整三速度回文的linearity/curvature诊断，且没有补值。lone early 0.1 map只进入chronology/time-aware model。
- Test A在同浓度refresh边界观察到force/contact-state不连续，但它与elapsed time和block order共变，不能量化纯refresh因果；map1–3因此不用于无干预速度因果的primary结论。
- Test B虽是目前最强的速度证据，但force随时间强烈、非线性衰减；block斜率不一致，100–200 nm又对baseline发生符号变化。当前结论是 **hyd slope不可唯一识别，不能据此做hyd subtraction或U→0 equilibrium recovery**。
- 20 nm的time-aware斜率方向与hyd相反，50 nm则强烈依赖时间模型；同时fixed-retract pull-off随approach speed改变。这支持额外的approach-history/contact-state systematic，而不是把速度差解释成单一drainage项。
- 当前输出不做PB/zeta/Debye拟合；先判定哪些距离的U→0 force可辨识，之后才应进入表面电势拟合。

## 输出

- `map_inventory_QC.csv`, `pixel_QC.csv`, `map_force_curves.csv`: 29-map raw reconstruction and QC.
- `pair_inventory.csv`, `palindrome_pair_curves.csv`: explicit missingness and early/late decomposition.
- `block_velocity_contrasts.csv`, `test_velocity_summary.csv`: block-level slopes and provisional U→0 intercepts.
- `primary_time_aware_speed_models.csv`, `primary_time_aware_model_predictions.csv`: map-level time/speed sensitivity models and fitted chronology.
- `pair_contact_state_QC.csv`, `pair_contact_state_velocity_contrasts.csv`: paired InvOLS/load/contact/snap/retract diagnostics.
- `liquid_refresh_discontinuity.csv`, `chronological_trends.csv`, `baseline_sensitivity_slices.csv`: intervention/time/baseline systematics.
- `fixed_pixel_row3_col3_slices.csv`, `pixel_force_slices_20_50_100_200nm.csv`: same-position examples and spatial target slices.
- `figures/`: chronology, block pair curves, velocity slopes, zero-speed intercepts, refresh jump, midpoint residual, fixed pixel, and 8×8 far-slope maps.
- `provenance.json`, `artifact_manifest.sha256`: complete raw hashes, parameters, software, and artifact identities.

# 27-08-26：不做hyd/time修正的apparent PB参数粗拟合

## 直接结论

这里按要求暂时忽略hydrodynamic force、acquisition-history relaxation和contact-state变化，把每条finite-speed F–D曲线直接当成equilibrium curve拟合。结果中的Debye length和surface potential因此是 **apparent/model-conditioned parameters**，不是已经测得的bulk Debye length或平衡surface potential。

**直接回答：存在非常强的apparent时间依赖；当前没有可辨识的速度依赖。**

在无换液Test B中，line-corrected primary fit从acquisition #13到#29（128.6 min），apparent λD由 21.62 降到 7.36 nm（-66.0%），|ψ|由 76.0 降到 41.5 mV（-45.4%）。这个趋势在两种far-field baseline和10–250/20–200 nm两个window下均保持。

相反，Test B的五种time-adjusted map模型中没有一个速度系数的HC3 95% CI排除零；回文Fsym的0.2−0.05 µm/s参数差在map3/map4/map5之间变号，n=3 paired t-test也不显著。因此速度对apparent PB参数的独立影响至多是当前误差下的几个百分点。

- `linear_drift_corrected`：acquisition #1→#29的apparent λD为 16.78→7.36 nm，|ψ|为 71.5→41.5 mV。
- `far_constant_referenced`：acquisition #1→#29的apparent λD为 16.78→7.36 nm，|ψ|为 70.0→41.5 mV。

下面把时间趋势、time-adjusted speed coefficient以及回文Fsym的0.2−0.05 µm/s block contrast分开报告。若速度结果随baseline或time model变号，就不能称为稳定速度依赖。

## 采用的旧公式

同材料silica sphere–plane、equal constant-potential nonlinear PB Derjaguin模型：

`F(D)=2πR ε (kBT/e)^2 κ G(κD,e|ψ|/kBT) − A_H R/(6D²) + b`。

固定 `R=4.546849 µm`, `A_H=2.4e-21 J`, `εr=78.5`, `T=25.6 °C`；每条曲线只自由拟合 `λD`, `|ψ|` 和常数baseline `b`。silica按负号约定输出 `ψ=-|ψ|`，但force本身不能从该同表面模型判定电势符号。primary window沿用之前的10–250 nm，同时保留20–200 nm sensitivity branch。

## 时间和速度检验

| baseline | parameter | all-29 time Spearman ρ / p | Test-B time ρ / p | time-model 0.2 vs 0.05 effect range | HC3 CI excludes 0 | palindrome blocks 3/4/5 high-vs-low | block t p |
|:---|:---|:---|:---|:---|---:|:---|---:|
| linear_drift_corrected | Apparent Debye length | -0.484 / 0.00781 | -0.999 / 2.95e-23 | +1.6%…+2.2% | 0/5 | +0.0%/+6.5%/-3.1% | 0.742 |
| linear_drift_corrected | Apparent surface-potential magnitude | -0.505 / 0.00518 | -0.987 / 2.37e-13 | -1.9%…-1.2% | 0/5 | +0.0%/-6.9%/+3.4% | 0.720 |
| far_constant_referenced | Apparent Debye length | -0.507 / 0.00504 | -1.000 / 0 | +3.1%…+3.4% | 0/5 | -3.1%/+6.5%/+0.0% | 0.742 |
| far_constant_referenced | Apparent surface-potential magnitude | -0.484 / 0.00787 | -0.989 / 7.58e-14 | -3.5%…-3.0% | 0/5 | +4.4%/-7.0%/-3.5% | 0.587 |

- time-model range来自Test B的17张map，分别使用global一至三次时间多项式以及block-fixed线性/二次within-block时间模型；它仍是顺序数据的敏感性分析。
- palindrome contrast先拟合每个block/speed的Fsym曲线，再比较0.2与0.05 µm/s；三个百分数依次对应map3/map4/map5，实验重复数只有n=3。

## 拟合质量与边界

| dataset | baseline | window | fits | median R² | minimum R² | grid-boundary fits |
|:---|:---|:---|---:|---:|---:|---:|
| map | linear_drift_corrected | previous_primary_10_250nm | 29 | 0.9929 | 0.9774 | 0 |
| map | linear_drift_corrected | measurement_window_20_200nm | 29 | 0.9993 | 0.9961 | 0 |
| map | far_constant_referenced | previous_primary_10_250nm | 29 | 0.9902 | 0.9524 | 0 |
| map | far_constant_referenced | measurement_window_20_200nm | 29 | 0.9983 | 0.9771 | 0 |
| palindrome_Fsym | linear_drift_corrected | previous_primary_10_250nm | 14 | 0.9936 | 0.9614 | 0 |
| palindrome_Fsym | linear_drift_corrected | measurement_window_20_200nm | 14 | 0.9983 | 0.9858 | 0 |
| palindrome_Fsym | far_constant_referenced | previous_primary_10_250nm | 14 | 0.9819 | 0.9197 | 0 |
| palindrome_Fsym | far_constant_referenced | measurement_window_20_200nm | 14 | 0.9959 | 0.9676 | 0 |

每个distance bin的权重尺度为 `max(pixel IQR/1.349, 2 pN)`，是空间离散度而不是mean的standard error；相邻D bins也高度相关。因此CSV中的profile ranges只能说明grid局部可辨识性，不能当正式95% CI。

## 解释边界

- 这个计算回答的是：如果强行把每条曲线解释成同一个PB+vdW equilibrium模型，拟合参数怎样随测量变化。它不回答真实λD或ψ是否随时间改变。
- 同一真实表面在约数小时内出现大幅apparent参数变化，首先说明未建模force relaxation被PB参数吸收；不能优先解释成纯水离子强度或silica化学真的同步改变。
- 速度依赖只有在不同block、baseline、fit window和time model下同号且量级稳定时才可信；当前表格直接显示这一稳定性。
- `|ψ|`与`λD`在有限窗口内相关，且同材料force不识别ψ符号；负号仅来自silica convention。

## 数值检查

- 共享PB网格与原 `total_equilibrium_force_pN` 的最大绝对差为 `5.457e-12 pN`，最大相对差为 `6.204e-15`。
- 网格为 181个λ点（1–300 nm）× 167个|ψ|点（1–250 mV）；所有模型值、fit scales、RSS和输出参数均检查finite。
- nonlinear PB evaluator沿用原脚本已经验证的small-potential linear limit和far-field asymptote；本脚本另做公式逐点identity检查。

## 输出

- `apparent_pb_fits.csv`: 29张map与14条palindrome Fsym曲线的两种baseline、两个window拟合。
- `apparent_pb_fit_predictions.csv`: primary 10–250 nm曲线、EDL/vdW/baseline分量和residual。
- `parameter_time_tests.csv`, `parameter_time_adjusted_speed_models.csv`: 时间相关与time-adjusted速度系数。
- `palindrome_parameter_contrasts.csv`, `palindrome_parameter_tests.csv`: n=3 block-level 0.2−0.05 µm/s比较。
- `parameter_summaries.csv`, `figures/`, `provenance.json`, `artifact_manifest.sha256`: 汇总、图、参数和身份记录。

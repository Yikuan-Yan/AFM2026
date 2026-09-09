# 31-08-26 map2–5 advanced same-pixel analysis

## 结论先行

现有数据确认20–50 nm曲线不只是整体幅值下降，归一化后的曲线形状也发生变化；但contact/separation零点、空间斑块、baseline和模型参数耦合仍足以显著改变拟合参数。因此这里报告的是 **apparent/model-conditioned** `lambda_D` 与 `|psi|`，不是已经验证的bulk Debye length或zeta potential。

Primary sphere–plane nonlinear-PB fit给出的map2→map5中位apparent `lambda_D`为 `23.78 → 11.29 nm`，`|psi|`为 `26.5 → 16.0 mV`。

Model-free endpoint-normalized force在50 nm的map5−map2 paired median为 `-0.225`，row-signflip `p=0.01562`，row-cluster 95% CI `[-0.296, -0.156]`；30–80 nm local log-slope得到的decay length变化为 `-6.87 nm`，但其row-signflip `p=0.2578`且cluster CI跨零，因此该window本身不构成稳健证据。

把所有曲线统一平移 `D0=-10…+10 nm` 后，map5−map2中位lambda变化仍处在 `-13.60…-12.27 nm`；对应potential变化 `-12.0…-5.2 mV`。但允许每条曲线从七档D0中自由选择时，`73.8%` 落在搜索边界，说明D0与PB参数不能由当前20–200 nm窗口稳定共同识别。

## Primary sphere–plane nonlinear-PB结果

模型为equal constant-potential silica sphere–plane nonlinear 1:1 PB Derjaguin + sphere–plane vdW + fixed ideal no-slip lubrication + per-curve constant baseline。固定 `R=4.546849 µm`, `A_H=2.4e-21 J`, `epsilon_r=78.5`, `T=25.6 °C`。

| map | time (min) | apparent lambda_D (nm), median [IQR] | apparent |psi| (mV), median [IQR] | RMSE median (pN) | R2 median |
|---:|---:|:---|:---|---:|---:|
| 2 | 0.00 | 23.78 [20.14, 28.76] | 26.5 [22.0, 34.0] | 10.56 | 0.987 |
| 3 | 27.16 | 17.32 [13.77, 21.62] | 25.0 [20.5, 32.5] | 14.34 | 0.967 |
| 4 | 53.33 | 13.44 [10.77, 18.60] | 22.0 [17.5, 28.0] | 13.70 | 0.934 |
| 5 | 78.66 | 11.29 [10.03, 13.02] | 16.0 [11.5, 22.0] | 7.95 | 0.925 |

### Same-pixel map2→map5统计

| metric | map2 median | map5 median | paired delta median [IQR] | row-signflip p | row-cluster 95% CI |
|:---|---:|---:|:---|---:|:---|
| apparent_lambda_D_nm | 23.8 | 11.3 | -12.7 [-18, -8.9] | 0.01562 | [-15.4, -9.11] |
| apparent_zeta_magnitude_mV | 26.5 | 16 | -10.5 [-13.5, -7.12] | 0.007812 | [-12, -9] |
| force_20nm_pN | 392 | 144 | -202 [-340, -107] | 0.007812 | [-259, -153] |
| force_50nm_pN | 136 | 26 | -97.2 [-177, -29.8] | 0.007812 | [-140, -73.1] |
| contact_height_nm | 6.1e+03 | 6.55e+03 | +427 [+407, +450] | 0.007812 | [+414, +440] |
| local_contact_InvOLS_nm_per_V | 71.5 | 67.9 | -1.24 [-3.3, +2.02] | 0.1406 | [-2.58, +0.364] |

## Contact-zero与参数可识别性

固定D0敏感性与per-curve离散D0 nuisance结果分别保存在 `contact_zero_fixed_offset_summary.csv` 与 `contact_zero_profiled_nuisance.csv`。固定共同D0不能消除map2→map5趋势；但per-curve D0 profile大量触边，不能把其最优D0直接解释为真实contact drift。

Map-median及逐曲线lambda–potential profile均使用far-noise或spatial-IQR scale；相邻distance bins相关，所以profile范围不是formal 95% CI。

尤其map5的map-median spatial-IQR profile在lambda上达到 `300 nm` grid上界、potential达到 `250 mV`上界；这说明用一条map median曲线做formal参数区间是不可靠的。逐曲线profiles通常更窄，但仍受相邻D bins相关影响。

## Joint/hierarchical parameter-sharing

在共享参数模型中，pseudo-BIC最低的是 `map_shared_lambda_and_zeta`；leave-one-physical-pixel-out shape RMSE最低的也是 `map_shared_lambda_and_zeta`（mean `63.92 pN`）。两项均要求map-specific lambda和potential；一个跨四图完全共享的参数对不足以描述数据。

每条曲线仍保留自己的constant baseline。CV每次留出同一physical pixel的四条曲线，用其150–200 nm只估baseline，再预测20–145 nm shape；因此没有把held-out lambda/zeta重新拟合回来。

## 空间结构与QC关联

共 `70` 个map/difference空间检验中，Moran permutation在BH校正后 `q<0.05` 的有 `58` 个。64个pixel因此不能普遍视为64个独立重复；报告以row-signflip和row-cluster bootstrap为主。

map2→map5的spatial correlation在20 nm force为 `+0.682`，到50 nm仅 `+0.062`；apparent lambda甚至为 `-0.475`，而potential仍为 `+0.808`。因此变化不是一个统一global scale factor：50 nm和lambda的空间pattern发生了重排。

对force/fit-parameter变化与InvOLS、far slope、noise、load、speed、pull-off、contact height变化的 `48` 个Spearman诊断中，BH `q<0.05` 的有 `11` 个；这些仍是association，不是因果校正。

其中50 nm force变化与far-field slope变化的关联最强：Spearman `rho=+0.873`，BH `q=2.47e-19`。这说明50 nm绝对force下降中有明显baseline/drift耦合，不能全部解释为EDL变化。

相反，endpoint-normalized 50 nm shape变化与far-field slope变化没有显示同样的单调关联：`rho=-0.119`，BH `q=0.568`。因此已观察到的normalized shape separation不能仅由“far-slope越变、50 nm force越变”这一单一关系解释，但这并不排除其他baseline curvature或history systematics。

Raw fitted contact-height scanner coordinate从map2到map5整体增加约427 nm，并保持很强空间相关；该坐标在每条曲线构造separation时已被逐条减去，因此它反映scanner/sample coordinate creep或offset，而不是未校正的427 nm物理gap变化，也不能直接等同于PB fit中的D0。

## Model comparison与residual

- `sphere_plane_nonlinear_PB_fixed_no_slip_hydrodynamics`: four-map median RMSE range `7.95–14.34 pN`; per-curve AIC winner fraction range `20.3%–65.6%`.
- `sphere_plane_nonlinear_PB_no_hydrodynamics`: four-map median RMSE range `7.78–14.49 pN`; per-curve AIC winner fraction range `12.5%–56.2%`.
- `sphere_plane_linear_HHF_fixed_no_slip_hydrodynamics`: four-map median RMSE range `7.91–14.41 pN`; per-curve AIC winner fraction range `4.7%–12.5%`.
- `equal_sphere_linear_HHF_fixed_no_slip_hydrodynamics`: four-map median RMSE range `7.78–14.78 pN`; per-curve AIC winner fraction range `9.4%–12.5%`.

Primary residual的lag-1 correlation中位数为 `+0.583`；naive runs test `p<0.05` 的曲线占 `71.5%`。这表明高R2不能代替residual结构检查。

## Calibration与time-law敏感性

D5 `k` repeatability SD对应common force scale约 `±7.0%`。在这两个端点下，四图中位lambda相对nominal的最大偏移为 `0.77 nm`；potential最大偏移 `1.5 mV`。InvOLS还会改变distance axis，不能由此common-scale branch完全代表。

四时点对linear、exponential-to-zero及exponential-with-asymptote共生成 `15` 个诊断拟合。四点不足以识别kinetic law；tau或asymptote只能作为探索性描述，不能用来证明液体浓度随时间变化。

在median gap speed下，固定ideal no-slip sphere–plane hydrodynamic force为20 nm `34.25 pN`、50 nm `13.70 pN`。它分别约占map2中位measured force的 `8.7%` / `10.1%`，到map5则约为 `23.7%` / `52.6%`。这些只是固定no-slip理论项相对line-corrected measured force的量级，不是数据已独立识别出的hyd比例。

## 数值与物理边界

- PB library formula identity max error `9.095e-13 pN`; prediction closure `0.000e+00 pN`; core checks PASS.
- Primary nonlinear PB避免了linear-HHF的small-potential限制，但仍假定symmetric 1:1 electrolyte、equal constant potential、Derjaguin、固定Hamaker和ideal no-slip。纯水中实际离子种类/CO2并未由force curve确定。
- Hydrodynamic term是理论固定项；当前same-speed四图不能从数据本身独立验证其幅度或slip boundary condition。
- AIC/BIC、profile和distance-wise p值均受同一曲线内distance correlation影响，已明确标为diagnostic。
- 接触零点、baseline、共同force scale、空间相关和模型选择分别做了sensitivity branch；没有把任一branch提升为bulk Debye length或zeta potential测量。

## 主要输出

- `model_free_*`: endpoint-normalized shape、local log-slope及paired tests。
- `contact_zero_*`: 七档fixed D0及per-curve discrete nuisance profile。
- `surface_model_comparison_*`, `sphere_plane_nonlinear_predictions.csv`: 四模型统一OLS比较和primary逐点预测。
- `same_pixel_*`, `spatial_*`, `global_scaling_vs_spatial_change.csv`: paired trajectories、空间统计和global-vs-patch变化。
- `per_curve_parameter_profile_intervals.csv`, `map_median_parameter_profiles.csv`: identifiability诊断。
- `hierarchical_*`: parameter-sharing、pseudo-information criteria和leave-one-pixel-out shape prediction。
- `residual_*`, `four_timepoint_time_law_diagnostics.csv`, `calibration_common_force_scale_sensitivity.csv`。
- `NUMERICS_AUDIT.md`, `provenance.json`, `artifact_manifest.sha256`, `figures/`。

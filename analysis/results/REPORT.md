# 三支 cantilever 标定报告

实验温度：**25.6 °C（298.75 K）**。

## 建议写入 JPK 的标定值

| Cantilever | Sensitivity / InvOLS (nm/V) | Spring constant (N/m) | f₀ (kHz) | Q |
|---:|---:|---:|---:|---:|
| 1 | 85.30 ± 0.96 | 0.2970 ± 0.0093 | 16.1101 ± 0.0016 | 68.91 ± 1.37 |
| 2 | 65.72 ± 0.62 | 0.2736 ± 0.0063 | 15.0774 ± 0.0007 | 64.75 ± 0.87 |
| 3 | 91.56 ± 0.43 | 0.2384 ± 0.0034 | 15.0675 ± 0.0011 | 67.62 ± 0.39 |

表中的 ± 是 5 次重复测量的 **repeatability sample SD**（cantilever 2 的错误导出峰也已从原始谱重拟合，因此仍为 5 次），不是可溯源的绝对标定不确定度。实际输入时可分别使用：

- Cantilever 1：`sensitivity = 85.30 nm/V`，`spring constant = 0.2970 N/m`。
- Cantilever 2：`sensitivity = 65.72 nm/V`，`spring constant = 0.2736 N/m`。
- Cantilever 3：`sensitivity = 91.56 nm/V`，`spring constant = 0.2384 N/m`。

## 谱选错峰的处理

`calibration/2/thermal-noise-data_vDeflection_2026.08.18-14.47.49.tnd` 的导出 header 写成 **97.67 kHz**。本分析没有沿用这个值，也没有删除整条谱。全谱 blind peak search 在原始 measured PSD 中先找到 15.039 kHz 与 97.668 kHz 两个显著峰；按“最低显著共振峰 = fundamental”选取前者，再做 SHO 拟合，得到 **f₀ = 15.0777 kHz, Q = 64.52**。97.67 kHz 峰保留为高阶模态证据，但不进入基频热标定。

峰选择不使用照片/厂家给出的 11–17 kHz 范围；该范围只在计算完成后用于 sanity check。三组重拟合的 f₀ 均落在照片标注范围内。

## 计算方法

1. 每个 `.tnd` 直接读取 `Frequency` 与 measured `average` PSD；导出的 `fit-data`、`parameter.f`、`parameter.Q`、`parameter.A` 不参与计算。对 1–200 kHz 的 log-PSD 去除宽尺度 median baseline，以 ≥0.8 decade prominence 搜峰，并取最低显著峰作为 fundamental。
2. 在选中峰的 ±3 kHz 内，以 Gamma/Whittle PSD likelihood 拟合

   `S_VV(f) = N + A² / [(1-(f/f₀)²)² + (f/(f₀Q))²]`。

   一侧共振面积为 `I_V = (π/2) A² Q f₀`，单位 V²。
3. 使用文件记录的 rectangular/dynamic correction factor `β = 0.8170` 与用户给定温度计算 voltage-domain thermal factor：`B_z = β k_B T / I_V`（单位 N·m/V²）。
4. `.jpk-force` 使用每条 approach segment 的 measuredHeight 和 vDeflection 原始 int32 数据及各自 metadata conversion 解码；在末端 50 nm 硬接触线性区拟合 `V = a + b z`，取 `Sensitivity = 1/|b|`。
5. 最终 vertical spring constant 为 `k_z = mean(B_z) / mean(Sensitivity)²`。这也正是 JPK 文件中 corrected vertical thermal factor 对 contact sensitivity 的换算关系。

## QC 与数值敏感性

- Cantilever 1：5 条 contact fit 的最低 R² = 0.999616；contact span 从 50 nm 改为 40/60 nm 时，k 最大移动 0.0024 N/m；thermal fit half-width 从 3 kHz 改为 2/4 kHz 时，k 最大移动 0.0017 N/m。
- Cantilever 2：5 条 contact fit 的最低 R² = 0.999659；contact span 从 50 nm 改为 40/60 nm 时，k 最大移动 0.0010 N/m；thermal fit half-width 从 3 kHz 改为 2/4 kHz 时，k 最大移动 0.0004 N/m。
- Cantilever 3：5 条 contact fit 的最低 R² = 0.999805；contact span 从 50 nm 改为 40/60 nm 时，k 最大移动 0.0016 N/m；thermal fit half-width 从 3 kHz 改为 2/4 kHz 时，k 最大移动 0.0005 N/m。
- 15 个 JPK force ZIP 均通过 CRC；15 条 thermal spectrum 均为 31,131 个递增 frequency bins，无 NaN/Inf，15 次独立 SHO optimizer 均成功。
- 照片标签给出的 CONT-W batch 范围为 f₀ = 11–17 kHz、C = 0.11–0.56 N/m；三个最终结果均在此范围。该比较仅作 ex-post physical sanity check。

## 解释边界

- Sensitivity 拟合假定 force curve 的接触基底相对 cantilever 足够刚。仅凭文件不能确认基底材料；若基底可压缩，InvOLS 会偏大，继而 k 会偏小。
- repeatability SD 不含 optical-lever spot position、硬接触几何、基底 compliance、0.817 correction model、scanner calibration 等 systematic uncertainty，因此不能当作 traceable absolute uncertainty。
- 这里报告的是 JPK vertical convention 的 `k_z`；不要把 thermal-factor header 中错误显示的 `N/m` 当成该中间量的单位，中间量实际为 N·m/V²。

## 输出

- `thermal_refits.csv`：每条原始谱的 blind candidates、独立 f₀/Q/A/area 与 thermal factor。
- `force_contact_fits.csv`：每条 force curve 的 contact slope、Sensitivity、R² 与 40/50/60 nm window sensitivity。
- `calibration_summary.csv`：最终三支 cantilever 标定值及 repeatability/window sensitivity。
- `figures/thermal_full_spectra_refits.png`、`thermal_fundamental_zoom.png`、`force_contact_fits.png`：原始 measured data 与重拟合。
- `provenance.json` 与 `artifact_manifest.sha256`：输入身份、参数、软件版本和输出 hashes。

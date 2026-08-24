# 0 wt%三张map的分布分离置信度

## 结论边界

三张map的采集顺序为 `2 um/s -> 1 um/s -> 4 um/s`，每个速度只有一张map。主分析使用三图共同有效的 `255` 个物理点位和 `16` 个4x4空间block。

因此，下述95%区间是已测10x10 um区域内、以空间block为重采样单位的条件置信度；它不能解释为速度处理在重复实验中的95%置信度。速度与map时间/顺序完全混杂，真正的速度实验单位每组n=1。

## 方法

在25-250 nm的每个5 nm距离bin，对相同物理pixel计算 `F_fast - F_slow`。主效应是配对中位差；分布级效应是同点位 `P(F_fast > F_slow)`。4x4 pixel block bootstrap共 10000 次。simultaneous band使用studentized max-T，一次覆盖三组比较和全部46个距离bin。

## 全窗口检验

| comparison | exact block max-T p | Holm p | simultaneous 95% band excludes 0 | 8x8-block sensitivity p |
|:---|---:|---:|:---|---:|
| 2-1 um/s | 3.0518e-05 | 9.1553e-05 | 30-105 nm (positive); 115-115 nm (positive) | 0.125 |
| 4-2 um/s | 3.0518e-05 | 9.1553e-05 | 25-250 nm (negative) | 0.125 |
| 4-1 um/s | 3.0518e-05 | 9.1553e-05 | 25-145 nm (negative); 180-190 nm (negative); 200-250 nm (negative) | 0.125 |

4x4结果量化区域内分离的一致性。8x8敏感性只有4个空间block，两侧exact sign-flip检验的最小可达p值受离散分辨率限制；它用于显示结论对更长空间相关尺度的依赖。

## 选定距离

数值为fast-minus-slow配对中位差；括号依次为pointwise 95% CI和全家族simultaneous 95% CI。

### 2-1 um/s

| D (nm) | median difference (pN) | pointwise 95% CI (pN) | simultaneous 95% CI (pN) | P(fast > slow) |
|---:|---:|:---|:---|---:|
| 25 | 22 | [3.92, 42.5] | [-14, 57.9] | 0.596 |
| 50 | 53 | [40.1, 62.9] | [27.5, 78.4] | 0.831 |
| 100 | 16.9 | [11.6, 21.3] | [8.04, 25.8] | 0.659 |
| 150 | 1.62 | [-3.37, 6.9] | [-9.13, 12.4] | 0.514 |
| 200 | -1.45 | [-5.02, 3.22] | [-9.76, 6.87] | 0.478 |
| 250 | -0.841 | [-5.15, 4.97] | [-9.77, 8.09] | 0.490 |

### 4-2 um/s

| D (nm) | median difference (pN) | pointwise 95% CI (pN) | simultaneous 95% CI (pN) | P(fast > slow) |
|---:|---:|:---|:---|---:|
| 25 | -94.1 | [-115, -79.5] | [-133, -55.5] | 0.192 |
| 50 | -129 | [-139, -118] | [-149, -108] | 0.043 |
| 100 | -44.7 | [-53.1, -36.8] | [-61.7, -27.7] | 0.161 |
| 150 | -19.6 | [-26.4, -14.4] | [-33.2, -6.05] | 0.322 |
| 200 | -18.6 | [-25.3, -12.5] | [-31.4, -5.82] | 0.318 |
| 250 | -15.3 | [-20.9, -7.47] | [-29.2, -1.31] | 0.329 |

### 4-1 um/s

| D (nm) | median difference (pN) | pointwise 95% CI (pN) | simultaneous 95% CI (pN) | P(fast > slow) |
|---:|---:|:---|:---|---:|
| 25 | -90.9 | [-106, -66.1] | [-132, -49.6] | 0.200 |
| 50 | -77.9 | [-85.6, -68.1] | [-95.9, -59.9] | 0.102 |
| 100 | -23.2 | [-26.9, -17.4] | [-32.7, -13.7] | 0.286 |
| 150 | -11.5 | [-18.4, -5.7] | [-24.4, 1.41] | 0.380 |
| 200 | -12.1 | [-17, -9.42] | [-19.3, -4.86] | 0.329 |
| 250 | -8.65 | [-10.8, -6.11] | [-13.2, -4.12] | 0.369 |

## 独立force curves为何不进入主置信度

0 wt%独立force curves在1/2/4 um/s分别只有3/2/0条，既不平衡也不能与map pixel同点位配对。它们继续参与raw hard-contact sensitivity共识，但若加入速度分布检验会把不同层级的观测混为独立重复。

## 文件

- `zero_wt_distribution_separation.csv`: 全部距离的效应、pointwise和simultaneous区间。
- `zero_wt_distribution_separation_selected.csv`: 六个选定距离。
- `zero_wt_distribution_separation_global.csv`: functional max-T检验和block-size敏感性。
- `figures/zero_wt_distribution_separation_confidence.png`: 置信带图。

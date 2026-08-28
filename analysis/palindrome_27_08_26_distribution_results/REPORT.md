# 27-08-26 回文平均后的速度分布与paired t-test

## 直接结论

回文平均后，若把64个pixel当样本，部分单独block和距离的速度分布确实给出很小的paired t-test p-value；但空间tile检验显著减少，而且map3/map4/map5的速度差会变号。以block作为真正实验重复时，Test B在20、50、100、200 nm均不显著；在完整20–200 nm扫描中最小raw p也只有 `0.0287`（80 nm），Holm后为 `1.0000`。全部592个可计算的n=3 block-level sensitivity tests中，没有一个Holm-adjusted p<0.05。因此当前证据是 **conditional within-map separation存在，但replicated velocity separation未建立**。

## 问题与统计单位

对每个block、速度和同一physical pixel先计算 `F_sym=(F_early+F_late)/2`。如果时间漂移在early/late之间近似线性，且各速度pair具有共同中心时刻，该平均会消除一阶时间项。t-test使用同位置差值 `F_sym(v_high)-F_sym(v_low)`，因此是paired而不是Welch independent-sample test。

结果分三层：

1. 64-pixel paired t-test：回答本次已测区域内同位置值是否系统变化，但pixels空间相关，会低估standard error。
2. 16个2×2 spatial-tile paired t-test：先对局部tile内差值取均值，作为较保守的区域内分离检验。
3. block-level paired t-test：每个block/速度先取64-pixel spatial median，再以三个block作为n=3实验重复。这一层才对应跨map的速度推断，但自由度只有2。

paired t-test检验的是mean paired difference是否为零，不等价于检验两个完整distribution的所有形状都不同。

## 回文在clock time上的对称程度

| block | complete speeds | pair-center spread (min) | largest early/late half-span (min) |
|---:|---:|---:|---:|
| map1 | 0.05, 0.1, 0.2 | 1.299 | 23.59 |
| map2 | 0.05, 0.1, 0.2 | 0.424 | 18.75 |
| map3 | 0.05, 0.1, 0.2 | 3.998 | 20.49 |
| map4 | 0.05, 0.1, 0.2 | 0.293 | 19.40 |
| map5 | 0.05, 0.2 | 0.179 | 9.54 |

map3的三种速度pair center相差约4 min，因此原始回文平均并不严格对应同一clock time。`linear_center_aligned` sensitivity branch用所有同block early/late差值估计一个共同线性drift，再把每个F_sym平移到该block的平均pair center；它只能检验一阶中心错位，不能消除已观察到的非线性relaxation。

## Primary Test B：0.2与0.05 µm/s

以下采用far-linear baseline和未经额外对齐的palindrome mean。ΔF为0.2−0.05 µm/s。pixel/tile的Holm p在每个block和contrast内跨20–200 nm的37个距离校正；最后一列以map3/map4/map5三个block spatial median作paired t-test。

| D | map3 Δmean / pixel Holm / tile Holm | map4 | map5 | n=3 block-median mean Δ / raw p / Holm p |
|---:|---:|---:|---:|---:|
| 20 nm | +57.70 pN / 0.2004 / 1.0000 | -75.25 pN / 1.5e-05 / 0.0415 | +7.37 pN / 1.0000 / 1.0000 | -13.92 pN / 0.7572 / 1.0000 |
| 50 nm | +19.71 pN / 0.0016 / 0.1129 | +9.93 pN / 0.1777 / 1.0000 | -5.51 pN / 0.1074 / 0.2619 | +8.58 pN / 0.3379 / 1.0000 |
| 100 nm | +0.42 pN / 1.0000 / 1.0000 | +3.65 pN / 1.0000 / 1.0000 | -0.42 pN / 1.0000 / 1.0000 | +0.92 pN / 0.4053 / 1.0000 |
| 200 nm | +0.24 pN / 1.0000 / 1.0000 | +0.41 pN / 1.0000 / 1.0000 | -0.52 pN / 1.0000 / 1.0000 | +0.16 pN / 0.4396 / 1.0000 |

逐block pixel/tile检验即使显著，也只说明该block、该区域的paired distributions具有非零mean shift。真正的block-level结果由三个Δ值组成；如果它们变号或离散很大，64 pixels不能增加独立速度重复数。n=3 exact sign test即使3/3同号，two-sided最小p也为0.25。

Test A的far-linear palindrome endpoint在20 nm给出raw block-level p=0.0103，但它跨越换液干预，而且在预先保留的全距离/全contrast Holm family中p=1.000；far-constant对应raw p=0.0045、Holm p=0.499。它最多是下一轮实验的候选距离，不能作为当前速度因果结论。

## 显著距离范围

下表只看Test B的0.2−0.05 endpoint、far-linear baseline。`+`表示high speed force更大，`−`相反；Holm correction跨37个距离。

| preprocessing | block | 64-pixel paired t | 16-tile paired t |
|---|---:|---|---|
| palindrome_mean | map3 | 25–60 nm + | none |
| palindrome_mean | map4 | 20–30 nm −; 55–80 nm + | 20–20 nm −; 70–80 nm + |
| palindrome_mean | map5 | 45–45 nm − | none |
| linear_center_aligned | map3 | 20–110 nm + | 25–100 nm + |
| linear_center_aligned | map4 | 20–30 nm −; 55–75 nm + | 20–20 nm −; 70–75 nm + |
| linear_center_aligned | map5 | 45–45 nm − | none |

## 线性中心时刻对齐和baseline敏感性

- 20 nm：map3–5的median absolute clock-center adjustment为 map3 14.97 pN, map4 2.43 pN, map5 2.53 pN。
- 50 nm：map3–5的median absolute clock-center adjustment为 map3 8.04 pN, map4 0.85 pN, map5 0.17 pN。
- 100 nm：map3–5的median absolute clock-center adjustment为 map3 1.63 pN, map4 0.02 pN, map5 0.06 pN。
- 200 nm：map3–5的median absolute clock-center adjustment为 map3 0.22 pN, map4 0.01 pN, map5 0.03 pN。

原始palindrome mean、linear-center alignment、far-linear和far-constant四种组合都保留在CSV。若显著区间或ΔF方向随这些合理处理改变，则不能称为robust velocity separation。特别是center alignment仍假设每个block存在一个共同线性drift；当前chronology显示relaxation具有曲率，而且该branch的t-test没有传播drift-estimation uncertainty，所以它是敏感性分析而不是修正后的ground truth。

## 结论

- 回文平均确实降低了同速度early/late的一阶时间偏差，并允许在同一block内做严格same-pixel paired comparison。
- pixel-level paired t-test回答的是已测空间区域内是否分离；它不能把64个pixel变成64次独立实验。2×2 tile结果更保守，但仍不是独立map replication。
- 速度因果的主统计单位是block。Test B的0.05–0.2 endpoint只有n=3；涉及0.1 µm/s的Test B block-level比较只有map3/map4两组，未执行正式t-test。
- 因pair center并非完全同时、relaxation又明显非线性，回文平均不能保证彻底消除时间效应。是否可主张velocity separation必须以n=3 block-level一致性和处理敏感性为准，而不能只看很小的pixel p-value。

## 输出

- `palindrome_pixel_paired_ttests.csv`: 64-pixel paired t-tests and multiplicity corrections.
- `palindrome_spatial_tile_paired_ttests.csv`: sixteen 2x2-tile paired t-tests.
- `palindrome_experimental_unit_ttests.csv`: n=3 block-level tests and n<3 exclusions.
- `pair_center_timing.csv`, `linear_center_alignment_QC.csv`: clock-time symmetry and alignment sensitivity.
- `palindrome_symmetric_pixel_forces.npz`: all F_sym pixel distributions used by the tests.
- `selected_target_summary.csv`, `figures/`: target-distance tables and visualizations.
- `provenance.json`, `artifact_manifest.sha256`: upstream identities, definitions, software and hashes.

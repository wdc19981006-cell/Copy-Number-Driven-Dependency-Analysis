# Copy Number–Driven Dependency Analysis

模块化 **R 4.5.0** 分析框架，用于检验 Gene A 拷贝数与 Gene B CRISPR dependency 的关联。当前正式案例为 [VPS4B → VPS4A](modules/VPS4B_VPS4A/README.md)。支持 Windows 和 32 GB RAM；数据库不会上传仓库。

用户提供的 [R 代码原文](docs/USER_SUPPLIED_ANALYSIS.R) 是统计核心。适配文件名、空白 ID 列和输出目录；十个核心函数体与原文的一致性由 `tests/test_r_core.R` 验证。新增 QC、mutation、CN covariation 和 TCGA 模块单独实现。

## Quick Start

先按 [DATA_SETUP](docs/DATA_SETUP.md) 放入数据，并使用 R 4.5.0 恢复包环境：

```r
renv::restore()
renv::status()
```

```bat
Rscript scripts/R/run_analysis.R --mode targeted_dependency --geneA VPS4B --geneB VPS4A
Rscript scripts/R/run_analysis.R --mode genomewide_dependency --geneA VPS4B --geneB VPS4A
Rscript modules/VPS4B_VPS4A/run_case.R
```

`--project` 指定项目根目录，`--bootstrap` 默认 1000，`--min_group_n` 默认 3。通过基因参数切换任意 Gene A / Gene B。主入口拒绝使用其他 R 版本。

| Mode | 分析 |
|---|---|
| qc | 文件、ModelID overlap、CN 分布和组样本数 |
| depmap_cn_expression | CN→同基因表达；表达→Gene B Chronos |
| genomewide_dependency | 所有提供的 Chronos 基因：Pearson、两组 Wilcoxon、BH-FDR、volcano、candidate rank |
| targeted_dependency | Pearson、Spearman、Wilcoxon；带括号、星号和 P 值的组合 PDF |
| lineage_dependency | OncotreeLineage 内中位数差、bootstrap CI 和 forest |
| adjusted_dependency | CN_log 与 CN-Low 的 lineage-adjusted 回归 |
| reverse_dependency | 交换基因后重新定义 CN 分组，比较两个方向 |
| mutation_dependency | damaging / hotspot 分别比较；缺失列不推断为 WT |
| cn_covariation | 所有提供的 CN 基因与 Gene A 的 Pearson 及 BH-FDR |
| tcga_cn_landscape | GDC current gene-level CN 按癌种分布 |
| tcga_cna_prevalence | PanCanAtlas reference GISTIC 五级 CNA prevalence |
| tcga_cn_expression | GDC current 同一样本 CN 与 log2(TPM+1) |
| full | 顺序执行模块，明确记录缺失数据或不足样本导致的跳过 |

## 数据与方法

DepMap **26Q1** 有 13 类本地导出。Chronos 为 1,208 × 18,531 基因，CN 为 1,118 × 18,613 基因，共享 **858 个 ModelID**。尺寸、原始名称、SHA256 和 active 状态见 [manifest](data/manifests/depmap_26Q1_manifest.csv)。文件名含 `subsetted`，尚不能宣称与官方完整 release 覆盖一致；筛选使用所有提供的基因列。

保留 `CN_relative`，计算 `CN_log=log2(CN_relative+1)`。CN-Low `<0.585`；Deep `<0.35`；Shallow `[0.35,0.585)`。这些是 **analysis-defined thresholds**，不是 DepMap 官方 GISTIC 分类。

按用户提供的统计核心，连续 dependency 相关与连续调整回归使用 **CN_log**；CN-expression 和 covariation 使用 **CN_relative**。Delta median = low − non-low，负值表示 CN-Low 更依赖目标。筛选分别对 Pearson/Wilcoxon 作 BH 校正，原始 Rank 按原代码保留；新增 **Eligible_Rank** 只在 Wilcoxon_FDR 非 NA 的基因中按既有排序编号。VPS4A 原始 Rank=276，Eligible_Rank=1，P/FDR/Delta 均未改变。方法和限制见 [ANALYSIS_METHODS](docs/ANALYSIS_METHODS.md)。

```text
data/raw/depmap/26Q1/                    原始导出，内容不变
data/raw/tcga/gdc_current_DR46/          公开 GDC DR46 原始文件
data/raw/tcga/pancanatlas_reference/     独立经典参考层
data/archive/duplicates/                重复导出归档，不删除
data/processed/                         ZSTD Parquet 与 DuckDB views
data/manifests/                         校验、源信息、选择规则、QC
scripts/R/                              分析入口与统计核心
scripts/download/                       发现和官方下载
scripts/preprocess/                     分块预处理
scripts/data_access.py                  单基因查询接口
modules/VPS4B_VPS4A/                     案例配置和入口
results/VPS4B_VPS4A/                     PDF、CSV、Summary、资源记录
docs/                                   设置、方法、来源和准备报告
```

仓库中的 GDC manifest、`config/tcga_layers.json` 和 [当前层报告](docs/TCGA_CURRENT_REPORT.md) 是注明时间的快照。本地下载→预处理→验证作业通过 `scripts/utils/complete_gdc_pipeline.py` 继续，实时状态在忽略目录 `data/raw/tcga/gdc_current_DR46/manifests/live/`，R 优先读取该目录的就绪状态。只有验证通过才标记 complete；current 不可用时不静默切换 reference。案例结果保持实际运行时的状态，经典 GISTIC 分类明确来自 reference。当前来源与缺失清单见 [DATA_SOURCES](docs/DATA_SOURCES.md)。

Python 使用 C 盘现有 Python 3.12，依赖放入项目 `.runtime`。不要调用 Windows Store 的 python/python3 stub。按 [DATA_SETUP](docs/DATA_SETUP.md) 执行发现、下载和预处理；不要并发启动两个写同一 manifest 的下载器。项目移动后重建 DuckDB views。

仅提交代码、配置、文档、manifest 和 VPS4B/VPS4A 示例结果，单文件不得超过 50 MB。代码使用 MIT License；数据使用条款独立适用。

GDC current **DR46 complete**：RNA 11,505/11,505、CN 11,339/11,339；85 项身份/数值/校验检查通过。完整 RNA/CN Parquet 与 DuckDB views 已生成。VPS4B/VPS4A 最新 TCGA 模块结果见 [案例 Summary](results/VPS4B_VPS4A/Summary/VPS4B_VPS4A_Summary.md)。

# 数据设置与复现

固定使用 **R 4.5.0**，运行 `renv::restore()` 和 `renv::status()`。包版本来自 renv.lock。基础 DepMap 分析只需普通 CRAN 包；R arrow/duckdb 为可选项，TCGA 通过 Python 单基因桥接读取。

在 [DepMap Portal](https://depmap.org/portal/download/all/) 自行取得 Public 26Q1 文件或合法导出，放在项目根目录后运行 `scripts/utils/organize_depmap.py`，或按下表放到 `data/raw/depmap/26Q1/`。不提供浏览器验证后的受保护链接。raw 内容不编辑；移动、重命名前后比较 SHA256。同类导出依据模型/基因覆盖选择 active，其余保留到 `data/archive/duplicates/`。

| 内容 | 规范文件名 |
|---|---|
| 模型 | Model.csv |
| 模型条件 / profiles | ModelCondition.csv / OmicsProfiles.csv |
| CRISPR / sequence mapping | CRISPRScreenMap.csv / ScreenSequenceMap.csv |
| WGS relative CN | CopyNumber_WGS_26Q1.csv |
| CRISPR Chronos gene effect | CRISPR_Chronos_26Q1.csv |
| Gene dependency probability | CRISPR_GeneDependency_26Q1.csv |
| Short-read expression | Expression_26Q1.csv |
| Damaging mutation | Mutation_Damaging_26Q1.csv |
| Hotspot mutation | Mutation_Hotspot_26Q1.csv |
| Omics signatures | OmicsSignatures_26Q1.csv |
| Molecular subtypes | MolecularSubtypes_26Q1.csv |

矩阵为 ModelID × gene-symbol 列；空白表头由 fread 读成 V1，在内存识别为 ModelID。矩阵 ModelID 必须唯一。使用所有提供的基因，但 subsetted 导出不代表官方完整覆盖。缺失基因或模型不填作零、WT 或 Diploid。Chronos gene effect 与 dependency probability 区分；expression 沿用导出原始尺度。

TCGA current 在 `data/raw/tcga/gdc_current_DR46/`。从 [GDC Portal](https://portal.gdc.cancer.gov/) / [官方 API](https://docs.gdc.cancer.gov/API/Users_Guide/Downloading_Files/) 查询 TCGA、open、released 的 STAR Counts、Gene Level CN、CN Segments 和 Masked Somatic Mutation。先核查磁盘估计，再下载。逐文件比较 publisher MD5 并计算 SHA256；UUID 子目录避免同名覆盖。

官方 gdc-client 首先尝试，本机其 legacy endpoint TLS 失败后使用官方多 UUID POST API 备用下载器。较大的包曾截断，RNA/CN 使用小批次。只解包清单中的数据成员；已成功 raw 不覆盖，失败和中断会明确记录。

参考层在 `data/raw/tcga/pancanatlas_reference/`，来自完整 Xena CN、thresholded GISTIC、Toil norm-count、phenotype 和 MC3 补充表。参考层与 GDC DR46 独立。Toil 表达是 log2(norm_count+1)，GDC 提供 STAR TPM 和派生 log2(TPM+1)。MC3 Xena 表含 FILTER=PASS，不能当成完整原 MAF。

Git Bash 示例（Windows Python 使用原生绝对脚本路径）：

```bash
project_root='D:/CodexProjects/Copy Number–Driven Dependency Analysis'
/c/Python312/python.exe -m pip install --target "$project_root/.runtime" -r "$project_root/requirements.txt"
/c/Python312/python.exe "$project_root/scripts/download/discover_gdc.py"
/c/Python312/python.exe "$project_root/scripts/download/download_gdc_bundles.py" --batch-size 32
/c/Python312/python.exe "$project_root/scripts/preprocess/preprocess_gdc_clinical.py"
/c/Python312/python.exe "$project_root/scripts/preprocess/preprocess_gdc.py"
/c/Python312/python.exe "$project_root/scripts/preprocess/preprocess_local_depmap.py"
```

processed 由脚本生成，不与 raw 混放。单基因查询读取 requested columns / row groups；DuckDB 仅保存 views。没有 raw 数据不要声称复现完成；新的合法导出可能改变本例样本数和结果。

本项目已有完整 current 数据。继续现有清单时只运行下面的协调入口，不再执行发现或重下全量 RNA/CN。协调入口保留已校验 raw，补齐缺失文件，完成 ETL 和 validation；拒绝竞争协调器。

```bash
/c/Python312/python.exe "$project_root/scripts/utils/complete_gdc_pipeline.py"
```

确认 RNA/CN verified == selected 且 validation 全部通过后，发布快照并仅重跑需要更新的 TCGA 模块：

```bash
/c/Python312/python.exe "$project_root/scripts/utils/publish_gdc_snapshot.py"
/d/R/R-4.5.0/bin/Rscript.exe --version
/d/R/R-4.5.0/bin/Rscript.exe scripts/R/run_analysis.R --project . --mode tcga_cn_landscape --geneA VPS4B --geneB VPS4A
/d/R/R-4.5.0/bin/Rscript.exe scripts/R/run_analysis.R --project . --mode tcga_cn_expression --geneA VPS4B --geneB VPS4A
/d/R/R-4.5.0/bin/Rscript.exe scripts/R/run_analysis.R --project . --mode tcga_cna_prevalence --geneA VPS4B --geneB VPS4A
```

前两项使用 GDC current DR46；最后一项继续使用明确标记的 PanCanAtlas GISTIC reference。无需重跑已完成的 DepMap 分析。

"""Write the project guide and factual case documentation from saved outputs."""
from pathlib import Path
import csv,json,shutil
ROOT=Path(__file__).resolve().parents[2]

def write(path,text):
 p=ROOT/path;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text.strip()+'\n',encoding='utf-8')

def main():
 old=ROOT/'docs/INITIAL_PREPARATION_README.md'
 if not old.exists():shutil.copy2(ROOT/'README.md',old)
 write('README.md','''
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
| cn_threshold_sensitivity | CN_log 0.585/0.50/0.40/0.35 的分组敏感性；连续相关只算一次 |
| expression_dependency | 独立表达连续相关与 bottom 10%/20% Wilcoxon |
| genomewide_expression_dependency | 全 Chronos 的表达定义筛选；仅显式运行 |
| genomewide_adjusted_dependency | 全 Chronos 的 CN_log/CN-Low + lineage 回归；仅显式运行 |
| full | 顺序执行模块，明确记录缺失数据或不足样本导致的跳过 |

## 数据与方法

DepMap **26Q1** 有 13 类本地导出。Chronos 为 1,208 × 18,531 基因，CN 为 1,118 × 18,613 基因，共享 **858 个 ModelID**。尺寸、原始名称、SHA256 和 active 状态见 [manifest](data/manifests/depmap_26Q1_manifest.csv)。文件名含 `subsetted`，尚不能宣称与官方完整 release 覆盖一致；筛选使用所有提供的基因列。

保留 `CN_relative`，计算 `CN_log=log2(CN_relative+1)`。CN-Low `<0.585`；Deep `<0.35`；Shallow `[0.35,0.585)`。这些是 **analysis-defined thresholds**，不是 DepMap 官方 GISTIC 分类。

按用户提供的统计核心，连续 dependency 相关与连续调整回归使用 **CN_log**；CN-expression 和 covariation 使用 **CN_relative**。Delta median = low − non-low，负值表示 CN-Low 更依赖目标。筛选分别对 Pearson/Wilcoxon 作 BH 校正。用户默认查看 **Eligible_Rank / Eligible_N**：只在 Wilcoxon_FDR 非 NA 的基因中按 FDR、Delta_median 升序排序。VPS4A Eligible_Rank=1/18,256；历史 Rank=276 保留用于 provenance，P/FDR/Delta 均未改变。方法和限制见 [ANALYSIS_METHODS](docs/ANALYSIS_METHODS.md)。

`full` 默认加入 CN threshold sensitivity 和 expression dependency；两个新增 genome-wide 模块仅在显式指定 mode 时运行。它们读取现有 processed DepMap Parquet，独立于原始 CSV 统计核心。TCGA CN-expression 保留 overall Pearson/Spearman，另输出 N>=20 的每癌种相关、Fisher-z CI，以及原尺度/标准化 cancer-adjusted 回归。**Pan-cancer overall correlation may be influenced by between-cancer differences.** Summary 同时报告 cancer-adjusted CN beta。

新增模块的 synthetic tests 与原始 core 函数测试由 [GitHub Actions CI](.github/workflows/ci.yml) 在 R 4.5.0 下执行；CI 仅恢复 renv 锁定的测试依赖，不下载 TCGA/DepMap 数据。[官方 R Actions](https://github.com/r-lib/actions) 提供 R 环境设置。

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
''')
 write('docs/DATA_SETUP.md','''
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
''')
 write('docs/ANALYSIS_METHODS.md','''
# Statistical methods

The user-supplied R source is retained verbatim, with SHA256 and section mapping in R_CORE_PROVENANCE.json. tests/test_r_core.R compares ten core function bodies to that source. Adapters recognize V1 as the blank CSV ID, route folders and restore state after reverse analysis. Formulae, thresholds, correlations, Wilcoxon calculations, BH families and ranking remain unchanged.

CN_relative is the supplied WGS value; CN_log=log2(CN_relative+1). Analysis-defined low is <0.585, deep <0.35, shallow [0.35,0.585), non-low >=0.585. These are not DepMap official GISTIC states. Finite observations join by exact ModelID. OncotreeLineage defines stratification/adjustment; disease, subtype and code remain available.

CN-expression uses relative CN and original expression. Targeted CN-dependency uses CN_log, matching the supplied code. Pearson and Spearman (exact=FALSE) are reported. Two-sided Wilcoxon uses exact=FALSE and R default continuity correction. Chronos is gene effect, with more negative scores indicating greater dependency. Three-group comparison requires all groups >=3; optional pairwise tests report raw P and separate BH-FDR.

The genome-wide function tests every supplied Chronos gene column. Pearson requires >=10 finite models, uses the correlation t statistic (df=n−2), and retains the supplied undefined P at absolute r=1. Group effects require each group >=min_group_n. Delta mean/median are low minus non-low. Pearson and Wilcoxon have separate complete-screen BH families. Rank sorts Wilcoxon FDR then delta median, retaining the supplied NA ordering for provenance. Eligible_Rank is the default user-facing rank: only non-NA Wilcoxon_FDR rows, sorted by FDR ascending then Delta_median ascending (original Rank breaks exact ties). Eligible_N records that complete eligible family size. Excluded rows have blank Eligible_Rank. The adapter adds these fields without changing any original numeric CSV tokens, verified against the initial committed CSVs. Neither rank is a causal-priority score.

Lineage groups each require >=3. Bootstrap independently resamples both groups and takes percentile 95% CI for median-low minus median-nonlow, default 1000 draws and seed 1234. BH applies to eligible lineages. Forest magnitude represents effect size, not P value.

Adjusted models exactly match the source: Chronos ~ CN_log + OncotreeLineage; Chronos ~ I(CN_binary == "CN-Low") + OncotreeLineage. All coefficients, estimates, SE, t and P are retained. Reverse swaps genes, rebuilds CN groups and applies the same targeted statistics before restoring the original gene/root state.

Additional mutation modes define positive finite mutation values as Mutant, zero as WT, and exclude missing values. Damaging and hotspot remain separate. Missing columns do not imply WT; groups below three skip. BH applies to eligible damaging/hotspot comparisons. CN covariation uses relative CN, Pearson (>=10 finite values), and BH across supplied CN genes; self-correlation remains in CSV with undefined r=1 P, excluded from the top-other-gene plot.

TCGA current uses source total gene CN and STAR log2(TPM+1). For sample-level CN-expression, one RNA aliquot is selected per exact sample UUID: prefer the chosen CN aliquot, then lexical RNA file UUID; selection and excluded alternatives are audited. Raw aliquots remain intact. Five-state prevalence uses explicitly labeled reference GISTIC, with denominators of tumor samples having a finite gene value. These are sample-level prevalences. Aliquot/sample/file provenance remains available; no patient-level average is silently created. Actual current validation and readiness are recorded in TCGA_CURRENT_REPORT.md and gdc_validation_report.json.

Observational associations do not establish a causal synthetic-lethal mechanism or clinical benefit. Lineage adjustment addresses measured lineage differences; CN covariation does not establish chromosomal adjacency or causality. Insufficient groups, absent genes and unfinished current TCGA are explicitly recorded.
''')
 write('LICENSE','''
MIT License

Copyright (c) 2026 Copy Number–Driven Dependency Analysis contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Dataset licenses and access terms are independent of this software license.
''')
 write('docs/R_ENVIRONMENT_NOTES.md','''
# R environment validation

R 4.5.0 is required. Renv is project-local; cache, sandbox and library paths are inside the project. Windows UTF-8 LC_CTYPE supports the en dash in the project name. Required package versions are in data/manifests/R_package_versions.csv and renv.lock.

The execution environment initially omitted PROCESSOR_ARCHITECTURE. Native cli shutdown called strcmp(getenv("PROCESSOR_ARCHITECTURE"), "ARM64") with NULL, causing Windows access violation 0xc0000005 in ucrtbase.dll after valid R outputs. Minimal probes and an owned-process debugger isolated the fault. Recompilation and a clean official R runtime did not remove the environment cause. Windows GetNativeSystemInfo returned architecture code 9 (AMD64). .Rprofile and the analysis entry point now fill this missing variable from the verified executable architecture before packages load. Original locked cli/rlang versions were restored; the final probe exited normally.

Some current CRAN Windows binaries report that they were built under R 4.5.3; the actual interpreter remains 4.5.0. Successful package loading, actual analysis exit codes and renv::status provide runtime validation. Earlier failed native probes/rebuilds are diagnostic history, not successful case runs.
''')
 summary=ROOT/'results/VPS4B_VPS4A/Summary/VPS4B_VPS4A_Summary.md'
 if summary.exists():
  write('modules/VPS4B_VPS4A/README.md','''
# VPS4B → VPS4A case

Research question: Does VPS4B copy-number loss associate with increased dependency on VPS4A?

R 4.5.0; user-supplied DepMap 26Q1 exports. The wrapper reads the main framework rather than duplicating statistical code.

```bat
Rscript modules/VPS4B_VPS4A/run_case.R
Rscript modules/VPS4B_VPS4A/run_case.R targeted_dependency
```

The workflow is QC → CN-expression → complete supplied Chronos screen → VPS4A candidate rank → targeted dependency → lineage forest → adjusted regression → reverse comparison. Mutation/covariation and explicitly separated TCGA layers are also supported. Results below come from saved outputs, with no preset conclusion.

'''+summary.read_text(encoding='utf-8'))
 write('modules/VPS4B_VPS4A/expected_outputs.md','''
# Output map

Results are under results/VPS4B_VPS4A/. Each core module has its named subdirectory:

- 00_QC/: CN_Distribution.pdf, CN_Group_Counts.csv, QC_Statistics.csv, model alignment.
- 01_CN_Expression/: VPS4B_CN_vs_Expression.pdf, CN_Expression_Statistics.csv.
- 02_GenomeWide_Dependency/: complete GenomeWide_Dependency.csv, volcano PDF, Top100, VPS4A_Candidate_Rank.csv.
- 03_Targeted_Dependency/: VPS4B_VPS4A_Statistics.csv, Group_Summary.csv, CellLines.csv and combined PDF; three-group plot only if every group >=3.
- 04_Lineage_Dependency/: VPS4B_VPS4A_Lineage_Dependency.csv and Lineage_Forest.pdf.
- 05_Adjusted_Dependency/: Continuous_CN_Adjusted.csv, CNLow_Adjusted.csv.
- 06_Reverse_Dependency/: Reverse_Dependency_Summary.csv and nested reverse-targeted outputs.
- 07_CN_Covariation/: VPS4B_CN_Covariation.csv and Top_CN_Covariation.pdf.
- 08_TCGA/: current CN landscape/expression if prepared; separately labeled reference GISTIC prevalence.
- 09_Mutation_Dependency/: eligible damaging/hotspot comparisons, or explicit skip reasons.
- 10_CN_Threshold_Sensitivity/: threshold table/PDF and once-only continuous statistics.
- 11_Expression_Dependency/: continuous and quantile statistics, matched models and three-panel PDF.
- 12_GenomeWide_Adjusted_Dependency/: full adjusted screen, optional candidate and volcano; explicit mode only.
- 13_GenomeWide_Expression_Dependency/: separate full expression screen; explicit mode only.
- Summary/: VPS4B_VPS4A_Summary.md, Key_Statistics.csv, Analysis_Metadata.json, Module_Runs.csv, Resource_Monitor.json and SessionInfo.txt.

Missing data or inadequate group sizes produce documented skips; they never produce fabricated values. Resource_Monitor records observed native process exits and approximate RSS, sampled every 0.5 s. No file above 50 MB is uploaded.
''')
 append_extensions_docs()

 layers=json.loads((ROOT/'config/tcga_layers.json').read_text(encoding='utf-8'))
 validation=ROOT/'data/manifests/gdc_validation_report.json'
 if layers.get('current_status')=='complete' and validation.exists():
  report=json.loads(validation.read_text(encoding='utf-8'))
  if report.get('all_passed'):
   path=ROOT/'README.md'
   with path.open('a',encoding='utf-8') as f:
    f.write(f"\nGDC current **DR46 complete**：RNA 11,505/11,505、CN 11,339/11,339；{len(report['checks'])} 项身份/数值/校验检查通过。完整 RNA/CN Parquet 与 DuckDB views 已生成。VPS4B/VPS4A 最新 TCGA 模块结果见 [案例 Summary](results/VPS4B_VPS4A/Summary/VPS4B_VPS4A_Summary.md)。\n")

def append_extensions_docs():
 methods='''

## Additional independent modules

CN threshold sensitivity retains the main analysis-defined 0.585 cutoff, and separately evaluates 0.50, 0.40 and 0.35 on CN_log=log2(relative CN+1). Low is strictly below each threshold; equality is non-low. Only finite CN/Chronos pairs enter. Means/medians and low-minus-nonlow effects are recorded even when groups are too small; Wilcoxon requires both N>=3, otherwise P/FDR are NA and eligible=false. BH is across eligible thresholds. Continuous CN_log Pearson/Spearman is computed once and stored separately. These are analysis-defined thresholds, not official DepMap GISTIC classifications.

Expression-defined dependency is independent of CN loss and mutation. Continuous correlation uses Gene A's supplied expression scale and Gene B Chronos. Targeted cutoffs are R type-7 quantiles of finite matched expression/Chronos pairs at 10% and 20%; low is expression <= cutoff, and missing expression is excluded rather than made non-low. Ties may increase the nominal percentile group size. Wilcoxon requires both groups N>=3, with BH across eligible quantiles. In an explicit genome-wide expression screen, cutoffs are fixed on the finite expression cohort aligned to all Chronos ModelIDs, then each gene's missing Chronos is excluded; continuous Pearson BH and the two grouped Wilcoxon BH families are separate. Candidate-specific missingness can therefore change its targeted vs genome-wide expression cutoff; no hidden cutoff tuning occurs.

TCGA CN-expression retains exactly the previously validated sample UUID join and overall results, using current GDC DR46 gene-level CN and log2(TPM+1). Per-cancer analysis requires N>=20 and variable CN/expression; ineligible cancers remain in CSV with NA tests. Pearson and Spearman each have a separate BH family across eligible cancers. Pearson 95% CI is tanh(atanh(r) +/- qnorm(0.975)/sqrt(N-3)); exact +/-1 uses a degenerate boundary interval. Forest plots show eligible cancer types sorted by Pearson r. Missing CancerType is excluded from per-cancer and adjusted analyses, with regression N retained.

Pan-cancer overall correlation may be influenced by between-cancer differences. Report overall correlation alongside the CN coefficient from Expression ~ CN + CancerType, with cancer treated as a factor. Standardized regression z(Expression) ~ z(CN) + CancerType uses global sample means/SDs on the same complete regression cohort; this is not within-cancer standardization. CSVs retain all estimable terms and N. These are expression associations, not CRISPR dependency or causal effects. GISTIC prevalence remains the separately labeled PanCanAtlas reference layer.

Explicit genome-wide lineage adjustment fits each Chronos gene on CN_log + OncotreeLineage and separately on I(CN_log<0.585) + OncotreeLineage. Exact ModelID alignment is enforced; missing outcome/CN/lineage is removed per gene. N>=20 and at least two observed lineage levels are required. Constant outcomes or aliased CN terms remain skipped/NA rather than assigned a test. Continuous and binary CN P values have separate complete-screen BH families. Eligible_Rank_adjusted ranks non-NA FDR_CN by FDR_CN, Beta_CN, then Gene ascending; Eligible_N_adjusted records that continuous family size. This rank belongs to the adjusted model and is distinct from the original Wilcoxon Eligible_Rank. Volcano displays Beta_CN vs -log10(FDR_CN), with Gene B highlighted when supplied. Both new genome-wide modules are excluded from default full.

Original outputs are retained in Git history and copied locally under .runtime/prior_case before the requested full rerun. tests/validate_preserved_results.py verifies original CSV tokens against the completed prior commit; tests/validate_eligible_rank.py additionally verifies the initial screen's original P/FDR/effects/Rank. Statistical core function bodies remain unchanged.
'''
 setup='''

## Extension modules and tests

No new TCGA or DepMap download is needed for these modules. New DepMap modes read existing Parquet under data/processed/depmap/26Q1; the original validated core retains its original CSV reader. TCGA extensions reuse existing current matched samples and processed DR46 matrices.

```bat
Rscript --vanilla scripts/R/verify_environment.R
Rscript --vanilla tests/run_synthetic_tests.R
Rscript --vanilla tests/validate_case_results.R VPS4B VPS4A
Rscript scripts/R/run_analysis.R --mode full --geneA VPS4B --geneB VPS4A
Rscript scripts/R/run_analysis.R --mode genomewide_adjusted_dependency --geneA VPS4B --geneB VPS4A
Rscript scripts/R/run_analysis.R --mode genomewide_expression_dependency --geneA VPS4B --geneB VPS4A
Rscript scripts/R/run_analysis.R --mode cn_threshold_sensitivity --geneA ENO1 --geneB ENO2 --output_case ENO1_ENO2_Method_Test
Rscript scripts/R/run_analysis.R --mode expression_dependency --geneA ENO1 --geneB ENO2 --output_case ENO1_ENO2_Method_Test
```

The two genome-wide extension commands are opt-in; full includes only the targeted sensitivity/expression extensions. ENO1/ENO2 is a method test without a required positive result. --output_case only accepts a single safe result-directory name. A genomewide_adjusted_dependency invocation without --geneB emits the complete screen/volcano without a candidate table.

The Windows GitHub Actions job installs R 4.5.0 and restores only data.table/dplyr/ggplot2 and their dependencies from the existing renv.lock. All six CI tests use synthetic fixtures or parse the supplied R function definitions; none reads raw or processed biological data. The CI environment intentionally omits optional arrow/duckdb because those libraries are only needed for actual local data analyses.
'''
 for name,value in [('docs/ANALYSIS_METHODS.md',methods),('docs/DATA_SETUP.md',setup)]:
  with (ROOT/name).open('a',encoding='utf-8') as f:f.write(value)

if __name__=='__main__':main()

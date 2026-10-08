# Copy Number–Driven Dependency Analysis

**R 4.5.0** 科研分析框架，检验 Gene A 拷贝数与 Gene B CRISPR dependency 的关联。普通分析 **DATA_MODE=local**，只读现有 processed Parquet；数据库维护与分析分开，不为重新分析一个基因扫描 raw 或下载数据库。

## 使用

准备好本地 processed 数据与 R 包环境后，在项目根目录运行：

```bat
Rscript scripts/R/run_analysis.R --workflow geneA_screen --geneA VPS4B
Rscript scripts/R/run_analysis.R --workflow geneA_geneB --geneA VPS4B --geneB VPS4A
```

只有 Gene A：TCGA CNA、CN 分布、逐癌种自身 CN–mRNA，接着 DepMap 全基因 dependency 与 CN covariation。A+B：先 TCGA，再 targeted 三图与三个补充分析，最后全基因筛选和 CN covariation。

缺少本地 processed 数据立即报错，不联网补齐。有效缓存默认跳过计算；`--force` 强制重算。默认不运行 reverse、mutation、expression dependency、genomewide expression 或 genomewide adjusted。已有 `--mode` 仅用于显式开发/单模块调用，普通用户使用上述 workflow。案例入口 `Rscript modules/VPS4B_VPS4A/run_case.R` 默认使用 A+B workflow。

## 结果

```text
results/<GeneA>[_<GeneB>]_Analysis/
  00_Analysis_Summary.txt
  Main_Results/       A-only 01–05；A+B 01–08：每个编号独立文件夹，PDF及对应CSV放在一起
                      A-only 37 PDF；A+B 40 PDF，均含33癌种 CN–mRNA
  Supplementary/      仅A+B：Lineage、Adjusted、CN Threshold Sensitivity；CSV在对应模块内
  Provenance/         数据来源、版本/哈希、sessionInfo、运行与验证记录
    Data_Audit/       公共TCGA样本、RNA代表选择和baseline明细；无集中Tables目录
```

[VPS4B/VPS4A 中文结果指南](results/VPS4B_VPS4A_Analysis/00_Analysis_Summary.txt) · [最终工作流与统计定义](docs/FINAL_WORKFLOWS.md)

## 数据与统计

TCGA 主图 01–03 默认使用本地 **current DR46 processed dataset**：`data/processed/tcga/gdc_DR46/`。01 为 current absolute gene-level CN 与 sample-specific baseline 的 analysis-defined 五分类；02 为 current gene-level CN；03 为 current CN 与 STAR `log2(TPM + 1)`。来源为 NCI GDC，Release 为 DR46 / 46.0，完整记录于 Provenance。PanCanAtlas/Xena 仅保留为 historical/reference，默认 workflow 不使用。

当前 processed 数据没有可用的实测 ploidy；baseline 取每个样本常染色体整数 CN 的众数（并列取较小值），分块读取并缓存样本级摘要，独立抽样验证。五分类为 CN=0、0<CN<baseline、CN=baseline、baseline<CN<2×baseline、CN≥2×baseline，依次标为 Deep Deletion、Shallow Deletion、Diploid、Gain、Amplification；它们不是官方 GISTIC 五级值，baseline 也不是实测 ploidy。

所有 TCGA 主图按 Gene A finite continuous CN 的癌种 median 从小到大排列，视觉从上到下；CN landscape 展示每个有限样本点和实际 N，保留极端值。33癌种分别分析自身 CN–mRNA，固定 N<20 或常量变量时保留图并不编造相关统计。

DepMap 使用全部提供的 Chronos 基因列，用户导出相对官方完整 release 的完整性未验证。保留原始十个函数体，测试对照 [用户提供的统计核心](docs/USER_SUPPLIED_ANALYSIS.R)。运行时 I/O 适配器只读 Parquet，最终 workflow 复用统计前缀并独立制图。

`CN_log=log2(relative CN+1)`，CN-low `<0.585`，Deep `<0.35`，均为分析定义，不是官方 GISTIC 分类。`Delta_median=median(low)-median(nonlow)`：负值表示 CN-low 组依赖更强。全筛选 Wilcoxon/BH-FDR/效应与历史数值保持一致；主排名 `Eligible_Rank` 排除 undefined FDR，历史 Rank 仅用于复核。Targeted scatter 使用 relative CN，并明确保留表中的历史 CN_log Pearson/Spearman。**所有结果描述相关或依赖差异，不作因果解释。**

## 验证与维护

```bat
Rscript --vanilla tests/run_synthetic_tests.R
C:/Python312/python.exe tests/test_git_size_guard.py
C:/Python312/python.exe scripts/utils/validate_final_sources.py --case VPS4B_VPS4A_Analysis
```

普通工作流内置验证排序、五级计数、点数、相关独立复算、BH、组别和文件索引；每个模块缓存包括输入/代码/参数和输出哈希。数据库下载/重新预处理只通过 [独立维护命令](docs/FINAL_WORKFLOWS.md#独立数据库维护) 或既有 scripts/download、scripts/preprocess 执行。初始环境准备见 [DATA_SETUP](docs/DATA_SETUP.md)，数据库来源见 [DATA_SOURCES](docs/DATA_SOURCES.md)。

## GitHub

每次任务先 pull，完成代码或分析后必须 commit/push main；仅上传代码、文档、测试和小型结果。结果目录对所有基因开放跟踪。使用 `scripts/utils/git_size_guard.py stage` 在 staging 前拒绝 >50 MB、raw/processed/runtime/cache/数据库或凭据；pre-commit 检查实际 index。安装 hook：`git config core.hooksPath .githooks`。代码 push 后等待 [R synthetic validation](.github/workflows/ci.yml) success，并确认工作区 clean 与本地/远程 SHA 一致。纯结果更新不触发 CI。

代码为 MIT License；数据使用条款独立适用。

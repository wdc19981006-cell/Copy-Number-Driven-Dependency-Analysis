# 最终用户工作流

在项目根目录运行，要求 R **4.5.0** 和已经准备好的本地 processed 数据。普通分析的 `DATA_MODE=local`，不调用更新、下载、GDC API，也不读取 raw。缺少 processed 数据时直接报错：

```text
Local processed dataset unavailable.
Run the separate data-update/preparation workflow first.
```

```bat
Rscript scripts/R/run_analysis.R --workflow geneA_screen --geneA VPS4B
Rscript scripts/R/run_analysis.R --workflow geneA_geneB --geneA VPS4B --geneB VPS4A
```

`--force` 强制重新计算必要模块。`--project` 指定项目路径；`--output_case` 可指定安全的单层目录名。默认 `--min_group_n=3`、`--bootstrap=1000`，不为预期结果调整阈值。

```text
results/<GeneA>[_<GeneB>]_Analysis/
  00_Analysis_Summary.txt
  Main_Results/
    01_TCGA_CNA_Percentage/          PDF + TCGA_CNA_Percentage.csv
    02_TCGA_CopyNumber_Landscape/    PDF + TCGA_Cancer_Order.csv
    03_TCGA_CN_mRNA/                33 癌种 PDF + 汇总统计 + CNA counts
    04_DepMap_CN_vs_Dependency/     仅 A+B：PDF + targeted 统计 + matched cell lines
    05_DepMap_Dependency_Waterfall/ 仅 A+B：PDF + waterfall order
    06_DepMap_Dependency_CNlow_vs_Normal/ 仅 A+B：PDF + 已有 targeted 统计的分组字段
    04或07_DepMap_GenomeWide_Dependency/  PDF + 全筛选/候选 CSV + 编号_Analysis_Status.txt
    05或08_DepMap_CN_Covariation/    PDF + 全共变/Top CSV
  Supplementary/                 仅 A+B
    01_Lineage/
    02_Adjusted/
    03_CN_Threshold_Sensitivity/
  Provenance/                    来源、版本、代码和输入哈希、验证和运行记录
    Data_Audit/                  TCGA_Current_Samples.csv、TCGA_RNA_Representative_Selection.csv、TCGA_Sample_Baselines.csv
```

Gene A-only 共 **37** 张主 PDF；A+B 共 **40** 张主 PDF。默认不运行 reverse、mutation、expression dependency、genomewide expression 或 genomewide adjusted。研究反向依赖时交换 A/B 后运行同一个 A+B workflow。

每个主编号文件夹包含该图及其 CSV，Supplementary 的 CSV 也直接保存在所属模块内；不再生成集中 `Tables/`，不复制同一结果表。06 的分组表仅提取现有 targeted 结果字段，不新增检验。所有 PDF 的原文件名保持不变。

目录迁移使用 `scripts/R/migrate_result_layout.R` 的 `migrate_result_layout(previous_workflow, previous_file_sources)`，在已初始化的通用工作流环境调用；旧工作流源文件保存在结果包外。旧 TCGA 导出助手的文件快照通过具名列表 `list("scripts/utils/export_workflow_tcga.py"="旧源文件快照路径")` 传入，逐字核对只改变 baseline JSON 的输出目的地。迁移先验证旧缓存代码/输入签名与输出哈希，确认统计和绘图函数不变，记录 SHA256 后移动，再独立核对字节并更新缓存路径与签名。若旧签名不能验证，报告受影响模块并停止，不运行分析。`File_Index.csv` 按授权更新，其旧哈希保存在迁移审计 JSON；原始 PDF、结果/样本 CSV 和 Module_Runs.csv 的字节保持一致。原统计验证记录保留，另附结构迁移验证；不声称重新执行统计验证。

## TCGA 数据层与共同顺序

TCGA 主图 01–03 **统一使用本地 current DR46 processed dataset**，只读取 `data/processed/tcga/gdc_DR46/`。Database=TCGA，Source=NCI GDC，Release=DR46 / 46.0。01 保持 analysis-defined current CN categories；02 和 03 使用 `TCGA_Relative_CN_Change = CopyNumber / BaselineCN - 1`，直接复用现有 sample-specific baseline。03 的表达仍为 STAR `log2(TPM + 1)`。这是 analysis-derived relative copy-number change，不是 GISTIC score、log2 GISTIC 或 DepMap CN_relative。PanCanAtlas/Xena 仅作为 historical/reference 保留，不在默认 workflow 中。实际输入、prepared SHA256 与数据层记录于 Provenance。

01 使用所有有效分类的 current tumor samples 作为各癌种分母，五类颜色固定。当前没有可靠实测 ploidy，使用常染色体 gene-level CN 整数众数估算 sample-specific baseline（并列取较小值；不舍入 CN；不统一固定为 2）。按 CN=0、0<CN<baseline、CN=baseline、baseline<CN<2×baseline、CN≥2×baseline 分成 Deep Deletion、Shallow Deletion、Diploid、Gain、Amplification。此为项目分析定义，不是 PanCanAtlas GISTIC 五级值。baseline 在 bounded gene-column blocks 中计算，样本级缓存绑定源版本及哈希；普通单基因提取只读必要列。`TCGA_Sample_Baselines.csv` 和 `TCGA_Baseline_Method.json` 保留分布、众数支持度、并列与极端值。

02 使用全部 finite current CN tumor samples，箱线图叠加全部 jitter points，保留极端值，癌种标签显示实际 N。横轴以 0 为中心，增加 0 的垂直虚线：0 为样本自身 baseline，负值为 relative loss，正值为 relative gain。01/02 共同癌种顺序按 `median(Relative_CN_Change)` 升序（并列按癌种名）；`TCGA_Cancer_Order.csv` 是**视觉从上到下**的顺序，coord_flip 的 factor levels 为该表倒序。表中保留 Median_Absolute_CN 和 Median_Relative_CN_Change，历史 Median_CN 字段仍表示 absolute CN。Provenance 保存实际绘图 factor levels、点数和转换规则，自动验证与表一致。

`TCGA_Current_Samples.csv` 保留原 CopyNumber、BaselineCN、CNAState、RNA_TPM、Expression，并新增 TCGA_Relative_CN_Change 和导出别名 Relative_CN_Change。自动校验五类分别对应 -1、(-1,0)、0、(0,1)、[1,+Inf)，任何不一致均使 workflow 失败。此转换不改变原有 baseline、样本纳入或 CNA 分类。

03 精确按 SampleID 一对一匹配 current CN、analysis-defined CNA、STAR TPM，限定 tumor，使用三者都 finite 的样本。RNA 同样本多 aliquot 时先匹配 CN aliquot，再按 FileID 字典序选代表；检查 CaseID/ProjectID 一致并保留完整选择审计。横轴、回归线与 Pearson/Spearman 均使用 `TCGA_Relative_CN_Change`，增加 x=0 垂直虚线；表达仍为 `log2(TPM + 1)`。统计表记录 `CN_metric=Relative_CN_Change` 和同一匹配 cohort 的 absolute/relative CN median。逐癌种报告 N、Pearson/Spearman 和各自 eligible tests 内 BH FDR；固定 N<20 或常量变量时不编造统计量，仍生成全部33癌种图。左上角为真正 five-state pie chart，分母是同图匹配 cohort；小于8%的扇区省略文字，但图例始终保留5类。图上方预留空白，使 pie 和统计不遮挡散点。

## DepMap 统计定义

使用本地 26Q1 WGS CN 和 Chronos Parquet：单基因读必要列；genome-wide 和 CN covariation 各读取一次必要完整矩阵。用户提供 portal exports 的完整性相对官方全 release 未验证，使用其全部提供的 Chronos 基因列。

保留原始数值核心与 `CN_log=log2(relative CN+1)`、low `<0.585`。Genome-wide `Delta_median=median(low)-median(nonlow)`：负值表示 low 组更依赖，正值表示更弱。Wilcoxon P 和全筛选 BH FDR 不更改；历史 Rank 保留，主排名 Eligible_Rank 排除 undefined FDR。火山图标记每个方向 Top10，并额外高亮提供的 Gene B；候选表每方向 Top20。

Targeted scatter 的横轴和角落统计是 **CN_relative**；为保留历史定义，统计表另保留 CN_log Pearson/Spearman，不把两种 Pearson 默认为相同。Waterfall 按 Chronos 从高到低（弱到强）排序，并列按 ModelID；CN-Normal 为黑色，low 红色。Boxplot 展示全部点、N、组 median、Delta median、Wilcoxon P 和星号。

04 包含 y=0 灰色水平点线、CN-Low threshold=`2^0.585-1` 深灰垂直虚线、Deep CN Loss=`2^0.35-1` 黑色垂直点线。Deep 仅为视觉参考，分组仍只有 CN-Low/CN-Normal。只有散点图函数变化时直接复用已验证统计和05/06，只重新生成04。

Covariation 报告 positive Top20 与 negative Top20，排除 Gene A 自身；底部说明 Pearson r 是 Gene A 与各基因 CN 在 DepMap 细胞系中的线性相关系数，r > 0 倾向同向，r < 0 倾向反向，|r| 越接近 1 相关越强。panel 标题使用 usually，相关不解释为因果关系。TCGA 相对变化字段不进入任何 DepMap 分组或 dependency 统计。

Genome-wide 模块先以原有 CN_log 分组匹配 Chronos ModelID。如果任一组 n < min_group_n（默认 3），模块状态为 **SKIPPED**，不读取完整 Chronos 矩阵、不进行组间筛选；TCGA 和 CN covariation 继续完成。固定火山图文件保留跳过原因和样本量说明，统计/候选 CSV 仅保留表头；`Provenance/GenomeWide_Dependency_Status.json` 与 Module_Runs 记录状态、原因和原阈值 0.585。缓存仍校验这些输出，缓存命中不会把 SKIPPED 误记为已完成筛选。足量分组继续执行原统计核心。

## 缓存与验收

每个必要模块只调用一次。缓存键包含 Gene A/B、输入版本及已有 SHA256、文件大小/mtime、参数哈希、涉及的统计/绘图函数体和文件哈希。缓存记录各输出 SHA256；只有参数/输入/代码一致且所有输出存在、结构有效、内容哈希一致才跳过。修改某模块的绘图/统计函数只使相应模块失效，共用数据访问或阈值变动会使依赖它的模块失效。Provenance 的 `Cache_*.json` 是小型校验 metadata，不含中间矩阵。

运行后内置验证会检查五级计数、排序、全部点数、样本匹配、逐癌种相关独立重算/BH、targeted 统计、waterfall 和完整文件索引；额外的 processed-source 验证：

```bat
C:/Python312/python.exe scripts/utils/validate_final_sources.py --case VPS4B_VPS4A_Analysis
Rscript --vanilla tests/run_synthetic_tests.R
C:/Python312/python.exe tests/test_git_size_guard.py
```

旧 VPS4B/VPS4A 目录仅在新结果验证通过后归档/移除。默认普通分析不用 `full`，该 mode 仅供开发。

## 独立数据库维护

以下命令不被分析入口导入或调用。本轮不执行任何下载或重新预处理。

```bat
C:/Python312/python.exe scripts/maintain_data.py update_tcga
C:/Python312/python.exe scripts/maintain_data.py update_depmap
C:/Python312/python.exe scripts/maintain_data.py prepare_tcga_reference
C:/Python312/python.exe scripts/maintain_data.py prepare_depmap_local
C:/Python312/python.exe scripts/maintain_data.py validate_data
```

维护命令使用现有 downloader / preparer；数据库更新可能要求官方 release/catalog 等配置，遵循 DATA_SETUP。`validate_data` 属于数据库维护，可扫描 raw/hash；普通分析只做小 cohort 和结果验证。

## GitHub 跟踪

每次开始 `git pull --ff-only origin main`。结束使用显式路径 staging，先拒绝任何 >50,000,000 bytes、raw/processed/runtime/cache、数据库扩展名或凭据模式：

```bat
C:/Python312/python.exe scripts/utils/git_size_guard.py stage scripts tests config .github README.md docs/FINAL_WORKFLOWS.md results/<case>
C:/Python312/python.exe scripts/utils/git_size_guard.py audit
git commit -m "Analysis <case>"
git push origin main
git status --short
git rev-parse HEAD
git ls-remote origin refs/heads/main
```

`stage` 是有大小检查的统一 staging 入口；Git 本身没有 pre-stage hook，因此另安装 `git config core.hooksPath .githooks`，pre-commit 再检查实际 index blob，即使 worktree 在 staging 后变化也不会漏检。raw/processed、Parquet/DuckDB、临时 cache 不上传；超大表仅保存 top table 和 manifest。

CI 只由 scripts/tests/.github/renv.lock/.R-version/config 变化触发；纯 results 与分析文档更新不触发。代码 push 后检查 `R synthetic validation` success，确认本地与 remote main SHA 一致，工作区 clean；push 失败必须明确报告 `GitHub push failed`。

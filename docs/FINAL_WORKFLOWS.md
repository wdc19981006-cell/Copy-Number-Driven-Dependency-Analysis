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
    01_TCGA_<GeneA>_CNA_Percentage.pdf
    02_TCGA_<GeneA>_CopyNumber_Landscape.pdf
    03_TCGA_CN_mRNA/              33 癌种，各一个 PDF
    04..06_DepMap_...pdf          仅 A+B：scatter / waterfall / boxplot
    04或07_DepMap_<GeneA>_GenomeWide_Dependency_Volcano.pdf
    05或08_DepMap_<GeneA>_CN_Covariation.pdf
  Supplementary/                 仅 A+B
    01_Lineage/
    02_Adjusted/
    03_CN_Threshold_Sensitivity/
  Tables/                        所有统计 CSV、匹配样本与绘图输入
  Provenance/                    来源、版本、代码和输入哈希、验证和运行记录
```

Gene A-only 共 **37** 张主 PDF；A+B 共 **40** 张主 PDF。默认不运行 reverse、mutation、expression dependency、genomewide expression 或 genomewide adjusted。研究反向依赖时交换 A/B 后运行同一个 A+B workflow。

## TCGA 数据层与共同顺序

TCGA 主图 01–03 **统一使用现有 reference 层**，不与 DR46 主图混排。五级 CNA 来自 PanCanAtlas/Xena thresholded GISTIC2，连续 CN 为其匹配的 GISTIC2 all-data，表达是现有 Toil TCGA RSEM `log2(norm_count+1)`，不是 TPM。保留真实来源与现有 prepared SHA256 到 Provenance。当前 NCI GDC DR46 processed 数据保持独立可查询；普通主图不增加第二套 CN-expression 输出。

01 使用所有有效 GISTIC tumor samples 作为各癌种分母，五类颜色固定。02 使用全部 finite reference CN tumor samples，箱线图叠加全部 jitter points，保留极端值，癌种标签显示实际 N。共同癌种顺序按 Gene A reference continuous CN 的 median 升序（并列按癌种名）；`TCGA_Cancer_Order.csv` 是**视觉从上到下**的顺序，coord_flip 的 factor levels 为该表倒序。Provenance 保存实际绘图 factor levels、点数和转换规则，自动验证与表一致。

03 精确按 SampleID 一对一匹配 continuous CN、thresholded CNA、expression，限定 tumor，使用三者都 finite 的样本。逐癌种报告 N、Pearson/Spearman 和各自跨癌种 BH FDR；固定 N<20 或常量变量时不编造统计量，仍生成该癌种图。百分比 inset 的分母是同图匹配 cohort。癌种样本数与 CNA 频率分母不同的原因保留在样本表中。

## DepMap 统计定义

使用本地 26Q1 WGS CN 和 Chronos Parquet：单基因读必要列；genome-wide 和 CN covariation 各读取一次必要完整矩阵。用户提供 portal exports 的完整性相对官方全 release 未验证，使用其全部提供的 Chronos 基因列。

保留原始数值核心与 `CN_log=log2(relative CN+1)`、low `<0.585`。Genome-wide `Delta_median=median(low)-median(nonlow)`：负值表示 low 组更依赖，正值表示更弱。Wilcoxon P 和全筛选 BH FDR 不更改；历史 Rank 保留，主排名 Eligible_Rank 排除 undefined FDR。火山图标记每个方向 Top10，并额外高亮提供的 Gene B；候选表每方向 Top20。

Targeted scatter 的横轴和角落统计是 **CN_relative**；为保留历史定义，统计表另保留 CN_log Pearson/Spearman，不把两种 Pearson 默认为相同。Waterfall 按 Chronos 从高到低（弱到强）排序，并列按 ModelID；CN-Normal 为黑色，low 红色。Boxplot 展示全部点、N、组 median、Delta median、Wilcoxon P 和星号。

Covariation 报告 positive Top20 与 negative Top20，排除 Gene A 自身；同向与反向相关均不解释为因果关系。

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

# 精简后的本地数据库

两个正式入口保持不变：`geneA_screen` 和 `geneA_geneB`。数据版本固定为 TCGA GDC DR46 / 46.0 与用户提供的 DepMap 26Q1；不联网、不读取 raw、不回退到 PanCanAtlas。

## 默认分析实际依赖

| 模块 | 本地输入 |
|---|---|
| TCGA 五级 CNA、relative CN landscape、33 癌种 CN–mRNA | 完整 GeneLevel CN、CN 基因注释、完整 STAR TPM、TPM 基因注释、持久化 baseline、样本 metadata、RNA selection |
| 全基因组 dependency | 完整 OmicsCNGeneWGS、完整 CRISPRGeneEffect |
| CN covariation | 完整 OmicsCNGeneWGS |
| Gene A→Gene B 主图、threshold sensitivity | OmicsCNGeneWGS、CRISPRGeneEffect |
| Lineage、lineage-adjusted dependency | 上述两矩阵及完整 Model |

精确保留清单见 `config/analysis_storage.json`。既有矩阵本来就使用 ZSTD，保持原文件字节、double 精度、全部基因、样本顺序、缺失值与 Gene ID/Symbol 注释；没有重编码或有损压缩。按基因列投影读取，screen/covariation 仍读取各自完整矩阵。

表达量在 `get_tcga_expression()` 中从所请求基因的 TPM 计算 `log2(TPM + 1)`，不再依赖单独的整套 log 矩阵。`tcga_ready()` 只检查精简后的 processed 文件。

`TCGA_Sample_Baselines.parquet` 用原始完整常染色体 CN 算法重新核对后持久化，含 SampleID、BaselineCN、BaselineMethod、AutosomalFiniteGenes、ModalGeneCount、BaselineTies；并列取最小众数，不是固定 CN=2，不是实测 ploidy。sidecar 记录算法版本、输入 SHA256、文件 SHA256、原始方法报告。普通分析只加载已验证 baseline，不重新估算。

`TCGA_Sample_Metadata.parquet` 完整保留原 Sample_Map 的行与字段，并补充 CancerType、TumorNormal；保留 CaseID、ProjectID、SampleType、AliquotID、FileID，不新做 patient-level 去重。原 Sample_Map 仅记录为历史来源，不是运行依赖。

`TCGA_RNA_Selection.parquet` 保存全部 RNA 文件的 SampleID、FileID、AliquotID、对应 CN_AliquotID、Matches_CN_Aliquot、Selected、MatchStatus、SelectionBasis。规则完全沿用 matching CN aliquot 优先、再按 FileID 字典序；不依据某个基因的表达量/缺失值选样本。普通分析核对文件与样本身份并复用选择；现有审计 CSV 的列和数值不变。

持久化小表每次加载核对自身 SHA256 和输入大小/mtime；完整矩阵 SHA256 在准备、真实回归、删除前重复核验。更换核心数据后必须独立准备，不能静默复用过期小表。

## 验证与删除门槛

`scripts/utils/optimize_analysis_data.py prepare` 验证核心源 SHA256、完整基因注释、33 癌种、全部 TPM/log 数值及 Arrow null/NaN 状态、全部 CN 样本 baseline；只生成小表，不删除数据。已存在小表/报告时拒绝覆盖。

`scripts/utils/validate_analysis_ready.py` 在 `.runtime/database_optimization/core_only` 建立隔离项目：仅提供核心 Parquet、必要来源 sidecar 和已安装的运行环境。没有 raw、冗余 processed 或下载入口。Python audit hook 拒绝网络和 raw 读取，R 网络入口被阻止。分别执行正式 screen 和 A+B；A+B 必须复用 A 的全部三模块。逐列比较现有 MTAP→WDR77 的 main/supplementary/audit CSV，要求行顺序、字段、身份和缺失状态一致；数值 `rtol=5e-13, atol=0`，不以绝对容差接受极小 P 值变零。完整 TPM/log 转换比较 `rtol=2e-15, atol=0`。结果验证另外独立复算 Pearson、Spearman、P、BH、分组、Eligible_Rank，并按源列核验样本与 baseline。

R synthetic fixtures 覆盖特殊 CN baseline/分类、RNA 选择、常量和缺失数据、CN-Low N=0/1/2 及任一组不足时 SKIPPED。Python fixtures 覆盖输入过期、文件篡改、精度、缺失值、任意基因读取、缺少可选数据和删除越界。

`plan` 创建逐文件 SHA256/MD5、大小、mtime、类别清单，归档必要 DR46 来源清单；`cleanup` 在首个删除前核验所有目标和核心文件、代码哈希、回归证据以及活跃项目进程。仅逐项删除指定 data 目录内文件，不使用目录通配递归删除。逐项落盘删除 journal；任何失败立即停止并保留剩余文件。

所有真实结果保持原位。`migrate_analysis_storage_cache.R` 只为输入和输出仍与旧缓存严格相符的模块生成旧 key→新 key 的精确映射，不重跑统计、不改结果文件。运行时仍核验每个输出 SHA256；参数、数据、代码改变产生其他 key 时不能复用。映射绑定回归报告 SHA256。

## 溯源与不再支持的数据层

`data/manifests/analysis_ready_validation.json` 保存核心 SHA256/覆盖/全量验证；`analysis_ready_regression.json` 保存真实回归、隔离测试及代码哈希；`data_cleanup_report.json` 保存实际删除量和体积；`data_cleanup_files.jsonl` 是逐项删除清单，`data_cleanup_deleted.jsonl` 是删除 journal。GDC/DepMap 原始来源清单、历史 MD5/SHA256、预处理版本和 QC 继续保留。历史 raw QC 只说明当时验证通过，不表示文件现在还在。

`scripts/download/`、`scripts/preprocess/`、`verify_gdc_raw.py`、`validate_gdc_current.py`、`validate_database.py` 等属于历史/独立数据维护，不能用于精简后的默认分析或声称 raw 仍在。

原始 counts、FPKM、FPKM-UQ、segments、mutation、PanCanAtlas/Xena 及非核心 DepMap 数据删除后，对应的可选研究入口会明确报缺少数据，不会伪造结果或从其他数据源补充。今后开展这些研究需要另行准备数据。完整 Chronos 和 CN 矩阵、Model 信息均保留，两个正式工作流的基因/样本覆盖不变。

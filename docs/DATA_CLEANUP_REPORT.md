# 数据库精简验收结果（2026-10-20）

实际删除已完成，状态 PASS。未重新下载数据，TCGA 固定 DR46 / 46.0，DepMap 固定现有 26Q1。矩阵数值、精度、缺失状态、基因/样本覆盖、匹配规则、统计方法、阈值及图形定义保持原样。用户在删除前明确批准了逐文件清单。

## 实际体积

| 指标 | GiB | bytes |
|---|---:|---:|
| 精简前项目文件合计 | 103.746781 | 111397257737 |
| 精简后项目文件合计 | 4.502188 | 4834187143 |
| 项目净减少 | 99.244593 | 106563070594 |
| 删除 raw：45,667 个文件 | 86.695042 | 93088092365 |
| 删除 processed：35 个文件 | 12.612642 | 13542721246 |
| 删除文件合计：45,702 个 | 99.307684 | 106630813611 |
| 实测 D 盘可用空间增加 | 99.340538 | 106666090496 |

GiB = 2^30 bytes。项目总量使用精简前后相同的逐文件 byte 合计方法，包含环境、Git 和溯源文件；测量时隔离测试副本已经移除，Git 提交尚未创建。新增审计清单/报告使项目净减少略小于删除文件总量。D 盘可用空间增量是独立实测，NTFS 分配单元及同期小型文件变化可使其与文件 byte 合计不同；本轮执行了实际删除，没有把数据移到其他目录。

## 保留的核心库

全部采用既有或新建的 Parquet ZSTD；既有矩阵未重编码，未降为 float32，未降低小数精度。

| 文件 | MiB |
|---|---:|
| TCGA_GeneLevel_CN.parquet | 221.545 |
| TCGA_GeneLevel_CN.genes.parquet | 1.257 |
| TCGA_STAR_TPM.parquet | 2557.927 |
| TCGA_STAR_TPM.genes.parquet | 0.667 |
| TCGA_Sample_Baselines.parquet | 0.308 |
| TCGA_Sample_Metadata.parquet | 3.749 |
| TCGA_RNA_Selection.parquet | 1.409 |
| OmicsCNGeneWGS.parquet | 208.305 |
| CRISPRGeneEffect.parquet | 221.891 |
| Model.parquet | 0.133 |

TCGA：33 癌种，CN 60,623 基因 × 11,339 样本；TPM 60,660 基因 × 11,505 RNA 文件（11,370 SampleID）。CN/TPM 基因注释完整保留。Baseline 11,339 行已从完整 57,633 常染色体基因复核后持久化；metadata 44,673 行；RNA selection 保留全部 RNA 文件及代表选择标识。无新 patient-level 去重。

DepMap：CN 18,613 基因 × 1,118 模型；Chronos 18,531 基因 × 1,208 模型；Model 2,154 模型、49 字段。保留现有全部供给数据，并不新增对官方整套 release 完整性的断言。

## 已删除类别

| 类别 | 文件数 | GiB |
|---|---:|---:|
| DepMap 非核心 processed | 10 | 0.501316 |
| TCGA 其他冗余 processed 小表 | 4 | 0.027087 |
| TCGA_CN_Segments | 1 | 0.087926 |
| TCGA_Masked_Somatic_Mutation | 1 | 1.132931 |
| TCGA_STAR_FPKM | 2 | 2.354960 |
| TCGA_STAR_FPKMUQ | 2 | 2.365774 |
| TCGA_STAR_log2TPMplus1 | 2 | 3.254415 |
| TCGA_STAR_UnstrandedCounts | 2 | 0.867571 |
| TCGA 历史 reference processed | 11 | 2.020662 |
| DepMap 原始 CSV | 13 | 1.275087 |
| TCGA 下载临时文件 | 38 | 0.347759 |
| TCGA 原始 clinical/API 快照 | 48 | 0.602120 |
| TCGA 原始 Gene-Level CN | 11390 | 36.271648 |
| TCGA 原始 segments | 11189 | 0.067080 |
| TCGA 原始 RNA | 12138 | 45.433446 |
| TCGA 原始清单（必要溯源已保留） | 20 | 0.038055 |
| TCGA 原始 API 检索页 | 185 | 0.449458 |
| TCGA 原始 mutation | 10641 | 0.764360 |
| TCGA 历史 reference raw | 5 | 1.446029 |

原始清单、UUID、MD5/SHA256、来源版本及预处理记录继续保留。历史 raw 完整性记录只证明当时验证情况，不表示已删除文件现在仍存在。counts、segments、mutation、reference 或可选 DepMap expression 等数据已不可用；涉及这些层的未来研究需单独准备，默认工作流不会自动下载、换来源或伪造结果。

## 验证结论

- 核心十个 Parquet 删除前后 SHA256 不变；全部基因、样本、ID、double 精度和缺失状态保留。
- 全部 60,660 个 TPM 基因与旧 log2(TPM+1) 矩阵逐块核对，最大绝对差异 0。
- 全部 11,339 个 baseline 与既有真实审计一致；RNA/CN 代表匹配和 CaseID/ProjectID 一致。
- 仅提供核心文件的安全隔离项目完成 geneA_screen 和 geneA_geneB；raw 和待删除 processed 在环境中均不存在，Python/R 网络入口受到阻止，未调用数据维护/下载。
- MTAP→WDR77 的 21 份 main/supplementary/audit CSV 全部 SHA256 相同，字节级一致，最大数值差异 0，包括 Pearson/Spearman/P/FDR、CNA、relative CN、全基因筛选、Eligible_Rank、targeted、lineage/adjusted/sensitivity 和 CN covariation。
- A+B 验证复用 A 的三个模块；原始全基因组统计只运行一次。7 个有效 MTAP 模块缓存的精确新旧 key 映射通过独立复核；9 个原本失效的其他案例旧参数/布局缓存保持原样，未标记为有效。
- 删除前 R 4.5.0 九组 synthetic suites 和 Python 22 项测试通过。覆盖 CN-Low N=0/1/2 或任一组不足时 SKIPPED、缺失/常量、baseline 并列、RNA 选择、数据/输出篡改、删除越界、微小 P 值变零拒绝和活跃进程保护。
- 43 页 PDF 经 Poppler 渲染并完成视觉检查；固定 33 癌种 PDF 和 COAD/READ 独立覆盖通过。
- 删除后重新读取 MTAP TCGA、任意基因 TP53 表达及 858 个 MTAP→WDR77 匹配模型；完整 DepMap CN/Chronos/Model 矩阵读取通过；原 21 份结果 CSV 保持不变。跨语言源值复核以 round-trip CSV 解析和 rtol=5e-13、atol=0 处理序列化/机器精度，未修改分析值。
- 项目 raw 目录已清空移除，processed 只剩上述十个核心 Parquet。隔离测试矩阵别名已移除；代码、结果、配置、运行库和必要溯源文件保留。

## 可复核记录

[核心验证](../data/manifests/analysis_ready_validation.json) · [真实回归](../data/manifests/analysis_ready_regression.json) · [删除后验证](../data/manifests/analysis_ready_post_cleanup.json) · [精简报告 JSON](../data/manifests/data_cleanup_report.json) · [精确删除清单](../data/manifests/data_cleanup_files.jsonl) · [逐项删除记录](../data/manifests/data_cleanup_deleted.jsonl) · [缓存迁移元数据](../data/manifests/analysis_storage_cache_migration.json)

本地工作区开始时已有 96 个旧结果文件的未提交删除；本任务保留这些用户改动，不纳入本次 Git 提交。GitHub commit、push、CI 和远程 SHA 一致性见最终交付消息。

后续 CI 补充验收：Windows staging guard 的短/长路径根目录规范化已修复；新增路径别名、超限拒绝及越界拒绝测试后，Python 共 23 项通过。此修改仅影响 Git staging 安全检查，分析读取程序与统计代码未变，真实回归报告保留为当时执行的历史证据。

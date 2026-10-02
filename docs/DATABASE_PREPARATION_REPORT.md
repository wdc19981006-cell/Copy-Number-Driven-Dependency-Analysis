# 数据库准备报告

> 历史报告：仅描述 2026-10-02 16:28 的初始官方下载阶段。此后已接收并整理 13 类 DepMap 本地导出；当前状态见 [DATA_SOURCES.md](DATA_SOURCES.md) 和 [TCGA_CURRENT_REPORT.md](TCGA_CURRENT_REPORT.md)。

报告时间：2026-10-02T16:28:17+08:00。

最新 DepMap Public release：**26Q1**（官方目录自动识别）。按用户后续指示跳过 DepMap 实体下载；TCGA 下载、预处理及查询基础设施已执行。

## 下载成功文件

| 文件 | 字节数 | 大小 | SHA256 |
|---|---:|---:|---|
| `Gistic2_CopyNumber_Gistic2_all_data_by_genes.tsv.gz` | 468,113,188 | 446.43 MiB | `9f98810361c9d0e262c79da256ab5b663b0530f890348200b9ce0ed40dd9b1b6` |
| `Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.tsv.gz` | 84,674,099 | 80.75 MiB | `818e3a5ed5a9c0b8e7d7da73cacba48d79a75e0fa515694e381f4ae7e58ae6ce` |
| `TCGA_phenotype_denseDataOnlyDownload.tsv.gz` | 61,165 | 59.73 KiB | `ac3923bb155ada3391b653a811f6c8c9c4734593b0c2b0a1555728928d83e5a2` |
| `mc3.v0.2.8.PUBLIC.xena.tsv.gz` | 64,540,936 | 61.55 MiB | `50f6301defc8468fb4bdb7407236d75545edaade3bb73778da3c4bb954485cc6` |
| `tcga_RSEM_Hugo_norm_count.tsv.gz` | 935,272,500 | 891.95 MiB | `fdf320a8fb80ab9dc9b9ca230a3221b66be4023efde10fa6a71aea1569b4c023` |

## 磁盘空间

- DepMap raw：**0 B**。
- TCGA raw：**1.446 GiB**。
- processed 合计：**2.021 GiB**（Parquet、索引、sample order、provenance 和仅含 views 的 DuckDB）。
- `.runtime/` 软件依赖及文档不计入上述数据库空间。

## 未下载文件

| 数据库 | 文件 | 状态 / 原因 |
|---|---|---|
| DepMap | `CRISPRGeneDependency.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `CRISPRGeneEffect.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `CRISPRScreenMap.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `Model.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `ModelCondition.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `OmicsCNGene.csv` | 26Q1 官方目录未提供，不从旧 release 替代 |
| DepMap | `OmicsCNGeneWGS.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `OmicsExpressionTPMLogp1HumanAllGenesStranded.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `OmicsFusionFilteredSupplementary.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `OmicsGlobalSignatures.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `OmicsInferredMolecularSubtypes.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `OmicsProfiles.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `OmicsSomaticMutations.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `OmicsSomaticMutationsMatrixDamaging.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `OmicsSomaticMutationsMatrixHotspot.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| DepMap | `ScreenSequenceMap.csv` | 按用户要求跳过；官方 URL 需浏览器验证 |
| TCGA | `mc3.v0.2.8.PUBLIC.maf.gz` | 下载失败：官方 GDC 历史 MAF UUID 已失效，API 无此文件 |

精确失败信息见 manifest 的 error 列；单独清单见 `data/manifests/not_downloaded.csv`。TCGA 四个核心文件和 Xena MC3 补充表均成功；原始 MAF 可选文件失败 1 个。

## 完整性与预处理 QC

| 数据集 | 基因 | 样本 | 全部行解析 | 重复基因 | 缺失数值 |
|---|---:|---:|---|---:|---:|
| CN | 24,776 | 10,845 | True | 0 | 0 |
| GISTIC | 24,776 | 10,845 | True | 0 | 0 |
| Expression | 58,581 | 10,535 | True | 0 | 0 |

三个矩阵为 float64、ZSTD Parquet，一基因一 row group。所有列宽、TCGA ID、数值解析均检查；GISTIC 只接受 -2/-1/0/1/2。原始 gzip 文件完整读取并验证 SHA256。

### 样本与 phenotype overlap

| 数据集 | 矩阵样本 | phenotype 原表 overlap | 癌种缺失 |
|---|---:|---:|---:|
| CN | 10,845 | 10,845 | 0 |
| GISTIC | 10,845 | 10,845 | 0 |
| Expression | 10,535 | 10,534 | 1 |

统一 metadata 共 12,805 个样本，覆盖 33 个癌种。表达表的唯一 phenotype 缺失是 `TCGA-07-0249-20` 控制分析物，保留 Other、癌种 null；所有 Tumor/Normal 样本有癌种标注。患者 ID、样本类型、Tumor/Normal 已建立。

### DepMap ModelID overlap

**未测量**。Model.csv、CRISPR、CN、expression 均按用户指示未下载；不能报告 0 或声称通过。预处理脚本将来会检查 ACH-xxxxxx 格式、唯一映射及与 Model.csv 的 overlap，不要求各矩阵模型数量相同。

### MC3 补充表

完整 Xena 表已转换成 MC3.parquet：2,907,335 条记录，9,104 个样本，21,273 个基因；全部行解析且样本 ID 合法。

### 查询与安全验证

真实数据库验证通过：**True**；23 项 checks，0 项 errors。每个 TCGA 矩阵抽取 3 个基因，逐样本与原始值精确一致（float64），并检查 metadata。测试只验证数据读取，不做 Gene A/Gene B 分析。

5 个合成安全测试覆盖 release 数值排序、默认 profile 映射、歧义失败、gene column projection、断点续传、checksum skip 和拒绝覆盖异常 raw；全部通过。合成测试不等于真实 DepMap QC。

单基因查询使用显式 RowGroup 索引，仅读取一个基因向量；DuckDB views 另外保留供 SQL 使用。首次查询包括 footer/索引加载，后续查询只保留 footer 和基因索引缓存。

单基因查询耗时：

- CN / `ACAP3`：0.136 秒，10,845 样本。
- CN / `NUTM2D`：0.007 秒，10,845 样本。
- CN / `WASIR1|ENSG00000185203.7`：0.009 秒，10,845 样本。
- GISTIC / `ACAP3`：0.165 秒，10,845 样本。
- GISTIC / `NUTM2D`：0.006 秒，10,845 样本。
- GISTIC / `WASIR1|ENSG00000185203.7`：0.007 秒，10,845 样本。
- Expression / `CTD-2588J6.1`：0.337 秒，10,535 样本。
- Expression / `RP11-453O5.1`：0.006 秒，10,535 样本。
- Expression / `RP11-526P5.2`：0.006 秒，10,535 样本。

## 项目目录树

```text
Copy Number–Driven Dependency Analysis/
├─ data/
│  ├─ raw/
│  │  ├─ depmap/26Q1/             [空：本次跳过]
│  │  └─ tcga/pancan/             [5 个原始 gzip 文件]
│  ├─ processed/
│  │  ├─ depmap/26Q1/             [空：尚无 raw]
│  │  └─ tcga/                    [CN/GISTIC/Expression/MC3/Metadata + 索引 + DuckDB]
│  └─ manifests/                  [SHA256、缺失清单、源字典、QC]
├─ scripts/
│  ├─ download/
│  ├─ preprocess/
│  ├─ utils/
│  └─ data_access.py
├─ config/
├─ results/                       [无生物学分析]
├─ logs/
├─ docs/
├─ tests/
└─ README.md
```

## 后续更新

DepMap 定期检查新 release，并在明确选择后整套更新；旧版永久保留。TCGA 不需要定期重新下载。MC3 原始 MAF 仅在官方链接恢复或找到官方稳定版本后补充。DepMap 将来下载时需要含真实 URL 的同 release 官方目录。

自动审批曾拒绝排队的 DepMap 下载命令，理由是用户已要求跳过；该命令未执行。按后续指示完成其余已授权的 TCGA 工作。

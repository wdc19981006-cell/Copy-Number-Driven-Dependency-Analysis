# VPS4B → VPS4A

问题：VPS4B CN 降低是否伴随 VPS4A dependency 增加？这是观察性关联分析。

```bat
Rscript scripts/R/run_analysis.R --workflow geneA_geneB --geneA VPS4B --geneB VPS4A
Rscript modules/VPS4B_VPS4A/run_case.R
```

结果在 [中文结果指南](../../results/VPS4B_VPS4A_Analysis/00_Analysis_Summary.txt) 和 `results/VPS4B_VPS4A_Analysis/Main_Results`。40 张主 PDF 包含 33 癌种自身 CN–mRNA；targeted scatter、waterfall、分组 boxplot 是三个独立文件。补充结果只含 lineage、adjusted、CN threshold sensitivity。旧开发目录已在验收通过后移到本地忽略的 .runtime 归档。

本地 processed TCGA reference 和 DepMap 26Q1；NCI GDC DR46 保留为独立数据层。普通 workflow 不联网或扫描 raw。原始统计核心未改，所有 supplied Chronos 基因列进入筛选；full release 完整性未验证。

VPS4A 在 18,256 个 eligible tests 中 Eligible_Rank=1；matched CN-low/nonlow N=76/782。Delta median=-0.5607422623，Wilcoxon P=1.1533764853e-22，FDR=2.1056041116e-18。实际精确数值见 Tables，其他统计和方向解释见中文结果指南。

[最终工作流规范](../../docs/FINAL_WORKFLOWS.md) 说明本地数据边界、reference 样本匹配、癌种排序、缓存、独立验证和 GitHub 上传要求。反向研究应交换 A/B 重跑同一个 workflow。专项 `--mode` 仍可显式调用，普通用户不用 `full`。

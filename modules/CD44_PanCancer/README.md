# CD44 泛癌分析

本案例回答 CD44 的泛癌 CNA 频率、五级 CNA 与 RNA 的关系，以及低 CN 细胞模型的全基因 CRISPR 依赖。

结果：`results/CD44_PanCancer/Summary/CD44_泛癌分析报告.md`。输入为现有 TCGA GDC DR46、独立 PanCanAtlas/Xena reference，以及用户提供的 DepMap 26Q1 portal 导出；不宣称是截至今天最新的完整 DepMap release。

沿用项目阈值 `log2(relative CN+1)<0.585`。6 个模型符合阈值，其中 4 个匹配 CRISPR。18,531 个提供的 Chronos 基因进入筛选，17,787 个有合格 Wilcoxon 比较。AOC1、VWA5B2、PSMD10 有二元癌种校正回归信号，但没有基因同时通过 Wilcoxon 和校正回归两套全基因 FDR<0.05 标准。候选保留为探索线索。

参考层五级 CNA 频率分母 10,845，匹配 RNA 为 9,492；RNA 为 log2(norm_count+1)。当前层连续 CN–RNA 为 10,549 个肿瘤样本，RNA 为 log2(TPM+1)；按患者去重敏感性有 10,423 人。禁止混合两层单位。

在项目根目录的 Git Bash 中按顺序运行：

```bash
/c/Python312/python.exe "D:/CodexProjects/Copy Number–Driven Dependency Analysis/scripts/utils/run_r_case.py" --geneA CD44 --geneB CD44 --output_case CD44_PanCancer --modes qc depmap_cn_expression genomewide_dependency genomewide_adjusted_dependency tcga_cn_landscape tcga_cna_prevalence tcga_cn_expression cn_covariation
/c/Python312/python.exe "D:/CodexProjects/Copy Number–Driven Dependency Analysis/scripts/utils/prepare_cd44_supplement.py"
/d/R/R-4.5.0/bin/Rscript.exe --vanilla "D:/CodexProjects/Copy Number–Driven Dependency Analysis/scripts/R/cd44_supplement.R"
/d/R/R-4.5.0/bin/Rscript.exe --vanilla "D:/CodexProjects/Copy Number–Driven Dependency Analysis/tests/validate_cd44_results.R"
/c/Python312/python.exe "D:/CodexProjects/Copy Number–Driven Dependency Analysis/scripts/utils/report_cd44.py"
```

`geneB=CD44` 仅用于框架的默认候选位置；本任务从完整 genome-wide 结果选取其他依赖候选。`full` 还包含与本问题无关的定向／反向／突变模块，所以这里显式指定需要的模式。

专项脚本保存样本来源、所有五级分母、癌种内状态比较、癌种调整、逐一去除低 CN 模型、候选额外协变量和不能估计的原因。`tests/validate_cd44_results.R` 直接从 Parquet 重算代表候选并验证四套 BH 家族、RNA 对照和样本匹配。原统计核心未修改；自定义 output_case 排名导出路径修复由 `tests/test_output_case_finalize.R` 覆盖。

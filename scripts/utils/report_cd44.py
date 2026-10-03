"""Build a Chinese CD44 report directly from validated result tables."""
from pathlib import Path
import sys, json, math
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / '.runtime'))
import pandas as pd
base = ROOT / 'results/CD44_PanCancer'
out = base / '13_CD44_Supplement'
read = lambda n: pd.read_csv(out / (n + '.csv'), low_memory=False)
assert json.loads((out/'Independent_Validation.json').read_text())['all_passed']
def md(df, columns):
    rows = ['| ' + ' | '.join(label for _, label, _ in columns) + ' |',
            '| ' + ' | '.join('---' for _ in columns) + ' |']
    for _, row in df.iterrows():
        vals = []
        for col, _, fmt in columns:
            val = row[col]
            vals.append('NA' if pd.isna(val) else fmt(val))
        rows.append('| ' + ' | '.join(vals) + ' |')
    return '\n'.join(rows)
f = lambda x: f'{float(x):.3f}'
p = lambda x: f'{float(x):.3g}'
i = lambda x: str(int(x))
s = lambda x: str(x)
states = ['Deep deletion','Shallow deletion','Diploid','Gain','Amplification']
zh = dict(zip(states,['深缺失','浅缺失','中性／二倍体分类','低水平增益','高水平扩增']))
prev = read('Reference_CNA_Prevalence_Overall').set_index('CNA').loc[states].reset_index()
prev['中文分类'] = prev.CNA.map(zh)
freq = read('Reference_CNA_Prevalence_ByCancer')
groups = read('Reference_CNA_RNA_Overall_Groups').set_index('CNA').loc[states].reset_index()
groups['中文分类'] = groups.CNA.map(zh)
contrasts = read('Reference_CNA_RNA_CancerAdjusted')
contrasts = contrasts[contrasts.term.str.startswith('relevel')].copy()
contrasts['中文分类'] = [zh[next(st for st in states if st in t)] for t in contrasts.term]
current = read('Current_Tumor_CN_RNA_Overall').iloc[0]
by = read('Current_Tumor_CN_RNA_ByCancer')
positive = by[(by.Pearson_r>0)&(by.Pearson_FDR<.05)].sort_values('Pearson_r',ascending=False)
negative = by[(by.Pearson_r<0)&(by.Pearson_FDR<.05)]
negative_text = '；'.join(f'{r.CancerType} r={r.Pearson_r:.3f}、FDR={r.Pearson_FDR:.3g}' for r in negative.itertuples()) or '无'
coef = read('Current_Tumor_CN_RNA_Adjusted').set_index('term').loc['CN']
zcoef = read('Current_Tumor_CN_RNA_Adjusted_Standardized').set_index('term').loc['CN']
screen = read('CD44_AllGenes_CombinedScreen')
candidate_names = ['AOC1','VWA5B2','PSMD10','SNAP23','FOXE1']
candidates = screen.set_index('Gene').loc[candidate_names].reset_index()
loo = read('Candidate_LeaveOneLowModelOut')
cells = read('Candidate_Dependency_CellLines')
lowcells = cells[cells.CNLow].drop_duplicates('ModelID')[['ModelID','CellLineName','OncotreeLineage','CN_relative']]
cont = screen[(screen.Beta_CN>0)&(screen.FDR_CN<.05)].sort_values('FDR_CN')
stat = json.loads((out/'Supplement_Validation.json').read_text())
fig = lambda name: f'![{name}]({(out/(name+".png")).as_posix()})'
table_prev = md(prev,[('中文分类','CD44 拷贝数状态',s),('N','肿瘤样本数',i),('prevalence','比例',lambda x:f'{100*x:.2f}%')])
table_groups = md(groups,[('中文分类','状态',s),('N','匹配 RNA 样本数',i),('RNA_median','RNA 中位数',f)])
table_contrasts = md(contrasts,[('中文分类','相对中性分类',s),('estimate','癌种校正后 RNA 差值',f),('P_CNA_BH','四个状态比较的 BH-FDR',p)])
table_cand = md(candidates,[('Gene','基因',s),('N_low','低拷贝数 N',i),('N_nonlow','非低 N',i),('Delta_median','Chronos 中位数差（低−非低）',f),('Wilcoxon_FDR','全基因 Wilcoxon FDR',p),('Beta_CNLow','癌种校正低拷贝数系数',f),('FDR_CNLow','全基因校正回归 FDR',p)])
table_cancers = md(positive,[('CancerType','癌种',s),('N','N',i),('Pearson_r','Pearson r',f),('Pearson_FDR','BH-FDR',p)])
table_cells = md(lowcells,[('CellLineName','细胞系',s),('ModelID','ModelID',s),('OncotreeLineage','组织来源',s),('CN_relative','CD44 相对拷贝数',f)])
stable = []
for gene in candidate_names:
    a = loo[loo.Gene==gene]
    if len(a): stable.append(f'- {gene}：逐一移除低拷贝数模型后，{int((a.Delta_median<0).sum())}/{len(a)} 次中位数差仍为负，{int((a.AdjustedBeta<0).sum())}/{len(a)} 次校正系数仍为负。')
gains = freq.sort_values('GainAmplification_pct',ascending=False).head(5)
gain_text = '、'.join(f'{r.CancerType} {r.GainAmplification_pct:.1f}%' for r in gains.itertuples())
lines = f'''# CD44 泛癌拷贝数、RNA 表达及基因依赖分析

分析日期：2026-10-03（Asia/Shanghai）。数据库：本地已准备的 TCGA GDC DR46、独立 TCGA PanCanAtlas/Xena 参考层，以及 DepMap 26Q1 用户导出。

## 主要结论

CD44 的五级 CNA 中，浅缺失占 18.05%，深缺失占 0.23%，增益占 10.88%，扩增占 1.39%。在参考层中，缺失伴随较低 RNA，增益／扩增伴随较高 RNA；癌种校正后四种状态与中性分类的差异均显著。当前 GDC 肿瘤样本中的连续 CN–RNA 正相关很弱（Pearson r={current.Pearson_r:.3f}），表明拷贝数仅解释 CD44 表达差异的一小部分。

**CD44 低拷贝数的稳定依赖基因尚不能确定。** 主分析的严格低拷贝数组只有 4 个可匹配 CRISPR 的细胞模型；没有候选同时通过全基因 Wilcoxon FDR 和癌种校正回归 FDR。AOC1、VWA5B2、PSMD10 可以作为探索线索，不能据此写成已证实的 CD44 缺失合成致死靶点。

## 泛癌拷贝数变化

分母为有 CD44 GISTIC 值的 10,845 个肿瘤样本，覆盖 33 癌种。频率按样本计算，不把同一患者不同样本平均成一个观测值。GISTIC −2／−1／0／+1／+2 是推定的深缺失、浅缺失、中性、增益、扩增分类，不是绝对的 0／1／2／3／4 个拷贝。[cBioPortal 官方分类说明](https://docs.cbioportal.org/user-guide/faq/)

{table_prev}

浅缺失与深缺失合计 18.28%；增益与扩增合计 12.27%。缺失比例前列为 TGCT 66.0%、UCS 41.1%、BLCA 40.7%、PCPG 34.0%、OV 31.6%。这些主要是浅缺失，不能等同于 CD44 完全丢失。增益／扩增合计比例前列为 {gain_text}。

{fig('CD44_CNA_Prevalence')}

## 各种拷贝数状态与 RNA 表达

五级 CNA 比较使用同一参考层、同一样本 ID 的 9,492 个肿瘤样本。RNA 单位为 **log2(norm_count+1)**，不与当前 GDC 的 TPM 值混用。

{table_groups}

上表为合并癌种的描述统计。回归模型 `RNA ~ CNA状态 + CancerType` 以中性分类为参照，得到以下癌种校正差异；BH 校正家族为四个状态系数。

{table_contrasts}

差值属于上述 log RNA 尺度，不应直接称为 TPM 倍数变化。深缺失匹配组只有 22 个样本，癌种内分析常因不足 3 个样本而不做检验。全部癌种内状态比较与各组 N 已保存，不能将合并癌种的结果解释为每个癌种都相同。

{fig('CD44_CNA_RNA_Overall')}

分癌种五级表达箱线图见 [CD44_CNA_RNA_ByCancer.png]({(out/'CD44_CNA_RNA_ByCancer.png').as_posix()})；癌种内状态与中性分类的检验覆盖所有可检验的癌种×状态比较，并在整个家族中统一做 BH 校正。

## 当前 GDC 层的连续 CN–RNA 验证

RNA 使用 STAR log2(TPM+1)，CN 为源文件提供的总基因拷贝数。按同一样本 UUID 匹配，一份 RNA 代表样本；优先匹配 CN aliquot，其次按文件 UUID 排序。样本类型审计证实 10,549 个匹配样本均属肿瘤／癌症，无正常样本。

- 泛癌 Pearson r={current.Pearson_r:.4f}，P={current.Pearson_P:.3g}；Spearman rho={current.Spearman_rho:.4f}，P={current.Spearman_P:.3g}。Pearson r² 约 {100*current.Pearson_r**2:.2f}%，这是未调整模型的描述性方差解释比例。
- 癌种校正的 CN 系数={coef.estimate:.4f} log2(TPM+1)／一个源 CN 单位，P={coef['p.value']:.3g}；全样本标准化系数={zcoef.estimate:.4f}。
- 每位患者只保留一个肿瘤样本的敏感性分析有 10,423 人，Pearson r=0.0707；癌种校正系数=0.1105，结论方向不变。
- {len(positive)} 个癌种有 Pearson 正相关且跨癌种 BH-FDR<0.05，详见下表。显著负相关癌种为：{negative_text}。其余癌种不应默认存在同样关联。

{table_cancers}

{fig('CD44_Current_Tumor_CN_RNA_ByCancer')}

两套 TCGA 分析分别属于参考 GISTIC／归一化计数层和当前 GDC／TPM 层；它们的样本覆盖和单位不同，不能直接比较效应值大小。

## CD44 低拷贝数更依赖哪些基因

沿用项目已定义的阈值：`CN_log=log2(relative CN+1)`，`CN_log<0.585` 为低拷贝数，近似要求相对 CN<0.5。它是较严格的相对低拷贝数定义，不是 TCGA GISTIC −1／−2。相对 CN 对应细胞模型的整体倍性，不能直接当作绝对拷贝数。[DepMap 官方相对 CN 说明](https://forum.depmap.org/t/what-is-relative-copy-number-copy-number-ratio/104)

1,118 个 CN 模型中只有 6 个符合阈值；其中 4 个有匹配 CRISPR 数据，主筛选通常为 4 vs 854 个模型。每个靶基因还按自身的缺失值筛选，实际 N 见完整表。低拷贝数可匹配模型如下：

{table_cells}

任何单一组织来源的低拷贝数组都不足 3 个，因此本数据无法完成可靠的癌种内低组验证。Chronos 越负表示敲除越影响生长；下面的负中位数差表示低 CN 组相对更依赖。[Chronos 原始论文](https://link.springer.com/article/10.1186/s13059-021-02540-7)

18,531 个提供的 Chronos 基因被筛查；其中 17,787 个具备两组 Wilcoxon 检验条件，18,435 个具备癌种校正连续回归条件。依赖模型使用 `Chronos ~ CD44_CN_log + OncotreeLineage` 和 `Chronos ~ CD44_CNLow + OncotreeLineage`，各自的全基因 P 值独立 BH 校正。

{table_cand}

**没有同时满足两类全基因 FDR<0.05 的候选。** AOC1、VWA5B2、PSMD10 的校正二元回归有信号，但全基因 Wilcoxon FDR 均约 0.648；不能只报告更小的那一项 FDR。SNAP23 是原始 Wilcoxon Eligible_Rank 前列，原因是大量基因 FDR 并列时按中位数差排序；它不等于获得最强统计支持。

逐一移除低 CN 模型的描述性稳定性检查如下，P 值属于筛选后的敏感性，未作为新的独立发现：

{chr(10).join(stable)}

FOXE1 在移除一个低 CN 模型后方向会改变，不宜优先作为稳定候选。所有候选的自身 CN、CD44 RNA 额外协变量及排除 CN 最高 1% 的连续回归结果均保留，不足样本／不可估计时明确标记跳过。

连续癌种校正筛选有 {len(cont)} 个正 CN 系数且 FDR<0.05 的基因：{', '.join(cont.Gene)}。正系数表示随 CD44 CN 降低，Chronos 倾向更负；但这些基因没有通过严格低 CN 分组的两类联合标准。SMIM11、OR1L3、RAET1G 的效应常接近 0，不能据相关性直接称为强依赖；RPL6 覆盖只有 244 个模型、低组 2 个，未满足两组检验条件。

{fig('CD44_Dependency_Candidates')}

## 解释边界和复核

本任务检验关联，不证明 CD44 本身导致这些依赖。癌种已校正，肿瘤纯度、免疫／间质细胞比例、大片段共缺失及 CRISPR 屏幕特征尚未全部校正。RNA 是基因整体表达，未区分 CD44 剪接异构体；患者肿瘤表达关联和细胞系功能依赖是两类证据。

DepMap 输入是用户提供的 26Q1 portal 导出；已分析所有提供的基因列，但尚未验证其覆盖等于官方完整 release。本文不宣称使用截至今天最新的全部公开 DepMap 数据。

独立复核通过：源 Parquet 样本对齐、六个代表候选的统计重算、四个全筛选 BH 家族、Eligible_Rank、参考层分母和所有可检验的癌种内表达比较、当前层样本身份及癌种校正回归。原始十个统计核心函数体一致性测试通过。图表已检查可读性。

本次还修复了自定义输出目录的排名导出路径；两种目录的回归测试通过，统计公式未修改。早期 `Module_Runs.csv` 保留了导出失败的原始记录；失败步骤已按实际目录补齐，后续模块和独立验证全部完成。调用框架时 `geneB=CD44` 仅提供默认候选位置，本报告的依赖结果来自全基因筛选。

## 结果文件与重现

结果目录：`results/CD44_PanCancer/`。核心交付：本报告、5 张专项 PNG、原框架 PDF／CSV、完整全基因筛选 CSV 和独立验证 JSON。

- [完整联合筛选表]({(out/'CD44_AllGenes_CombinedScreen.csv').as_posix()})
- [癌种内 CNA–RNA 比较]({(out/'Reference_CNA_RNA_ByCancer_vsDiploid.csv').as_posix()})
- [五级 CNA 频率]({(out/'Reference_CNA_Prevalence_ByCancer.csv').as_posix()})
- [独立验证记录]({(out/'Independent_Validation.json').as_posix()})

在项目根目录的 Git Bash 中依次执行（默认目录已修复）：

```bash
/c/Python312/python.exe "D:/CodexProjects/Copy Number–Driven Dependency Analysis/scripts/utils/run_r_case.py" --geneA CD44 --geneB CD44 --output_case CD44_PanCancer --modes qc depmap_cn_expression genomewide_dependency genomewide_adjusted_dependency tcga_cn_landscape tcga_cna_prevalence tcga_cn_expression cn_covariation
/c/Python312/python.exe "D:/CodexProjects/Copy Number–Driven Dependency Analysis/scripts/utils/prepare_cd44_supplement.py"
/d/R/R-4.5.0/bin/Rscript.exe --vanilla "D:/CodexProjects/Copy Number–Driven Dependency Analysis/scripts/R/cd44_supplement.R"
/d/R/R-4.5.0/bin/Rscript.exe --vanilla "D:/CodexProjects/Copy Number–Driven Dependency Analysis/tests/validate_cd44_results.R"
/c/Python312/python.exe "D:/CodexProjects/Copy Number–Driven Dependency Analysis/scripts/utils/report_cd44.py"
```
'''
report = base/'Summary/CD44_泛癌分析报告.md'
report.write_text(lines, encoding='utf-8')
print('REPORT:',report.as_posix())
print('Positive current cancers:', positive[['CancerType','Pearson_r','Pearson_FDR']].to_dict('records'))
print('Selected candidates:',candidates[['Gene','Delta_median','Wilcoxon_FDR','FDR_CNLow']].to_dict('records'))

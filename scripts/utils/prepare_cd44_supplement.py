"""Export matched reference-layer CD44 data for five-state RNA analysis."""
from pathlib import Path
import sys, json
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.data_access import get_tcga_gistic, get_tcga_expression, get_tcga_cn, get_tcga_metadata
from scripts.data_access import get_depmap_cn, get_depmap_expression, get_depmap_model_metadata
import numpy as np

out = ROOT / 'results/CD44_PanCancer/13_CD44_Supplement'
out.mkdir(parents=True, exist_ok=True)
meta = get_tcga_metadata(layer='reference')
ref = get_tcga_gistic('CD44').merge(get_tcga_expression('CD44', layer='reference'), on='SampleID', validate='one_to_one')
ref = ref.merge(get_tcga_cn('CD44', layer='reference'), on='SampleID', validate='one_to_one')
ref = ref.merge(meta[['SampleID','CancerType','SampleType','TumorNormal']], on='SampleID', validate='one_to_one')
ref = ref.loc[(ref.TumorNormal == 'Tumor') & ref.CancerType.notna() & np.isfinite(ref.GISTIC) & np.isfinite(ref.Expression)].copy()
assert not ref.SampleID.duplicated().any()
assert set(ref.GISTIC.unique()).issubset({-2,-1,0,1,2})
ref['DataLayer'] = 'PanCanAtlas_Xena_reference'
ref['ExpressionUnit'] = 'log2(norm_count+1); not TPM'
ref.to_csv(out / 'TCGA_Reference_CNA_RNA_Samples.csv', index=False)
dm = get_depmap_cn('CD44').rename(columns={'CopyNumber':'CN_relative'})
dm = dm.merge(get_depmap_expression('CD44'), on='ModelID', how='left', validate='one_to_one')
dm = dm.merge(get_depmap_model_metadata(), on='ModelID', how='left', validate='one_to_one')
dm['CN_log'] = np.log2(dm.CN_relative + 1)
dm['CNLow'] = dm.CN_log < .585
dm.loc[dm.CNLow].to_csv(out / 'CD44_CNLow_All_Models.csv', index=False)
info = {'reference_matched_tumor_samples':len(ref), 'reference_cancers':int(ref.CancerType.nunique()),
        'depmap_CN_models':len(dm), 'CNLow_models':int(dm.CNLow.sum()),
        'reference_expression_unit':'log2(norm_count+1)', 'depmap_low_definition':'log2(relative CN+1)<0.585; analysis-defined',
        'preserved_core':'scripts/R/run_analysis.R unchanged; supplement is independent'}
(out / 'Supplement_Metadata.json').write_text(json.dumps(info, indent=2), encoding='utf-8')
print(json.dumps(info), flush=True)

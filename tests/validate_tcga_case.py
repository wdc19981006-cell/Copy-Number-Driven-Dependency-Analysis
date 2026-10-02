"""Check real current/reference separation and sample selection in case outputs."""
from pathlib import Path
import sys,json,csv
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.utils.common import ROOT,now,write_json
from scripts.utils.gdc import RAW
import pandas as pd
import numpy as np

def main():
 folder=ROOT/'results/VPS4B_VPS4A/08_TCGA'
 cn=pd.read_csv(folder/'TCGA_CN_Landscape.csv')
 rna=pd.read_csv(folder/'VPS4B_expression_current.csv')
 pairs=pd.read_csv(folder/'TCGA_CN_Expression_Samples.csv')
 stats=pd.read_csv(folder/'TCGA_CN_Expression_Statistics.csv').iloc[0]
 audit=pd.read_csv(folder/'TCGA_RNA_Representative_Selection.csv')
 selected=json.loads((RAW/'metadata/cn_selected.json').read_text(encoding='utf-8'))
 assert len(cn)==len(selected)==11339 and not cn.SampleID.duplicated().any()
 for frame in [cn,rna,audit]:assert frame.GDCRelease.eq('DR46').all() and frame.DataLayer.eq('gdc_DR46').all()
 assert pairs.GDCRelease_CN.eq('DR46').all() and pairs.GDCRelease_RNA.eq('DR46').all()
 assert pairs.DataLayer_CN.eq('gdc_DR46').all() and pairs.DataLayer_RNA.eq('gdc_DR46').all()
 actual=cn.set_index('FileID')[['SampleID','Workflow']].sort_index()
 expected=pd.DataFrame(selected).set_index('file_id')[['sample_id','workflow']].sort_index()
 assert actual.index.tolist()==expected.index.tolist() and actual.to_numpy().tolist()==expected.to_numpy().tolist()
 rna['matches_CN_aliquot']=rna.AliquotID.eq(rna.SampleID.map(cn.set_index('SampleID').AliquotID)).fillna(False)
 representatives=rna.sort_values(['SampleID','matches_CN_aliquot','FileID'],ascending=[True,False,True]).drop_duplicates('SampleID')
 flags=audit.selected_for_sample_analysis.astype(str).str.lower().isin(['true','1'])
 assert audit.loc[flags,'FileID'].sort_values().tolist()==representatives.FileID.sort_values().tolist()
 joined=cn[['SampleID','Value']].merge(representatives[['SampleID','Value']],on='SampleID',validate='one_to_one',suffixes=('_CN','_RNA'))
 joined=joined[np.isfinite(joined.Value_CN)&np.isfinite(joined.Value_RNA)].sort_values('SampleID')
 pairs=pairs.sort_values('SampleID')
 assert not pairs.SampleID.duplicated().any() and pairs.SampleID.tolist()==joined.SampleID.tolist()
 # R fwrite serializes doubles to 15 significant digits; Python's bridge CSV
 # retains more digits. Check that only this documented serialization rounding
 # differs, while sample/file identity remains exact above.
 actual_values=pairs[['Value_CN','Value_RNA']].to_numpy()
 expected_values=joined[['Value_CN','Value_RNA']].to_numpy()
 np.testing.assert_allclose(actual_values,expected_values,rtol=1e-14,atol=1e-14)
 max_csv_rounding=float(np.max(np.abs(actual_values-expected_values)))
 pearson=float(np.corrcoef(pairs.Value_CN,pairs.Value_RNA)[0,1])
 spearman=float(np.corrcoef(pairs.Value_CN.rank(method='average'),pairs.Value_RNA.rank(method='average'))[0,1])
 assert abs(pearson-stats.Pearson_r)<1e-12 and abs(spearman-stats.Spearman_rho)<1e-12 and int(stats.N)==len(pairs)
 modes=pd.read_csv(ROOT/'results/VPS4B_VPS4A/Summary/Module_Runs.csv').drop_duplicates('mode',keep='last').set_index('mode')
 assert all(modes.loc[k,'status']=='completed' for k in ['tcga_cn_landscape','tcga_cn_expression','tcga_cna_prevalence'])
 reference=pd.read_csv(folder/'VPS4B_gistic_reference.csv')
 prevalence=pd.read_csv(folder/'TCGA_GISTIC_CNA_Prevalence.csv')
 assert 'GDCRelease' not in reference and len(prevalence)>0 and prevalence.N.sum()==len(reference[(reference.TumorNormal=='Tumor') & reference.Value.notna() & reference.CancerType.notna()])
 assert prevalence.groupby('CancerType').prevalence.sum().between(.999999999,1.000000001).all()
 for name in ['TCGA_CN_Landscape.pdf','TCGA_CN_Expression.pdf','TCGA_GISTIC_CNA_Prevalence.pdf']:
  assert (folder/name).read_bytes().startswith(b'%PDF') and (folder/name).stat().st_size>1000
 write_json(ROOT/'data/manifests/tcga_case_validation.json',{'checked_at':now(),'all_passed':True,'current_release':'DR46','current_CN_samples':len(cn),'current_CN_finite_samples':int(cn.Value.notna().sum()),'matched_N':len(pairs),'Pearson_r':pearson,'Spearman_rho':spearman,'max_absolute_CSV_rounding_difference':max_csv_rounding,'CSV_rounding_tolerance':{'rtol':1e-14,'atol':1e-14},'unique_sample_mapping':True,'workflow_selection_matches_official_audit':True,'RNA_representative_selection_validated':True,'current_reference_separation':True,'latest_TCGA_modes_completed':True})
 print('PASS: real DR46 case mapping, original workflow selection, unique RNA representatives, independent Pearson/Spearman, latest module status and reference separation.',len(cn),len(pairs),pearson,spearman)

if __name__=='__main__':main()

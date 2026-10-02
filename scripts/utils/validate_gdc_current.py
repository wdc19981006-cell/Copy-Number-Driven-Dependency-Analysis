"""Validate current-layer coverage and raw/derived values without loading matrices."""
from pathlib import Path
import json,csv,time,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,now,write_json
from scripts.utils.gdc import RAW,PROCESSED,manifest_path,runtime_path,local_path
from scripts.preprocess.preprocess_gdc import selections,rna_frame,cn_frame
import numpy as np
import pyarrow.parquet as pq
import duckdb

def main():
 checks=[]
 def check(name,ok,detail=None):
  checks.append({'name':name,'passed':bool(ok),'detail':detail})
  if not ok:raise ValueError(name+': '+str(detail))
 with manifest_path().open(encoding='utf-8',newline='') as f:records=list(csv.DictReader(f))
 expected={k:len(selections(k)) for k in ['rna','cn','segments','mutation']}
 for kind,n in expected.items():
  rows=[r for r in records if r['kind']==kind]
  check(kind+' complete publisher checksums',len(rows)==n and all(r['status']=='verified' and r['MD5']==r['publisher_md5'] and len(r['SHA256'])==64 for r in rows),{'expected':n,'verified':sum(r['status']=='verified' for r in rows)})
 names=['TCGA_STAR_UnstrandedCounts','TCGA_STAR_TPM','TCGA_STAR_FPKM','TCGA_STAR_FPKMUQ','TCGA_STAR_log2TPMplus1','TCGA_GeneLevel_CN','TCGA_CN_Segments','TCGA_Masked_Somatic_Mutation','TCGA_Clinical','TCGA_Biospecimen','TCGA_Sample_Map']
 sizes={}
 for name in names:
  file=PROCESSED/f'{name}.parquet';pf=pq.ParquetFile(file)
  sizes[name]={'rows':pf.metadata.num_rows,'columns':len(pf.schema_arrow),'bytes':file.stat().st_size}
  check(name+' nonempty',pf.metadata.num_rows>0)
  if name.startswith('TCGA_STAR_'):check(name+' source file coverage',pf.metadata.num_rows==expected['rna'])
 cn=pq.read_table(PROCESSED/'TCGA_GeneLevel_CN.parquet',columns=['SampleID','FileID']).to_pandas()
 check('CN unique sample representative',len(cn)==expected['cn'] and not cn.SampleID.duplicated().any())
 for kind,read,names_map in [('rna',rna_frame,{'unstranded':'TCGA_STAR_UnstrandedCounts','tpm_unstranded':'TCGA_STAR_TPM','fpkm_unstranded':'TCGA_STAR_FPKM','fpkm_uq_unstranded':'TCGA_STAR_FPKMUQ'}),('cn',cn_frame,{'copy_number':'TCGA_GeneLevel_CN'})]:
  row=selections(kind)[0];source=read(local_path(kind,row));source=source[0] if kind=='rna' else source
  indices=np.linspace(0,len(source)-1,5,dtype=int);subset=source.iloc[indices]
  for field,dataset in names_map.items():
   columns=['FileID',*subset.gene_id.tolist()]
   pf=pq.ParquetFile(PROCESSED/f'{dataset}.parquet');table=pf.read_row_group(0,columns=columns).to_pandas();derived=table.loc[table.FileID==row['file_id'],subset.gene_id.tolist()].iloc[0].to_numpy(dtype=float)
   original=subset[field].to_numpy(dtype=float)
   check(dataset+' original values',np.allclose(derived,original,rtol=0,atol=1e-12,equal_nan=True))
   if field=='tpm_unstranded':
    table=pq.ParquetFile(PROCESSED/'TCGA_STAR_log2TPMplus1.parquet').read_row_group(0,columns=columns).to_pandas()
    vector=table.loc[table.FileID==row['file_id'],subset.gene_id.tolist()].iloc[0].to_numpy(dtype=float)
    check('log2(TPM+1) exact transformation',np.allclose(vector,np.log2(original+1),rtol=0,atol=1e-12))
 with duckdb.connect(str(ROOT/'data/processed/tcga/tcga_current.duckdb'),read_only=True) as con:
  con.execute("SET memory_limit='2GB'")
  n=con.execute('SELECT count(*) FROM TCGA_Masked_Somatic_Mutation WHERE aliquot_id != Tumor_Sample_UUID OR aliquot_barcode != Tumor_Sample_Barcode').fetchone()[0]
  check('MAF source tumor identity preserved',n==0)
  n=con.execute('SELECT count(*) FROM TCGA_Masked_Somatic_Mutation WHERE VAF < 0 OR VAF > 1').fetchone()[0]
  check('VAF range',n==0)
  n=con.execute('SELECT count(*) FROM TCGA_Clinical').fetchone()[0]
  check('clinical 11428 cases',n==11428)
  n=con.execute('SELECT count(DISTINCT project_id) FROM TCGA_Clinical').fetchone()[0]
  check('clinical 33 projects',n==33)
 write_json(runtime_path('gdc_validation_report.json'),{'created_at':now(),'checks':checks,'all_passed':all(c['passed'] for c in checks),'processed':sizes})
 layers=json.loads(runtime_path('current_layer_status.json').read_text(encoding='utf-8'))
 layers.update(current_status='complete',validated_at=now())
 write_json(runtime_path('current_layer_status.json'),layers)
 print('PASS',len(checks),'current-layer checks')

if __name__=='__main__':main()

"""Validate current-layer coverage and raw/derived values without loading matrices."""
from pathlib import Path
import json,csv,time,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,now,write_json,hashes
from scripts.utils.gdc import RAW,PROCESSED,manifest_path,runtime_path,local_path
from scripts.preprocess.preprocess_gdc import selections,rna_frame,cn_frame
from scripts.utils.verify_gdc_raw import verify_raw
from scripts.download.discover_gdc import select_cn
import numpy as np
import pyarrow.parquet as pq
import duckdb

def main():
 import logging
 logging.basicConfig(level=logging.INFO,format='%(asctime)s %(message)s')
 layers=json.loads(runtime_path('current_layer_status.json').read_text(encoding='utf-8'))
 layers.update(current_status='validating')
 write_json(runtime_path('current_layer_status.json'),layers)
 checks=[]
 def check(name,ok,detail=None):
  checks.append({'name':name,'passed':bool(ok),'detail':detail})
  if not ok:raise ValueError(name+': '+str(detail))
 with manifest_path().open(encoding='utf-8',newline='') as f:records=list(csv.DictReader(f))
 expected={k:len(selections(k)) for k in ['rna','cn','segments','mutation']}
 for kind,n in expected.items():
  rows=[r for r in records if r['kind']==kind]
  check(kind+' complete publisher checksums',len(rows)==n and all(r['status']=='verified' and r['MD5']==r['publisher_md5'] and len(r['SHA256'])==64 for r in rows),{'expected':n,'verified':sum(r['status']=='verified' for r in rows)})
 integrity=verify_raw()
 check('All selected raw actual MD5/SHA256 rechecked',integrity['all_passed'],integrity['kinds'])
 names=['TCGA_STAR_UnstrandedCounts','TCGA_STAR_TPM','TCGA_STAR_FPKM','TCGA_STAR_FPKMUQ','TCGA_STAR_log2TPMplus1','TCGA_GeneLevel_CN','TCGA_CN_Segments','TCGA_Masked_Somatic_Mutation','TCGA_Clinical','TCGA_Biospecimen','TCGA_Sample_Map']
 sizes={}
 for name in names:
  file=PROCESSED/f'{name}.parquet';pf=pq.ParquetFile(file)
  sizes[name]={'rows':pf.metadata.num_rows,'columns':len(pf.schema_arrow),'bytes':file.stat().st_size}
  provenance=file.with_suffix('.provenance.json')
  if provenance.exists():
   recorded=json.loads(provenance.read_text(encoding='utf-8'))
   check(name+' derived SHA256',hashes(file)[0]==recorded['SHA256'])
  check(name+' nonempty',pf.metadata.num_rows>0)
  if name.startswith('TCGA_STAR_'):check(name+' source file coverage',pf.metadata.num_rows==expected['rna'])
  if name.startswith('TCGA_STAR_') or name=='TCGA_GeneLevel_CN':
   check(name+' DR46 schema provenance',pf.schema_arrow.metadata.get(b'release')==b'GDC DR46.0')
   gene_ids=pq.read_table(PROCESSED/f'{name}.genes.parquet',columns=['gene_id'])['gene_id'].to_pylist()
   sizes[name]['genes']=len(gene_ids)
   check(name+' complete gene annotation/columns',[c for c in pf.schema_arrow.names if c.startswith('ENSG')]==gene_ids)
 library=pq.ParquetFile(PROCESSED/'TCGA_STAR_Library_Summary.parquet')
 check('STAR N_* library summary coverage',library.metadata.num_rows==expected['rna'])
 cn=pq.read_table(PROCESSED/'TCGA_GeneLevel_CN.parquet',columns=['SampleID','FileID']).to_pandas()
 check('CN unique sample representative',len(cn)==expected['cn'] and not cn.SampleID.duplicated().any())
 with (RAW/'metadata/tcga_cn_all_metadata.tsv').open(encoding='utf-8',newline='') as source:
  original_selection=select_cn(list(csv.DictReader(source,delimiter='\t')))
 check('CN workflow priority and aliquot/file tie break',[(r['sample_id'],r['file_id'],r['workflow']) for r in selections('cn')]==[(r['sample_id'],r['file_id'],r['workflow']) for r in original_selection])
 identity_fields={'SampleID':'sample_id','SampleBarcode':'sample','AliquotID':'aliquot_id','AliquotBarcode':'aliquot','CaseID':'case_id','CaseBarcode':'case','ProjectID':'project_id','FileID':'file_id'}
 bio=pq.read_table(PROCESSED/'TCGA_Biospecimen.parquet',columns=['aliquot_id','sample_id','case_id','aliquot_submitter_id','sample_submitter_id']).to_pylist()
 biospecimen={r['aliquot_id']:r for r in bio if r['aliquot_id']}
 for kind,dataset in [('rna','TCGA_STAR_TPM'),('cn','TCGA_GeneLevel_CN')]:
  ids=pq.read_table(PROCESSED/f'{dataset}.parquet',columns=list(identity_fields)).to_pylist()
  selection=selections(kind)
  check(kind+' processed identity matches selected metadata',len(ids)==len(selection) and all(all(a[k]==b[v] for k,v in identity_fields.items()) for a,b in zip(ids,selection)))
  check(kind+' official biospecimen hierarchy identity',all((r['aliquot_id'] in biospecimen and biospecimen[r['aliquot_id']]['sample_id']==r['sample_id'] and biospecimen[r['aliquot_id']]['case_id']==r['case_id'] and biospecimen[r['aliquot_id']]['aliquot_submitter_id']==r['aliquot'] and biospecimen[r['aliquot_id']]['sample_submitter_id']==r['sample']) for r in selection))
 for kind,read,names_map in [('rna',rna_frame,{'unstranded':'TCGA_STAR_UnstrandedCounts','tpm_unstranded':'TCGA_STAR_TPM','fpkm_unstranded':'TCGA_STAR_FPKM','fpkm_uq_unstranded':'TCGA_STAR_FPKMUQ'}),('cn',cn_frame,{'copy_number':'TCGA_GeneLevel_CN'})]:
  selected=selections(kind)
  sample_indices=set(np.linspace(0,len(selected)-1,5,dtype=int).tolist())
  if kind=='cn':
   first_workflow={}
   for i,row in enumerate(selected):first_workflow.setdefault(row['workflow'],i)
   sample_indices.update(first_workflow.values())
  sources=[]
  for index in sorted(sample_indices):
   row=selected[index];source=read(local_path(kind,row));source=source[0] if kind=='rna' else source
   indices=set(np.linspace(0,len(source)-1,5,dtype=int).tolist())
   indices.update(np.flatnonzero(source.gene_name.to_numpy()=='VPS4B').tolist())
   sources.append((index,row,source.iloc[sorted(indices)]))
  for field,dataset in names_map.items():
   pf=pq.ParquetFile(PROCESSED/f'{dataset}.parquet')
   log_file=pq.ParquetFile(PROCESSED/'TCGA_STAR_log2TPMplus1.parquet') if field=='tpm_unstranded' else None
   for index,row,subset in sources:
    columns=['FileID',*subset.gene_id.tolist()]
    offset=0
    for group in range(pf.num_row_groups):
     if offset+pf.metadata.row_group(group).num_rows>index:break
     offset+=pf.metadata.row_group(group).num_rows
    table=pf.read_row_group(group,columns=columns).to_pandas();derived=table.loc[table.FileID==row['file_id'],subset.gene_id.tolist()].iloc[0].to_numpy(dtype=float)
    original=subset[field].to_numpy(dtype=float)
    check(dataset+' original values '+row['file_id'],np.allclose(derived,original,rtol=0,atol=1e-12,equal_nan=True))
    if log_file:
     table=log_file.read_row_group(group,columns=columns).to_pandas()
     vector=table.loc[table.FileID==row['file_id'],subset.gene_id.tolist()].iloc[0].to_numpy(dtype=float)
     check('log2(TPM+1) exact transformation '+row['file_id'],np.allclose(vector,np.log2(original+1),rtol=0,atol=1e-12))
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
 write_json(runtime_path('gdc_validation_report.json'),{'created_at':now(),'release':'DR46','input_manifest_sha256':hashes(manifest_path())[0],'checks':checks,'all_passed':all(c['passed'] for c in checks),'processed':sizes})
 layers=json.loads(runtime_path('current_layer_status.json').read_text(encoding='utf-8'))
 layers.update(current_status='complete',validated_at=now())
 write_json(runtime_path('current_layer_status.json'),layers)
 print('PASS',len(checks),'current-layer checks')

if __name__=='__main__':
 try:main()
 except Exception as error:
  layers=json.loads(runtime_path('current_layer_status.json').read_text(encoding='utf-8'))
  layers.update(current_status='validation_failed',validation_error=str(error))
  write_json(runtime_path('current_layer_status.json'),layers)
  write_json(runtime_path('gdc_validation_report.json'),{'created_at':now(),'all_passed':False,'error':str(error)})
  raise

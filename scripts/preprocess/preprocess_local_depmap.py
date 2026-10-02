"""Stream verified local exports to ZSTD Parquet; keep raw untouched."""
from pathlib import Path
import csv,json,re,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,hashes,now,write_json,setup
import pyarrow as pa
import pyarrow.csv as pcsv
import pyarrow.parquet as pq
MAP={'CopyNumber_WGS_26Q1.csv':'OmicsCNGeneWGS','CRISPR_Chronos_26Q1.csv':'CRISPRGeneEffect',
 'CRISPR_GeneDependency_26Q1.csv':'CRISPRGeneDependency','Expression_26Q1.csv':'ExpressionProteinCoding',
 'Mutation_Damaging_26Q1.csv':'MutationDamaging','Mutation_Hotspot_26Q1.csv':'MutationHotspot',
 'OmicsSignatures_26Q1.csv':'GlobalSignatures','MolecularSubtypes_26Q1.csv':'MolecularSubtypes'}
MATRICES={'OmicsCNGeneWGS','CRISPRGeneEffect','CRISPRGeneDependency','ExpressionProteinCoding','MutationDamaging','MutationHotspot'}

def main():
 setup('preprocess_local_depmap');folder=ROOT/'data/processed/depmap/26Q1';folder.mkdir(parents=True,exist_ok=True)
 with (ROOT/'data/manifests/depmap_26Q1_manifest.csv').open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
 qc=[]
 for row in rows:
  if row['active'].lower()!='true':continue
  source=ROOT/row['path'];sha,md5=hashes(source)
  if sha!=row['sha256']:raise ValueError('Local export SHA256 changed; no overwrite')
  dataset=MAP.get(row['canonical_name'],Path(row['canonical_name']).stem);target=folder/f'{dataset}.parquet';meta=folder/f'{dataset}.parquet.provenance.json'
  with source.open(encoding='utf-8-sig',newline='') as f:original=next(csv.reader(f))
  names=list(original)
  if dataset in MATRICES:names[0]='ModelID'
  if len(set(names))!=len(names):raise ValueError('Duplicate columns in raw export')
  types={n:pa.float64() if dataset in MATRICES and n!='ModelID' else pa.string() for n in names}
  if not (target.exists() and meta.exists() and json.loads(meta.read_text())['raw_sha256']==sha and hashes(target)[0]==json.loads(meta.read_text())['SHA256']):
   reader=pcsv.open_csv(source,read_options=pcsv.ReadOptions(column_names=names,skip_rows=1,block_size=32*1024*1024),convert_options=pcsv.ConvertOptions(column_types=types,strings_can_be_null=True))
   temp=target.with_suffix('.parquet.part')
   with pq.ParquetWriter(temp,reader.schema,compression='zstd',compression_level=6) as writer:
    for batch in reader:writer.write_batch(batch)
   temp.replace(target)
   write_json(meta,{'release':'26Q1','raw_sha256':sha,'SHA256':hashes(target)[0],'raw_source':row['path'],'created_at':now(),'source_type':'user_supplied_portal_export','full_release_completeness':'unverified','orientation':'models x gene symbols' if dataset in MATRICES else 'source metadata table'})
  pf=pq.ParquetFile(target);entry={'dataset':dataset,'rows':pf.metadata.num_rows,'columns':len(pf.schema_arrow),'raw_SHA256':sha}
  if dataset in MATRICES:
   ids=pq.read_table(target,columns=['ModelID'])['ModelID'].to_pylist()
   if len(ids)!=len(set(ids)) or any(not v.startswith('ACH-') for v in ids):raise ValueError('Ambiguous/invalid ModelID in export')
   entry.update(unique_models=len(set(ids)),genes=len(names)-1)
  qc.append(entry);print(json.dumps(entry),flush=True)
 write_json(ROOT/'data/manifests/depmap_26Q1_local_qc.json',{'created_at':now(),'datasets':qc})

if __name__=='__main__':main()

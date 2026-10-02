"""Inspect user-supplied exports; move without changing a byte; keep duplicates."""
from pathlib import Path
import csv, hashlib, json, re, shutil, datetime
ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / 'data/raw/depmap/26Q1'
ARCHIVE = ROOT / 'data/archive/duplicates'
NAMES = {
 'Copy_Number_WGS': 'CopyNumber_WGS_26Q1.csv',
 'CRISPR_(DepMap_Public': 'CRISPR_Chronos_26Q1.csv',
 'CRISPR_Gene_Dependency': 'CRISPR_GeneDependency_26Q1.csv',
 'Expression_(Short-read)': 'Expression_26Q1.csv',
 'Damaging_Mutations': 'Mutation_Damaging_26Q1.csv',
 'Hotspot_Mutations': 'Mutation_Hotspot_26Q1.csv',
 'Omics_Signatures': 'OmicsSignatures_26Q1.csv',
 'Inferred_Molecular_Subtypes': 'MolecularSubtypes_26Q1.csv',
 **{n:n+'.csv' for n in ['Model','ModelCondition','OmicsProfiles','CRISPRScreenMap','ScreenSequenceMap']}}

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()

def inspect(path, canonical):
 ids=set(); rows=0
 with path.open(encoding='utf-8-sig',newline='') as f:
  reader=csv.reader(f); header=next(reader)
  idx=header.index('ModelID') if 'ModelID' in header else 0
  for r in reader:
   if not r:continue
   if len(r)!=len(header):raise ValueError(f'Ragged CSV {path.name}, row {rows+2}')
   rows+=1
   if r[idx].startswith('ACH-'):ids.add(r[idx])
 genes=len(header)-1 if canonical.startswith(('CopyNumber_','CRISPR_Chronos','CRISPR_GeneDependency','Expression_','Mutation_')) else 0
 return dict(canonical_name=canonical,original_name=path.name,path='',size_MB=round(path.stat().st_size/1e6,3),
  rows=rows,columns=len(header),model_count=len(ids),gene_count=genes,sha256=sha(path),active=False,
  modified_at=datetime.datetime.fromtimestamp(path.stat().st_mtime).isoformat(),
  notes='User-supplied 26Q1 portal export. Filename subsetted does not establish whole-genome completeness.' if genes else 'User-supplied metadata/export.',
  _source=path,_ids=ids,_header=header)

def main():
 RAW.mkdir(parents=True,exist_ok=True);ARCHIVE.mkdir(parents=True,exist_ok=True)
 groups={}
 for p in ROOT.glob('*.csv'):
  canonical=next((v for k,v in NAMES.items() if p.stem==k or p.name.startswith(k+'_') or p.name.startswith(k+'(')),None)
  if canonical:
   row=inspect(p,canonical);groups.setdefault(canonical,[]).append(row)
   print(p.name,row['rows'],row['columns'],row['model_count'],row['gene_count'],flush=True)
 if not groups:raise SystemExit('No root exports found; existing manifest retained.')
 out=[];ids={};headers={}
 for name,rows in groups.items():
  rows.sort(key=lambda x:(x['model_count']*max(x['gene_count'],1),x['model_count'],x['gene_count'],x['rows'],x['columns'],x['modified_at']),reverse=True)
  for i,row in enumerate(rows):
   active=i==0;source=row.pop('_source');model_ids=row.pop('_ids');header=row.pop('_header')
   target=RAW/name if active else ARCHIVE/source.name
   if target.exists():raise FileExistsError(f'Will not overwrite {target}')
   if not target.resolve().is_relative_to(ROOT.resolve()):raise ValueError('Target outside project')
   shutil.move(str(source),str(target))
   if sha(target)!=row['sha256']:raise ValueError('SHA changed after move')
   row.update(active=active,path=target.relative_to(ROOT).as_posix())
   if len(rows)>1:row['notes']+=' Active ranked by model/gene coverage; all duplicates preserved.'
   out.append(row)
   if active:ids[name]=model_ids;headers[name]=header
 manifest=ROOT/'data/manifests/depmap_26Q1_manifest.csv'
 with manifest.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
 matrices={k:v for k,v in ids.items() if len(v)>0}
 overlaps=[dict(file_a=a,file_b=b,overlap=len(matrices[a]&matrices[b])) for a in matrices for b in matrices if a<b]
 report={'files':out,'pairwise_model_overlap':overlaps,'headers':headers,
  'duplicate_files':sum(len(v)-1 for v in groups.values()),'note':'Analyze all supplied columns; completeness of portal subsets relative to official full releases is not asserted.'}
 (ROOT/'docs/DEPMAP_INSPECTION.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
 config=ROOT/'config/current_release.json';c=json.loads(config.read_text(encoding='utf-8'))
 c.update(installation_status='local_exports_available',download_date=None,imported_at=datetime.datetime.now().isoformat(),note='User-supplied 26Q1 exports; original download date and export URLs not provided. SHA256 recorded; contents unchanged.')
 config.write_text(json.dumps(c,ensure_ascii=False,indent=2),encoding='utf-8')

if __name__=='__main__':main()

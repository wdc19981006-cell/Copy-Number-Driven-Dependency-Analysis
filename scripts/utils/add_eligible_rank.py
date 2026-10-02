"""Decorate existing screen CSVs without recomputing or changing numeric tokens."""
from pathlib import Path
import argparse,csv,json,sys
from decimal import Decimal
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT

def read(path):
 with path.open(encoding='utf-8',newline='') as source:
  reader=csv.DictReader(source);return list(reader.fieldnames),list(reader)

def write(path,fields,rows):
 temp=path.with_suffix(path.suffix+'.tmp')
 with temp.open('w',encoding='utf-8',newline='') as output:
  writer=csv.DictWriter(output,fieldnames=fields);writer.writeheader();writer.writerows(rows)
 temp.replace(path)

def decorate(case,candidate):
 folder=ROOT/'results'/case/'02_GenomeWide_Dependency'
 if not folder.resolve().is_relative_to((ROOT/'results').resolve()):raise ValueError('Result path outside project')
 path=folder/'GenomeWide_Dependency.csv';fields,rows=read(path)
 original=[{k:v for k,v in r.items() if k not in {'Eligible_Rank','Eligible_N'}} for r in rows]
 if len({r['Gene'] for r in rows})!=len(rows):raise ValueError('Duplicate screen gene')
 ordered=sorted(rows,key=lambda r:int(r['Rank']))
 if [int(r['Rank']) for r in ordered]!=list(range(1,len(rows)+1)):raise ValueError('Original Rank is not complete/unique')
 eligible=[r for r in rows if r['Wilcoxon_FDR'].strip().lower() not in {'','na','nan'}]
 eligible.sort(key=lambda r:(Decimal(r['Wilcoxon_FDR']),Decimal(r['Delta_median']),int(r['Rank'])))
 n=len(eligible);mapping={r['Gene']:str(i+1) for i,r in enumerate(eligible)}
 for r in rows:r['Eligible_Rank']=mapping.get(r['Gene'],'');r['Eligible_N']=str(n)
 if 'Eligible_Rank' not in fields:fields.append('Eligible_Rank')
 if 'Eligible_N' not in fields:fields.append('Eligible_N')
 write(path,fields,rows)
 _,saved=read(path)
 if original!=[{k:v for k,v in r.items() if k not in {'Eligible_Rank','Eligible_N'}} for r in saved]:raise AssertionError('Existing fields changed')
 for name in [f'{candidate}_Candidate_Rank.csv','GenomeWide_Dependency_Top100.csv']:
  target=folder/name
  if not target.exists():continue
  columns,subset=read(target)
  for row in subset:row['Eligible_Rank']=mapping.get(row['Gene'],'');row['Eligible_N']=str(n)
  if 'Eligible_Rank' not in columns:columns.append('Eligible_Rank')
  if 'Eligible_N' not in columns:columns.append('Eligible_N')
  if name==f'{candidate}_Candidate_Rank.csv':
   preferred=['Gene','Eligible_Rank','Eligible_N','Rank','Wilcoxon_FDR','Wilcoxon_P','Delta_median','Pearson_r','Pearson_FDR']
   columns=preferred+[c for c in columns if c not in preferred]
  write(target,columns,subset)
 result=[r for r in rows if r['Gene']==candidate]
 print(json.dumps({'case':case,'eligible_tests':n,'candidate':[{'Gene':r['Gene'],'Rank':r['Rank'],'Eligible_Rank':r['Eligible_Rank']} for r in result]},indent=2))

def main():
 p=argparse.ArgumentParser();p.add_argument('--case',default='VPS4B_VPS4A');p.add_argument('--candidate',default='VPS4A');a=p.parse_args();decorate(a.case,a.candidate)

if __name__=='__main__':main()

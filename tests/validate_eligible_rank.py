"""Compare rank decoration to the committed original CSV fields, token for token."""
from pathlib import Path
import argparse,csv,io,subprocess,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.utils.common import ROOT,write_json,now

def parse(text):return list(csv.DictReader(io.StringIO(text)))

def main():
 p=argparse.ArgumentParser();p.add_argument('--base',default='0e09076c276e9f988966676dfead3d2865f8302f');a=p.parse_args()
 checks=[];folder='results/VPS4B_VPS4A/02_GenomeWide_Dependency/'
 for name in ['GenomeWide_Dependency.csv','VPS4A_Candidate_Rank.csv','GenomeWide_Dependency_Top100.csv']:
  path=folder+name
  before=parse(subprocess.check_output(['git','show',a.base+':'+path],cwd=ROOT).decode('utf-8'))
  after=parse((ROOT/path).read_text(encoding='utf-8'))
  clean=[{k:v for k,v in r.items() if k!='Eligible_Rank'} for r in after]
  assert before==clean,'Existing CSV fields changed: '+name
  checks.append({'file':path,'rows':len(after),'all_original_field_tokens_unchanged':True})
  if name=='GenomeWide_Dependency.csv':
   eligible=[r for r in after if r['Wilcoxon_FDR'].lower() not in {'','na','nan'}]
   assert [int(r['Eligible_Rank']) for r in eligible]==list(range(1,len(eligible)+1))
   assert all(r['Eligible_Rank']=='' for r in after if r['Wilcoxon_FDR'].lower() in {'','na','nan'})
   candidate=next(r for r in after if r['Gene']=='VPS4A')
   assert candidate['Rank']=='276' and candidate['Eligible_Rank']=='1'
 write_json(ROOT/'data/manifests/eligible_rank_validation.json',{'checked_at':now(),'base_commit':a.base,'all_passed':True,'checks':checks,'candidate':{'Gene':'VPS4A','Rank':276,'Eligible_Rank':1}})
 print('PASS: all original CSV field tokens preserved; eligible ordering validated; VPS4A Rank 276 / Eligible_Rank 1.')

if __name__=='__main__':main()

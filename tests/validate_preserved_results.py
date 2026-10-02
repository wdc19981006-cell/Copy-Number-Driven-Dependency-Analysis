"""Preserve old outputs and verify original statistics against the completed commit."""
from pathlib import Path
import argparse,csv,io,json,subprocess,sys,shutil
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.utils.common import ROOT,now,write_json
BASE='d07ad561ed91064182949ca99f185d8e156818b3'
def main():
 p=argparse.ArgumentParser();p.add_argument('--capture',action='store_true');a=p.parse_args()
 folder='results/VPS4B_VPS4A/'
 if a.capture:
  target=ROOT/'.runtime/prior_case'/BASE
  if not target.exists():shutil.copytree(ROOT/folder,target)
  print('Preserved all prior case outputs:',target);return
 names=subprocess.check_output(['git','ls-tree','-r','--name-only',BASE,folder],cwd=ROOT).decode().splitlines()
 checks=[]
 for name in names:
  if not name.endswith('.csv') or '/Summary/' in name:continue
  old=list(csv.DictReader(io.StringIO(subprocess.check_output(['git','show',f'{BASE}:{name}'],cwd=ROOT).decode('utf-8'))))
  with (ROOT/name).open(encoding='utf-8',newline='') as f:new=list(csv.DictReader(f))
  # New rank columns/header order may differ; original field tokens may not.
  assert len(old)==len(new),name
  for before,after in zip(old,new):
   assert all(after.get(k)==v for k,v in before.items()),'Original field changed: '+name
  checks.append({'file':name,'rows':len(old),'all_original_field_tokens_unchanged':True})
 write_json(ROOT/'data/manifests/original_statistics_preservation.json',{'checked_at':now(),'baseline_commit':BASE,'all_passed':True,'checks':checks,'prior_outputs_preserved_locally':str((ROOT/'.runtime/prior_case'/BASE).relative_to(ROOT))})
 print(f'PASS: {len(checks)} prior CSVs retain all original row/field tokens.')
if __name__=='__main__':main()

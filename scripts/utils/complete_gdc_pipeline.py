"""Continue authorized download -> ETL -> validation as a local, resumable job.

Live status stays under ignored raw/manifests/live. Repository manifests are
explicitly timestamped snapshots; publishing does not silently claim readiness.
"""
from pathlib import Path
import sys,subprocess,time,json,traceback
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,now,write_json
from scripts.utils.gdc import runtime_path,report_path

def run(script,args=()):
 log=ROOT/'logs/gdc_pipeline.log'
 with log.open('a',encoding='utf-8') as output:
  p=subprocess.run([sys.executable,'-u',str(ROOT/script),*args],cwd=ROOT,stdout=output,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 if p.returncode:raise RuntimeError(f'{script}: exit {p.returncode}; see logs/gdc_pipeline.log')

def main():
 status=runtime_path('pipeline_status.json')
 try:
  for attempt in range(1,4):
   write_json(status,{'updated_at':now(),'phase':'downloading','attempt':attempt,'note':'Current TCGA unavailable until ETL and validation complete.'})
   try:run('scripts/download/download_gdc_bundles.py',['--workers','8','--batch-size','32']);break
   except RuntimeError:
    if attempt==3:raise
  report=json.loads(report_path('gdc_download_summary.json').read_text())
  if any(v['verified']!=v['files'] for v in report['counts'].values()):raise ValueError('Incomplete download; no ETL readiness claim')
  write_json(status,{'updated_at':now(),'phase':'preprocessing'})
  run('scripts/preprocess/preprocess_gdc.py')
  write_json(status,{'updated_at':now(),'phase':'validating'})
  run('scripts/utils/validate_gdc_current.py')
  write_json(status,{'updated_at':now(),'phase':'complete','note':'Raw downloaded, publisher MD5/SHA256 verified, current ETL and validation complete. Repository/case snapshots remain timestamped; rerun desired TCGA modes to update example results.'})
 except Exception as e:
  write_json(status,{'updated_at':now(),'phase':'failed','error':str(e),'traceback':traceback.format_exc(),'note':'Success is not claimed. Verified raw and partial downloads retained. Rerun this script to resume.'})
  raise

if __name__=='__main__':main()

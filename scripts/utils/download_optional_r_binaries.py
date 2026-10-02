from pathlib import Path
import sys,zipfile,json
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,curl,hashes,remote_info,write_json
folder=ROOT/'.runtime/R-optional';folder.mkdir(parents=True,exist_ok=True)
report=[]
for p,v in [('arrow','25.0.1'),('duckdb','1.5.5')]:
 url=f'https://cloud.r-project.org/bin/windows/contrib/4.5/{p}_{v}.zip';target=folder/f'{p}_{v}.zip';part=target.with_suffix('.zip.part')
 headers=remote_info(url);size=int(headers['content-length']) if headers.get('content-length') else None
 if not target.exists():
  curl(['--max-time','600','--retry','5','--continue-at','-','--output',part,url])
  if size and part.stat().st_size!=size:raise ValueError(f'{p}: incomplete binary ZIP')
  with zipfile.ZipFile(part) as z:
   if z.testzip():raise ValueError('ZIP CRC failed')
  part.rename(target)
 report.append({'package':p,'version':v,'url':url,'bytes':target.stat().st_size,'SHA256':hashes(target)[0]})
 print(p,target.stat().st_size,'ZIP CRC passed',flush=True)
write_json(ROOT/'data/manifests/R_optional_downloads.json',report)

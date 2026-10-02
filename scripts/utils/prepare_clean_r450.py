"""Extract an official R 4.5.0 runtime locally, without running its installer."""
from pathlib import Path
import sys,zipfile,subprocess,json
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,curl,fetch_text,hashes,write_json
folder=ROOT/'tools/R450';folder.mkdir(parents=True,exist_ok=True)
urls={'R-4.5.0-win.exe':'https://cran.r-project.org/bin/windows/base/old/4.5.0/R-4.5.0-win.exe',
 'innoextract.zip':'https://constexpr.org/innoextract/files/innoextract-1.9-windows.zip'}
for name,url in urls.items():
 target=folder/name
 if not target.exists():curl(['--retry','3','--output',target,url])
md5=fetch_text('https://cran.r-project.org/bin/windows/base/old/4.5.0/md5sum.R-4.5.0.txt').split()[0]
if hashes(folder/'R-4.5.0-win.exe')[1]!=md5:raise ValueError('Official R installer MD5 mismatch')
with zipfile.ZipFile(folder/'innoextract.zip') as z:
 for member in z.infolist():
  target=(folder/member.filename).resolve()
  if not target.is_relative_to(folder.resolve()):raise ValueError('Unsafe extractor ZIP')
  z.extract(member,folder)
exe=next(folder.rglob('innoextract.exe'))
result=subprocess.run([str(exe),'--extract','--output-dir',str(folder/'runtime'),str(folder/'R-4.5.0-win.exe')],capture_output=True)
(ROOT/'logs/R450_extract.log').write_bytes(result.stdout+result.stderr)
print('Extract exit',result.returncode)
if result.returncode:raise SystemExit('See logs/R450_extract.log; existing system R remains untouched')
candidate=next((folder/'runtime').rglob('Rscript.exe'))
write_json(ROOT/'config/R_runtime.json',{'version':'4.5.0','rscript':candidate.relative_to(ROOT).as_posix(),'source':urls['R-4.5.0-win.exe'],'official_md5':md5,'installer_sha256':hashes(folder/'R-4.5.0-win.exe')[0]})
print(candidate)

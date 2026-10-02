from pathlib import Path
import psutil,json,zipfile
ROOT=Path(__file__).resolve().parents[2]
for p in psutil.process_iter(['pid','name','cmdline','memory_info','cpu_times']):
 try:
  if not any(str(ROOT).lower() in str(v).lower() or 'scripts/' in str(v) for v in p.info['cmdline'] or []):continue
  if p.info['name'].lower() not in ('python.exe','rterm.exe','rscript.exe'):continue
  print(json.dumps({'pid':p.pid,'name':p.info['name'],'args':p.info['cmdline'][-6:],'RSS_GiB':round(p.info['memory_info'].rss/2**30,3),'cpu_seconds':round(sum(p.info['cpu_times'][:2]),1)}))
 except psutil.Error:pass
print('Available memory GiB',round(psutil.virtual_memory().available/2**30,3))
for p in (ROOT/'.runtime/renv-root/binary').rglob('*.zip'):
 if p.name.startswith(('arrow_','duckdb_')):print(p.name,p.stat().st_size,'ZIP',zipfile.is_zipfile(p))

"""Inspect exactly staged files before the authorized public initial commit."""
from pathlib import Path
import subprocess,sys,json,re
ROOT=Path(__file__).resolve().parents[2]
names=subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode('utf-8').split('\0')
names=[name for name in names if name]
bad=[];largest=[]
for name in names:
 p=ROOT/name;size=p.stat().st_size;largest.append((size,name))
 if size>50_000_000:bad.append((name,'larger than 50 MB'))
 if name.startswith(('data/raw/','data/processed/','data/archive/','logs/','.runtime/','tools/','renv/library/','renv/staging/')):bad.append((name,'excluded data/runtime directory'))
 if p.suffix.lower() in {'.parquet','.duckdb','.rds','.rdata','.exe','.zip'}:bad.append((name,'database/binary extension'))
 if name.startswith('results/') and not name.startswith('results/VPS4B_VPS4A/'):bad.append((name,'non-case result'))
 if p.suffix.lower() in {'.py','.r','.json','.yaml','.md','.txt'}:
  value=p.read_text(encoding='utf-8',errors='replace')
  if re.search(r'gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|sk-[A-Za-z0-9]{30,}|-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----',value):bad.append((name,'credential pattern; inspect locally'))
print(json.dumps({'staged_files':len(names),'total_bytes':sum(v[0] for v in largest),'largest_files':sorted(largest,reverse=True)[:8],'issues':bad},indent=2))
if bad:raise SystemExit('Staged audit failed; do not commit/push')

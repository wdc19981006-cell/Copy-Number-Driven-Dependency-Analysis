"""Compact progress snapshot; presence/size counts are not publisher MD5 QC."""
from pathlib import Path
import json
import shutil
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT
from scripts.utils.gdc import RAW,local_path

out={}
for kind in ['rna','cn','segments','mutation']:
    rows=list({r['file_id']:r for r in json.loads((RAW/f'metadata/{kind}_selected.json').read_text(encoding='utf-8'))}.values())
    done=[r for r in rows if local_path(kind,r).exists() and local_path(kind,r).stat().st_size==r['file_size']]
    parts=[local_path(kind,r).with_name(r['file_name']+'.curl.part') for r in rows]
    partial_complete=sum(p.exists() and p.stat().st_size==r['file_size'] for p,r in zip(parts,rows))
    out[kind]={'expected':len(rows),'final_files_present':len(done),'batch_files_waiting_hash':partial_complete,
               'final_bytes':sum(r['file_size'] for r in done)}
from scripts.utils.gdc import report_path
summary=report_path('gdc_download_summary.json')
if summary.exists(): out['verified']=json.loads(summary.read_text(encoding='utf-8'))['counts']
out['disk_free_GiB']=round(shutil.disk_usage(ROOT).free/1024**3,3)
print(json.dumps(out,indent=2))

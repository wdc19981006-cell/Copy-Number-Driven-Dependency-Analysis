"""Recheck pinned raw size/MD5/SHA256; never download, delete or modify raw."""
from pathlib import Path
import argparse,csv,json,sys,hashlib,logging
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,now,write_json,hashes
from scripts.utils.gdc import RAW,manifest_path,runtime_path,local_path

def verify_raw(kinds=('rna','cn','segments','mutation'),workers=4):
 with manifest_path().open(encoding='utf-8',newline='') as source:
  records={r['file_id']:r for r in csv.DictReader(source)}
 report_path=runtime_path('gdc_raw_integrity.json')
 report=json.loads(report_path.read_text(encoding='utf-8')) if report_path.exists() else {'kinds':{}}
 for kind in kinds:
  selected=json.loads((RAW/f'metadata/{kind}_selected.json').read_text(encoding='utf-8'))
  rows=list({r['file_id']:r for r in selected}.values())
  fingerprint=hashlib.sha256();inputs=[]
  for row in rows:
   record=records[row['file_id']];path=local_path(kind,row);stat=path.stat()
   if record['status']!='verified' or stat.st_size!=int(row['file_size']):raise ValueError(f'{kind} raw incomplete: {row["file_id"]}')
   if record['MD5']!=row['md5sum'] or len(record['SHA256'])!=64:raise ValueError('Invalid pinned checksum')
   fingerprint.update(f"{row['file_id']}:{stat.st_size}:{stat.st_mtime_ns}:{record['SHA256']}:{row['md5sum']}\n".encode())
   inputs.append((path,row,record,stat.st_mtime_ns))
  digest=fingerprint.hexdigest();previous=report['kinds'].get(kind,{})
  if previous.get('all_passed') and previous.get('raw_fingerprint')==digest:
   logging.info('Raw integrity audit unchanged; reuse %s (%s files)',kind,len(inputs));continue
  def validate(item):
   path,row,record,mtime=item;sha,md5=hashes(path)
   if sha!=record['SHA256'] or md5!=row['md5sum'] or path.stat().st_mtime_ns!=mtime:
    raise ValueError(f'Raw checksum/identity changed: {kind} {row["file_id"]}; preserved unchanged')
   return int(row['file_size'])
  total=0
  with ThreadPoolExecutor(max_workers=workers) as pool:
   for n,size in enumerate(pool.map(validate,inputs),1):
    total+=size
    if n%1000==0:logging.info('Raw MD5/SHA256 rechecked %s %s/%s',kind,n,len(inputs))
  report['kinds'][kind]={'checked_at':now(),'files':len(inputs),'bytes':total,'all_passed':True,'raw_fingerprint':digest,'method':'Actual raw MD5 and SHA256; reuse only when all pinned checksum/size/mtime fingerprints remain unchanged'}
  write_json(report_path,report)
 report['updated_at']=now();report['all_passed']=all(report['kinds'].get(k,{}).get('all_passed',False) for k in ['rna','cn','segments','mutation'])
 write_json(report_path,report)
 return report

def main():
 p=argparse.ArgumentParser();p.add_argument('--kinds',nargs='+',choices=['rna','cn','segments','mutation'],default=['rna','cn','segments','mutation']);a=p.parse_args()
 logging.basicConfig(level=logging.INFO,format='%(asctime)s %(message)s')
 print(json.dumps(verify_raw(a.kinds),indent=2))

if __name__=='__main__':main()

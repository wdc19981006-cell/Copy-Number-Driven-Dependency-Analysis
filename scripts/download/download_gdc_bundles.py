"""Official GDC multi-UUID POST archives; selected UUIDs, per-file MD5 + SHA256.

Uses official documented bundling to reduce per-request latency/compress text.
Only explicitly selected regular archive members are copied; no extractall.
Existing complete raw remains immutable. Interrupted per-file downloads remain.
"""
from pathlib import Path
import argparse, csv, json, subprocess, sys, tarfile, logging, time, hashlib, shutil
from concurrent.futures import ThreadPoolExecutor,as_completed
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,setup,now,curl_binary,fetch_text
from scripts.utils.gdc import RAW,API,local_path,manifest_path
from scripts.download.download_gdc_api import base,save_manifest,summary,finalize

def transfer(kind,batch,number):
 directory=RAW/'_transfer';directory.mkdir(parents=True,exist_ok=True)
 payload=directory/f'{kind}_{number:06d}.json';archive=directory/f'{kind}_{number:06d}.tar.gz.part'
 payload.write_text(json.dumps({'ids':[r['file_id'] for r in batch]}),encoding='ascii')
 began=time.monotonic()
 log=ROOT/f'logs/gdc_bundle_{kind}_{number:06d}.log'
 with log.open('a',encoding='utf-8') as output:
  result=subprocess.run([curl_binary(),'--fail','--location','--silent','--show-error','--connect-timeout','30',
   '--max-time','1800','--retry','4','--retry-delay','3','--retry-all-errors','--request','POST','--header','Content-Type: application/json',
   '--data-binary','@'+payload.relative_to(RAW).as_posix(),'--output',archive.relative_to(RAW).as_posix(),API+'/data'],
   cwd=str(RAW),stdout=output,stderr=subprocess.STDOUT)
 if result.returncode:raise RuntimeError(f'{kind} archive {number}: curl exit {result.returncode}; retained {archive.name}')
 expected={r['file_id']:r for r in batch};found={};results=[]
 with tarfile.open(archive,'r|gz') as tar:
  for member in tar:
   if not member.isfile():continue
   parts=Path(member.name).parts
   ids=[v for v in parts if v in expected]
   if not ids:continue # archive metadata / manifest
   if len(ids)!=1:raise ValueError('Ambiguous archive UUID')
   row=expected[ids[0]]
   if Path(member.name).name!=row['file_name']:continue
   if row['file_id'] in found:raise ValueError('Duplicate archive data member')
   if member.size!=row['file_size']:raise ValueError('Archive source size mismatch')
   target=local_path(kind,row);target.parent.mkdir(parents=True,exist_ok=True)
   if not target.resolve().is_relative_to(RAW.resolve()):raise ValueError('Target outside current layer')
   if target.exists():raise FileExistsError('Concurrent raw write refused')
   partial=target.with_name(target.name+'.bundle.part');h=hashlib.sha256();md5=hashlib.md5()
   stream=tar.extractfile(member)
   with partial.open('wb') as out:
    for block in iter(lambda:stream.read(4*1024*1024),b''):out.write(block);h.update(block);md5.update(block)
   if md5.hexdigest()!=row['md5sum']:raise ValueError(f'Publisher checksum failed {row["file_id"]}; partial retained')
   partial.rename(target);found[row['file_id']]=True
   results.append({**base(kind,row),'status':'verified','SHA256':h.hexdigest(),'MD5':md5.hexdigest(),'verified_at':now(),'error':''})
 if set(found)!=set(expected):raise ValueError(f'Archive missing {len(set(expected)-set(found))} files')
 archive_bytes=archive.stat().st_size
 # Remove only generated, successfully consumed transfer artifacts within RAW.
 for generated in (archive,payload):
  if not generated.resolve().is_relative_to((RAW/'_transfer').resolve()):raise ValueError('Unsafe transfer cleanup')
  generated.unlink()
 return results,time.monotonic()-began,archive_bytes

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--kinds',nargs='+',default=['segments','mutation','rna','cn']);parser.add_argument('--workers',type=int,default=8);parser.add_argument('--batch-size',type=int,default=256);parser.add_argument('--limit',type=int)
 args=parser.parse_args();setup('download_gdc_bundles')
 current=json.loads(fetch_text(API+'/status'))
 if current.get('data_release_version',{}).get('major')!=46:raise ValueError('Release changed')
 records={}
 manifest=manifest_path()
 if manifest.exists():
  with manifest.open(encoding='utf-8',newline='') as f:records={r['file_id']:r for r in csv.DictReader(f)}
 tasks=[]
 for kind in args.kinds:
  rows=json.loads((RAW/f'metadata/{kind}_selected.json').read_text(encoding='utf-8'))
  rows=list({r['file_id']:r for r in rows}.values());pending=[]
  for row in rows:
   path=local_path(kind,row)
   previous=records.get(row['file_id'],{})
   curl_part=path.with_name(path.name+'.curl.part')
   if path.exists() or (curl_part.exists() and curl_part.stat().st_size==int(row['file_size'])):
    if path.exists() and previous.get('status')=='verified' and path.stat().st_size==int(row['file_size']):continue
    record=finalize(kind,row,records)
    if record['status']!='verified':raise ValueError('Existing raw integrity failure; will not overwrite')
   else:pending.append(row)
  if args.limit:pending=pending[:args.limit]
  for number,start in enumerate(range(0,len(pending),args.batch_size)):tasks.append((kind,pending[start:start+args.batch_size],number))
 if sum(sum(r['file_size'] for r in batch) for _,batch,_ in tasks)>shutil.disk_usage(ROOT).free*.8:raise ValueError('Insufficient disk reserve')
 save_manifest(records);summary(records)
 failed=[]
 with ThreadPoolExecutor(max_workers=args.workers) as pool:
  jobs={pool.submit(transfer,*task):task for task in tasks}
  for future in as_completed(jobs):
   kind,batch,number=jobs[future]
   try:
    result,seconds,compressed=future.result()
    for row in result:records[row['file_id']]=row
    logging.info('%s bundle %s: %s files MD5 verified in %.1fs; compressed %.1f MB',kind,number,len(result),seconds,compressed/1e6)
   except Exception as e:
    failed.append(str(e));logging.error('%s',e)
    # Successful members before an archive failure are rechecked and checkpointed.
    for row in batch:
     if local_path(kind,row).exists():finalize(kind,row,records)
   save_manifest(records);summary(records)
 logging.info('Complete; bundle failures=%s',len(failed));print(json.dumps(summary(records)['counts'],indent=2))
 return bool(failed)

if __name__=='__main__':raise SystemExit(main())

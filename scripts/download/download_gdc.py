"""Download selected official manifests with GDC DTT; all UUIDs pinned to DR46."""
from pathlib import Path
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import csv
import json
import logging
import subprocess
import shutil
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,setup,hashes,write_json,now,fetch_text
from scripts.utils.gdc import RAW,API,manifest,local_path,runtime_path


def selected(kind):
    return json.loads((RAW/f'metadata/{kind}_selected.json').read_text(encoding='utf-8'))


def scan(kind,rows):
    out=[]
    for row in {r['file_id']:r for r in rows}.values():
        path=local_path(kind,row)
        state='not_downloaded'
        sha,md5='',''
        if path.exists():
            sha,md5=hashes(path)
            state='verified' if path.stat().st_size==row['file_size'] and md5==row['md5sum'] else 'checksum_failed_preserved'
        out.append({**row,'kind':kind,'original_url':f'{API}/data/{row["file_id"]}',
                    'relative_path':path.relative_to(ROOT).as_posix(), 'SHA256':sha,'MD5':md5,
                    'status':state,'verified_at':now() if state=='verified' else ''})
    return out


def run_shard(executable,kind,number,rows,connections):
    pending=ROOT/f'logs/gdc_{kind}_shard_{number:02d}.manifest.tsv'
    manifest(pending,rows)
    directory=local_path(kind,rows[0]).parents[1]
    directory.mkdir(parents=True,exist_ok=True)
    log=ROOT/f'logs/gdc_{kind}_shard_{number:02d}.log'
    args=[str(executable),'download','-m',str(pending),'-d',str(directory),'-n',str(connections),
          '--retry-amount','5','--wait-time','2','--no-related-files','--no-annotations','--color_off']
    logging.info('DTT %s shard %s started: %s files',kind,number,len(rows))
    with log.open('a',encoding='utf-8') as output:
        result=subprocess.run(args,stdout=output,stderr=subprocess.STDOUT,check=False)
    logging.info('DTT %s shard %s exit=%s',kind,number,result.returncode)
    return kind,number,result.returncode


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--kinds',nargs='+',default=['rna','cn','segments','mutation'],choices=['rna','cn','segments','mutation'])
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--limit',type=int,help='Download a small schema-validation sample first; default downloads complete manifests.')
    args=parser.parse_args()
    setup('download_gdc')
    summary=json.loads((RAW/'metadata/discovery_summary.json').read_text(encoding='utf-8'))
    if not summary['safe_to_download']:
        raise SystemExit('STOP: discovery disk threshold failed; no files downloaded')
    status=json.loads(fetch_text(API+'/status'))
    if status.get('data_release_version',{}).get('major')!=46:
        raise ValueError('Current release changed; refusing to label another release DR46')
    installed=json.loads((ROOT/'config/gdc_client.json').read_text(encoding='utf-8'))
    executable=ROOT/installed['executable']
    jobs=[]
    all_rows=[]
    for kind in args.kinds:
        rows=selected(kind)
        unique=list({r['file_id']:r for r in rows}.values())
        if args.limit:
            unique=unique[:args.limit]
        verified=scan(kind,unique)
        all_rows.extend(verified)
        bad=[r for r in verified if r['status']=='checksum_failed_preserved']
        if bad:
            raise ValueError(f'{len(bad)} existing raw files fail MD5; preserved, no overwrite')
        pending=[r for r in verified if r['status']=='not_downloaded']
        slots=min(args.workers,len(pending))
        for i in range(slots):
            jobs.append((kind,i,pending[i::slots]))
    pending_bytes=sum(r['file_size'] for _,_,rows in jobs for r in rows)
    free=shutil.disk_usage(ROOT).free
    if pending_bytes>free*0.8:
        raise SystemExit('STOP: remaining raw files exceed current free-space threshold')
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_shard,executable,kind,i,rows,1) for kind,i,rows in jobs]
        codes=[f.result() for f in as_completed(futures)]
    result=[]
    for kind in args.kinds:
        result.extend(scan(kind,selected(kind)))
    columns=list(dict.fromkeys(k for r in result for k in r))
    target=runtime_path('gdc_data_manifest.csv')
    with target.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(result)
    counts={kind:{'files':sum(r['kind']==kind for r in result),
                  'verified':sum(r['kind']==kind and r['status']=='verified' for r in result),
                  'verified_bytes':sum(r['file_size'] for r in result if r['kind']==kind and r['status']=='verified')}
            for kind in args.kinds}
    write_json(runtime_path('gdc_download_summary.json'),{'created_at':now(),'counts':counts,
               'probe_only':bool(args.limit),'DTT_exit_codes':codes,
               'failures':[{'kind':r['kind'],'file_id':r['file_id'],'status':r['status']} for r in result if r['status']!='verified']})
    print(json.dumps(counts,indent=2))
    return not args.limit and any(r['status']!='verified' for r in result)


if __name__=='__main__':
    raise SystemExit(main())

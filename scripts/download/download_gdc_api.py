"""Manifest-driven fallback using the official GDC data API and curl parallel.

The official DTT was attempted first; this avoids its legacy endpoint TLS failure.
UUID + publisher MD5 pin each resource. Only selected open/active TCGA files are
requested. Complete raw files are immutable; partial files support safe resumption.
"""
from pathlib import Path
import argparse
import csv
import json
import logging
import shutil
import subprocess
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,setup,hashes,now,write_json,curl_binary,fetch_text
from scripts.utils.gdc import RAW,API,local_path,runtime_path,manifest_path

FIELDS=['kind','release','file_id','file_name','file_size','publisher_md5','MD5','SHA256',
        'project_id','case_id','case','sample_id','sample','aliquot_id','aliquot','workflow',
        'original_url','relative_path','status','verified_at','error']


def quote(value):
    return '"'+str(value).replace('\\','\\\\').replace('"','\\"')+'"'


def save_manifest(records):
    path=runtime_path('gdc_data_manifest.csv')
    temp=path.with_suffix('.csv.tmp')
    with temp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,extrasaction='ignore');w.writeheader();w.writerows(records.values())
    temp.replace(path)


def base(kind,row):
    path=local_path(kind,row).resolve()
    if not path.is_relative_to(RAW.resolve()):
        raise ValueError('File target leaves the authorized current-layer directory')
    return {**row,'kind':kind,'release':'46.0','publisher_md5':row['md5sum'],
            'original_url':API+'/data/'+row['file_id'],
            'relative_path':path.relative_to(ROOT).as_posix()}


def finalize(kind,row,records):
    path=local_path(kind,row)
    partial=path.with_name(path.name+'.curl.part')
    candidate=path if path.exists() else partial
    result=base(kind,row)
    if not candidate.exists() or candidate.stat().st_size!=row['file_size']:
        result.update(status='failed' if partial.exists() else 'not_downloaded',error='Absent or incomplete download',SHA256='',MD5='')
    else:
        sha,md5=hashes(candidate)
        if md5.lower()!=row['md5sum'].lower():
            result.update(status='checksum_failed_preserved',error='Publisher MD5 mismatch; file preserved',SHA256=sha,MD5=md5)
        else:
            if candidate!=path:
                candidate.rename(path)
            result.update(status='verified',SHA256=sha,MD5=md5,verified_at=now(),error='')
    records[row['file_id']]=result
    return result


def summary(records):
    counts={}
    for kind in ['rna','cn','segments','mutation']:
        rows=[r for r in records.values() if r['kind']==kind]
        good=[r for r in rows if r['status']=='verified']
        counts[kind]={'files':len(rows),'verified':len(good),
                      'verified_bytes':sum(int(r['file_size']) for r in good)}
    report={'created_at':now(),'method':'Official GDC data API: multi-UUID POST bundles / per-file curl fallback; DTT attempted first',
            'counts':counts,'failures':[{'kind':r['kind'],'file_id':r['file_id'],'status':r['status'],'error':r.get('error','')}
                                      for r in records.values() if r['status']!='verified']}
    write_json(runtime_path('gdc_download_summary.json'),report)
    return report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--kinds',nargs='+',default=['rna','cn','segments','mutation'],choices=['rna','cn','segments','mutation'])
    parser.add_argument('--parallel',type=int,default=16)
    parser.add_argument('--batch-size',type=int,default=512)
    args=parser.parse_args()
    setup('download_gdc_api')
    discovery=json.loads((RAW/'metadata/discovery_summary.json').read_text(encoding='utf-8'))
    if not discovery['safe_to_download']: raise SystemExit('STOP: space threshold failed')
    current=json.loads(fetch_text(API+'/status'))
    if current.get('data_release_version',{}).get('major')!=46: raise ValueError('Release changed')
    path=manifest_path()
    records={}
    if path.exists():
        with path.open(encoding='utf-8',newline='') as f: records={r['file_id']:r for r in csv.DictReader(f)}
    selections={}
    for kind in ['rna','cn','segments','mutation']:
        rows=json.loads((RAW/f'metadata/{kind}_selected.json').read_text(encoding='utf-8'))
        selections[kind]=list({r['file_id']:r for r in rows}.values())
        for row in selections[kind]:
            if row['file_id'] not in records:
                records[row['file_id']]={**base(kind,row),'status':'not_downloaded'}
    for kind in args.kinds:
        pending=[]
        for row in selections[kind]:
            file=local_path(kind,row)
            if file.exists():
                result=finalize(kind,row,records)
                if result['status']!='verified':
                    save_manifest(records);summary(records)
                    raise ValueError('Existing complete raw failed integrity; it will not be overwritten')
            else: pending.append(row)
        remaining=sum(r['file_size'] for r in pending)
        if remaining>shutil.disk_usage(ROOT).free*0.8:
            raise SystemExit('STOP: remaining dataset exceeds available space threshold')
        save_manifest(records);summary(records)
        for start in range(0,len(pending),args.batch_size):
            batch=pending[start:start+args.batch_size]
            config=ROOT/f'logs/gdc_api_{kind}_{start:06d}.curl.conf'
            lines=[]
            for i,row in enumerate(batch):
                file=local_path(kind,row)
                file.parent.mkdir(parents=True,exist_ok=True)
                partial=file.with_name(file.name+'.curl.part')
                if i: lines.append('next')
                lines+=['fail','location','silent','show-error','connect-timeout = 30',
                        'max-time = 900','retry = 5','retry-delay = 2','retry-all-errors',
                        'continue-at = "-"','url = '+quote(API+'/data/'+row['file_id']),
                        'output = '+quote(partial.relative_to(RAW).as_posix())]
            config.write_text('\n'.join(lines)+'\n',encoding='utf-8')
            log=ROOT/f'logs/gdc_api_{kind}_{start:06d}.log'
            began=time.monotonic()
            with log.open('a',encoding='utf-8') as output:
                result=subprocess.run([curl_binary(),'--parallel','--parallel-immediate','--parallel-max',str(args.parallel),
                    '--config',str(config)],cwd=str(RAW),stdout=output,stderr=subprocess.STDOUT,check=False)
            failed=0
            for row in batch:
                result=finalize(kind,row,records)
                failed+=result['status']!='verified'
            save_manifest(records)
            report=summary(records)
            logging.info('%s batch %s/%s complete in %.1fs; MD5 verified=%s, failed=%s',kind,
                         min(start+len(batch),len(pending)),len(pending),time.monotonic()-began,
                         report['counts'][kind]['verified'],failed)
    report=summary(records)
    print(json.dumps(report['counts'],indent=2))
    return any(r['status']!='verified' for r in records.values())


if __name__=='__main__':
    raise SystemExit(main())

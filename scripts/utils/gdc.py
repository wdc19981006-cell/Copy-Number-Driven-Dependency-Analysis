"""Official GDC queries and pinned layer paths. No protected/raw sequencing files."""
from __future__ import annotations
import csv
import json
import logging
from pathlib import Path
import sys
import urllib.parse
import hashlib
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT, curl, fetch_text, write_json, now

API = 'https://api.gdc.cancer.gov'
RAW = ROOT / 'data/raw/tcga/gdc_current_DR46'
PROCESSED = ROOT / 'data/processed/tcga/gdc_DR46'
PRIORITY = {'ABSOLUTE LiftOver':0, 'ASCAT3':1, 'AscatNGS':2, 'ASCAT2':3}
LIVE = RAW / 'manifests/live'

def runtime_path(name):
    LIVE.mkdir(parents=True,exist_ok=True)
    return LIVE / name

def manifest_path():
    live=runtime_path('gdc_data_manifest.csv')
    return live if live.exists() else ROOT/'data/manifests/gdc_data_manifest.csv'

def report_path(name):
    live=runtime_path(name)
    return live if live.exists() else ROOT/'data/manifests'/name
FILE_FIELDS = ','.join([
    'file_id','file_name','file_size','md5sum','state','access','data_category','data_type','data_format',
    'analysis.workflow_type','associated_entities',
    'cases.case_id','cases.submitter_id','cases.project.project_id',
    'cases.samples.sample_id','cases.samples.submitter_id','cases.samples.sample_type',
    'cases.samples.portions.portion_id','cases.samples.portions.analytes.analyte_id',
    'cases.samples.portions.analytes.aliquots.aliquot_id',
    'cases.samples.portions.analytes.aliquots.submitter_id',
])


def init_dirs():
    for d in ['manifests','metadata','expression/raw_files','copy_number/raw_files',
              'copy_number_segments/raw_files','mutation/raw_files','clinical']:
        (RAW / d).mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)


def eq(field, value):
    return {'op':'in','content':{'field':field, 'value':value if isinstance(value,list) else [value]}}


def filters(category=None, datatype=None, workflow=None):
    fs = [eq('cases.project.program.name','TCGA'),eq('access','open'),eq('state','released')]
    for field, value in [('data_category',category),('data_type',datatype),('analysis.workflow_type',workflow)]:
        if value:
            fs.append(eq(field,value))
    return {'op':'and','content':fs}


def request(endpoint, params):
    # GET is the official search API; curl avoids local requests/proxy incompatibility.
    url = API + '/' + endpoint + '?' + urllib.parse.urlencode({
        k:json.dumps(v,separators=(',',':')) if isinstance(v,(dict,list)) else str(v)
        for k,v in params.items()})
    return json.loads(fetch_text(url))


def paginated(endpoint, query, name, page_size=500):
    folder = RAW / ('clinical' if endpoint=='cases' else 'metadata') / f'{name}_pages'
    folder.mkdir(parents=True, exist_ok=True)
    rows, offset, total = [], 0, None
    while total is None or offset < total:
        response = request(endpoint,{**query,'size':page_size,'from':offset,'sort':'case_id:asc' if endpoint=='cases' else 'file_id:asc'})
        data = response.get('data',{})
        hits = data.get('hits')
        if hits is None:
            raise ValueError(f'Unexpected GDC response: {response}')
        total_now = data['pagination']['total']
        if total is not None and total_now != total:
            raise ValueError('GDC total changed during pagination; query snapshot not consistent')
        total = total_now
        write_json(folder / f'{offset:06d}.json', response)
        if not hits and offset < total:
            raise ValueError('GDC pagination ended early')
        rows.extend(hits)
        offset += len(hits)
        logging.info('%s %s/%s',name,offset,total)
    idfield = 'case_id' if endpoint=='cases' else 'file_id'
    if len({r[idfield] for r in rows}) != len(rows):
        raise ValueError('Duplicate entity IDs across API pages')
    return rows


def entities(file):
    """Only file-associated samples; never attach every sample from a case."""
    matches = []
    associated = file.get('associated_entities',[])
    associated_ids = {a.get('entity_id') for a in associated}
    associated_barcodes = {a.get('entity_submitter_id') for a in associated}
    for case in file.get('cases',[]):
        project = case.get('project',{}).get('project_id','')
        if not project.startswith('TCGA-'):
            raise ValueError('Non-TCGA case in file query')
        for sample in case.get('samples',[]):
            for portion in sample.get('portions',[]):
                for analyte in portion.get('analytes',[]):
                    for aliquot in analyte.get('aliquots',[]):
                        if (not associated or aliquot.get('aliquot_id') in associated_ids
                                or aliquot.get('submitter_id') in associated_barcodes):
                            matches.append({'project_id':project,'case_id':case['case_id'],
                                'case':case.get('submitter_id',''),'sample_id':sample['sample_id'],
                                'sample':sample.get('submitter_id',''),'sample_type':sample.get('sample_type',''),
                                'aliquot_id':aliquot.get('aliquot_id',''),'aliquot':aliquot.get('submitter_id','')})
    # Some data types associate directly to a sample instead of an aliquot.
    if not matches:
        for case in file.get('cases',[]):
            for sample in case.get('samples',[]):
                if (not associated or sample.get('sample_id') in associated_ids
                        or sample.get('submitter_id') in associated_barcodes):
                    matches.append({'project_id':case.get('project',{}).get('project_id',''),
                        'case_id':case['case_id'],'case':case.get('submitter_id',''),
                        'sample_id':sample['sample_id'],'sample':sample.get('submitter_id',''),
                        'sample_type':sample.get('sample_type',''),'aliquot_id':'','aliquot':''})
    return matches


def file_row(file, entity=None):
    return {**(entity or {}), 'file_id':file['file_id'],'file_name':file['file_name'],
            'file_size':int(file['file_size']), 'md5sum':file['md5sum'],
            'workflow':file.get('analysis',{}).get('workflow_type',''),
            'access':file.get('access',''), 'state':file.get('state',''),
            'data_type':file.get('data_type','')}


def target_entities(file):
    candidates=entities(file)
    if len(candidates)<=1:
        return [{**e,'mapping_method':'unique API-linked biospecimen'} for e in candidates]
    name=file['file_name']
    explicit=[e for e in candidates if (e.get('aliquot_id') and e['aliquot_id'] in name)
              or (e.get('aliquot') and e['aliquot'] in name)]
    if len(explicit)==1:
        return [{**explicit[0],'mapping_method':'aliquot UUID/barcode explicitly embedded in official filename'}]
    # Somatic CN/MAF derivations include both input tumor and matched normal in
    # API case relationships. The output belongs to the uniquely linked tumor.
    tumor=[e for e in candidates if e.get('sample','')[13:15] in {'01','02','03','04','05','06','07','08','09','40'}]
    if file.get('data_type') in {'Gene Level Copy Number','Allele-specific Copy Number Segment','Masked Somatic Mutation'} and len(tumor)==1:
        return [{**tumor[0],'mapping_method':'unique tumor in API-linked somatic tumor-normal pair',
                 'matched_normal_aliquot_ids':';'.join(e['aliquot_id'] for e in candidates if e not in tumor)}]
    return candidates


def write_tsv(path, rows, columns=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',encoding='utf-8',newline='') as f:
        w = csv.DictWriter(f,fieldnames=columns,delimiter='\t',extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def manifest(path, rows):
    write_tsv(path,[{'id':r['file_id'],'filename':r['file_name'],'md5':r['md5sum'],
                   'size':r['file_size'],'state':r['state']} for r in rows],
              ['id','filename','md5','size','state'])


def local_path(kind, row):
    return RAW / {'rna':'expression','cn':'copy_number','segments':'copy_number_segments','mutation':'mutation'}[kind] / 'raw_files' / row['file_id'] / row['file_name']

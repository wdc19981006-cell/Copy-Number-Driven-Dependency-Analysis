"""Pin DR46, discover all requested TCGA active/open files, select CN, budget disk."""
from pathlib import Path
import json
import logging
import shutil
import sys
from collections import Counter,defaultdict
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT, setup, now, fetch_text, write_json, hashes
from scripts.utils.gdc import (RAW,API,FILE_FIELDS,PRIORITY,filters,paginated,target_entities,file_row,
                               init_dirs,write_tsv,manifest)

QUERIES = {
 'rna':('Transcriptome Profiling','Gene Expression Quantification','STAR - Counts'),
 'cn':('Copy Number Variation','Gene Level Copy Number',None),
 'segments':('Copy Number Variation','Allele-specific Copy Number Segment',None),
 'mutation':('Simple Nucleotide Variation','Masked Somatic Mutation',None),
}


def select_cn(rows):
    groups = defaultdict(list)
    for row in rows:
        if row['workflow'] not in PRIORITY:
            raise ValueError(f'Unknown CN workflow {row["workflow"]}; explicit policy required')
        if not row.get('sample_id'):
            raise ValueError('CN file lacks uniquely mapped sample')
        groups[row['sample_id']].append(row)
    selected = []
    for sample, candidates in groups.items():
        ordered = sorted(candidates,key=lambda r:(PRIORITY[r['workflow']],r['aliquot'],r['file_id']))
        best = ordered[0]
        tied = [r for r in ordered if r['workflow']==best['workflow']]
        # A sample can have several aliquots. Pick one deterministically without
        # inventing an extra patient; retain all ties and reason in the audit table.
        best = {**best,'selected_workflow':best['workflow'],
                'available_workflows':';'.join(sorted({r['workflow'] for r in candidates},key=lambda x:PRIORITY[x])),
                'same_workflow_candidates':len(tied),
                'selection_rule':'workflow priority; then aliquot submitter_id; then file UUID',
                'candidate_file_ids':';'.join(r['file_id'] for r in ordered)}
        selected.append(best)
    return sorted(selected,key=lambda r:(r['project_id'],r['sample']))


def select_segments(rows, cn):
    groups = defaultdict(list)
    for row in rows:
        if row.get('sample_id'):
            groups[row['sample_id']].append(row)
    selected, audit = [], []
    for gene in cn:
        candidates = groups.get(gene['sample_id'],[])
        exact = [r for r in candidates if r['workflow']==gene['workflow'] and r['aliquot_id']==gene['aliquot_id']]
        same = [r for r in candidates if r['workflow']==gene['workflow']]
        supported = [r for r in candidates if r['workflow'] in PRIORITY]
        pool = exact or same or supported
        if pool:
            best = sorted(pool,key=lambda r:(PRIORITY.get(r['workflow'],99),r['aliquot'],r['file_id']))[0]
            reason = 'same workflow and aliquot' if exact else ('same workflow' if same else 'fallback: ABSOLUTE liftover segments are not published')
            selected.append({**best,'gene_cn_workflow':gene['workflow'],'workflow_match':best['workflow']==gene['workflow'],
                             'segment_selection_reason':reason})
            audit.append({**gene,'segment_file_id':best['file_id'],'segment_workflow':best['workflow'],
                          'segment_selection_reason':reason})
        else:
            audit.append({**gene,'segment_file_id':'','segment_workflow':'','segment_selection_reason':'no current open allele-specific segment for this sample'})
    return sorted({r['file_id']:r for r in selected}.values(),key=lambda r:r['file_id']),audit


def main():
    setup('discover_gdc')
    init_dirs()
    status = json.loads(fetch_text(API+'/status'))
    version = status.get('data_release_version',{})
    if (version.get('major'),version.get('minor')) != (46,0):
        raise ValueError(f'Official API is not DR46.0: {status}; do not relabel another release as DR46')
    path = RAW / 'metadata/gdc_status.json'
    if not path.exists():
        write_json(path,status)
    write_json(RAW / 'metadata/query_context.json',{'queried_at':now(),'api_url':API,
                  'api_version':status.get('version'),'api_tag':status.get('tag'),
                  'data_release':status['data_release'],'active_files':'state=released; current /files endpoint',
                  'program':'TCGA','access':'open'})
    catalogs = {}
    for kind,(category,datatype,workflow) in QUERIES.items():
        cached = RAW / f'metadata/{kind}_files.json'
        if cached.exists():
            files = json.loads(cached.read_text(encoding='utf-8'))
        else:
            files = paginated('files',{'filters':filters(category,datatype,workflow),'fields':FILE_FIELDS},kind)
            write_json(cached,files)
        rows = []
        unmapped = []
        for file in files:
            if file.get('access')!='open' or file.get('state')!='released':
                raise ValueError('Unexpected non-open or inactive result')
            associated = target_entities(file)
            if kind!='mutation' and len(associated)!=1:
                unmapped.append({'file_id':file['file_id'],'associated_count':len(associated),
                                 'associated_entities':file.get('associated_entities',[])})
            for entity in associated or [{}]:
                rows.append(file_row(file,entity))
        if unmapped:
            write_json(RAW / f'metadata/{kind}_mapping_issues.json',unmapped)
            raise ValueError(f'{kind}: {len(unmapped)} files have ambiguous/missing sample mapping; inspect metadata before download')
        catalogs[kind] = rows
        write_tsv(RAW / f'metadata/tcga_{kind}_all_metadata.tsv',rows)
        logging.info('%s files=%s workflows=%s',kind,len(files),dict(Counter(r['workflow'] for r in rows)))
    cn = select_cn(catalogs['cn'])
    segments,audit = select_segments(catalogs['segments'],cn)
    # Current masked MAF workflow should be ensemble, not duplicate per-caller MAFs.
    workflows = Counter(r['workflow'] for r in catalogs['mutation'])
    preferred = [r for r in catalogs['mutation'] if 'Ensemble' in r['workflow'] and 'Mask' in r['workflow']]
    mutation = preferred or catalogs['mutation']
    mutation = list({r['file_id']:r for r in mutation}.values())
    if not preferred and len(workflows)>1:
        raise ValueError(f'Mutation workflow policy requires review: {workflows}')
    selections = {'rna':catalogs['rna'],'cn':cn,'segments':segments,'mutation':mutation}
    names = {'rna':'tcga_star_counts','cn':'tcga_gene_level_cn_selected',
             'segments':'tcga_cn_segments_selected','mutation':'tcga_masked_somatic_mutation'}
    for kind,rows in selections.items():
        unique = list({r['file_id']:r for r in rows}.values())
        manifest(RAW / f'manifests/{names[kind]}_manifest.tsv',unique)
        write_json(RAW / f'metadata/{kind}_selected.json',rows)
        write_tsv(RAW / f'metadata/{names[kind]}_metadata.tsv',rows)
    write_tsv(RAW/'metadata/tcga_star_counts_sample_sheet.tsv',selections['rna'])
    write_tsv(RAW/'metadata/tcga_gene_level_cn_selection.tsv',cn)
    write_tsv(RAW/'metadata/tcga_segment_selection_audit.tsv',audit)
    projects = sorted({r['project_id'] for rows in selections.values() for r in rows if r.get('project_id')})
    estimates = []
    for kind,rows in selections.items():
        files = {r['file_id']:r for r in rows}
        estimates.append({'data_type':kind,'files':len(files),
                          'bytes':sum(r['file_size'] for r in files.values()),
                          'GiB':round(sum(r['file_size'] for r in files.values())/1024**3,3)})
    raw_bytes = sum(e['bytes'] for e in estimates)
    # Budget all STAR metrics in float64 (five matrices), CN and segment/mutation,
    # plus processing headroom. Estimates deliberately include more than download bytes.
    rna_count = len(selections['rna'])
    processed_estimate = rna_count*65000*8*5 + len(cn)*65000*8*2 + 5*1024**3
    # Temporary memory-mapped raw metrics during bounded ETL are additional to
    # final Parquet; include their maximum simultaneous footprint before download.
    working_reserve = rna_count*65000*8*4 + 10*1024**3
    free = shutil.disk_usage(ROOT).free
    required = raw_bytes+processed_estimate+working_reserve
    safe = required <= free*0.8
    estimates += [{'data_type':'processed_conservative_estimate','files':'','bytes':processed_estimate,'GiB':round(processed_estimate/1024**3,3)},
                  {'data_type':'working_reserve','files':'','bytes':working_reserve,'GiB':round(working_reserve/1024**3,3)},
                  {'data_type':'total_required','files':'','bytes':required,'GiB':round(required/1024**3,3)},
                  {'data_type':'free_space','files':'','bytes':free,'GiB':round(free/1024**3,3)}]
    write_tsv(ROOT/'data/manifests/gdc_download_size_estimate.csv',estimates)
    # Use genuine comma-separated CSV for the user-facing size estimate.
    import csv
    for target in [ROOT/'data/manifests/gdc_download_size_estimate.csv', RAW/'manifests/gdc_download_size_estimate.csv']:
        with target.open('w',encoding='utf-8',newline='') as f:
            w=csv.DictWriter(f,fieldnames=['data_type','files','bytes','GiB']);w.writeheader();w.writerows(estimates)
    summary={'created_at':now(),'gdc_release':'46.0','api_version':status.get('version'),
             'api_tag':status.get('tag'),'projects':projects,'project_count':len(projects),
             'file_estimates':estimates,'raw_bytes':raw_bytes,'required_bytes':required,
             'free_bytes':free,'safe_to_download':safe,
             'cn_workflow_distribution':dict(Counter(r['workflow'] for r in cn)),
             'segment_missing_samples':sum(not r['segment_file_id'] for r in audit),
             'segment_workflow_mismatches':sum(bool(r['segment_workflow']) and r['segment_workflow']!=r['workflow'] for r in audit)}
    write_json(RAW/'metadata/discovery_summary.json',summary)
    print(json.dumps(summary,indent=2))
    if not safe:
        raise SystemExit('STOP: required_space > free_space * 0.8; no bulk files downloaded')


if __name__=='__main__':
    main()

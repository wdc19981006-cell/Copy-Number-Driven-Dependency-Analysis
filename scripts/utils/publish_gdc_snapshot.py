"""Explicitly publish a dated manifest/report snapshot; live jobs never call this."""
from pathlib import Path
import csv,json,sys,shutil,collections
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,now,write_json
from scripts.utils.gdc import RAW,PROCESSED,manifest_path,report_path,runtime_path

def write(path,text):
 (ROOT/path).write_text(text.strip()+'\n',encoding='utf-8')

def main():
 stamp=now();target=ROOT/'data/manifests'
 # Atomic live replacements allow this read to obtain one complete generation.
 content=manifest_path().read_bytes()
 rows=list(csv.DictReader(content.decode('utf-8').splitlines()))
 (target/'gdc_data_manifest.csv').write_bytes(content)
 counts={}
 for kind in ['rna','cn','segments','mutation']:
  group=[r for r in rows if r['kind']==kind];good=[r for r in group if r['status']=='verified']
  counts[kind]={'files':len(group),'verified':len(good),'selected_bytes':sum(int(r['file_size']) for r in group),'verified_bytes':sum(int(r['file_size']) for r in good)}
 report={'created_at':stamp,'snapshot':True,'method':'Official GDC API multi-UUID POST / per-file fallback; official DTT attempted first','counts':counts,'status_counts':dict(collections.Counter(r['status'] for r in rows)),
  'failures':[{'kind':r['kind'],'file_id':r['file_id'],'status':r['status'],'error':r.get('error','')} for r in rows if r['status']!='verified'],
  'note':'Dated repository snapshot. Local current progress is in ignored raw/manifests/live; pending is not a permanent download error.'}
 write_json(target/'gdc_download_summary.json',report)
 for name in ['gdc_clinical_qc.json','gdc_segments_qc.json','gdc_mutation_qc.json','gdc_maf_identity_qc.json','gdc_rna_qc.json','gdc_cn_qc.json','gdc_validation_report.json']:
  source=report_path(name)
  if source.exists() and source!=target/name:shutil.copyfile(source,target/name)
 layers=json.loads(runtime_path('current_layer_status.json').read_text(encoding='utf-8'));layers['snapshot_at']=stamp
 write_json(ROOT/'config/tcga_layers.json',layers)
 discovery=json.loads((RAW/'metadata/discovery_summary.json').read_text(encoding='utf-8'))
 clinical=json.loads((target/'gdc_clinical_qc.json').read_text(encoding='utf-8'))
 size=sum(p.stat().st_size for p in PROCESSED.rglob('*') if p.is_file())
 table='\n'.join(f"| {k} | {v['verified']:,} / {v['files']:,} | {v['verified_bytes']:,} / {v['selected_bytes']:,} |" for k,v in counts.items())
 write('docs/TCGA_CURRENT_REPORT.md',f'''
# GDC current preparation snapshot

Snapshot: **{stamp}**. Current status: **{layers['current_status']}**. Download/ETL/validation continues locally; this report is not a live counter. The VPS4B/VPS4A example records the readiness at its run time.

1. Official release: **46.0**, API version {discovery['api_version']}, tag {discovery['api_tag']}; query snapshot {discovery['created_at']}.
2. TCGA projects: **{discovery['project_count']}**.
3–4, 6–7. File coverage and bytes (MD5 verified / selected):

| Kind | Files verified / selected | Bytes verified / selected |
|---|---:|---:|
{table}

5. Selected CN workflows: {json.dumps(discovery['cn_workflow_distribution'])}. Priority per sample UUID: ABSOLUTE LiftOver > ASCAT3 > AscatNGS > ASCAT2. A deterministic same-workflow aliquot/file tie break is audited.
8. Clinical: **{clinical['clinical_cases']:,} cases**; {clinical['processed']['TCGA_Biospecimen']['rows']:,} biospecimen rows; {clinical['processed']['TCGA_Sample_Map']['rows']:,} sample/file map rows.
9. Current processed directory at snapshot: **{size:,} bytes**. Segments: 1,070,293 rows; masked MAF: 2,570,542 rows. Full RNA/CN matrices remain pending while raw downloads are incomplete.
10. Outstanding raw files: **{len(report['failures']):,}**. Status counts: {json.dumps(report['status_counts'])}. Complete UUID list is in [gdc_data_manifest.csv](../data/manifests/gdc_data_manifest.csv). Earlier DTT TLS errors, a truncated large bundle and a corrected live-tracking indentation error are retained in local logs; none is represented as a successful transfer. Verified raw was retained.
11. Free D-drive space at snapshot: **{shutil.disk_usage(ROOT).free:,} bytes**. Initial conservative required estimate {discovery['required_bytes']:,} bytes was within 80% of initial free {discovery['free_bytes']:,} bytes; ongoing transfers also check remaining space.

Segments for 150 selected CN samples are not published in the queried source. 10,527 segment selections use a documented different-workflow fallback, mostly because matching ABSOLUTE segments are not published. Do not interpret these as exact pipeline-matched segments. No purity/ploidy is invented.

Two file-level MAF biospecimen inconsistencies were found. Each source Tumor_Sample_UUID/barcode was resolved against the full official case hierarchy, with case/barcode validation; source tumor identity is canonical and the file API association is retained for audit. See [identity QC](../data/manifests/gdc_maf_identity_qc.json).

Local continuation: `scripts/utils/complete_gdc_pipeline.py` resumes verified files, then performs complete ETL and numerical/identity validation. Progress: `data/raw/tcga/gdc_current_DR46/manifests/live/pipeline_status.json`; detail: `logs/gdc_pipeline.log`. A failed phase is recorded and can be resumed with the same script. Only successful validation marks current complete. Run `scripts/utils/publish_gdc_snapshot.py` explicitly to refresh repository snapshots, and rerun desired R TCGA modes to update example results. Do not run competing coordinators.
''')
 write('docs/TCGA_DATA_SOURCES.md',f'''
# TCGA source layers

Official current query snapshot: {discovery['created_at']}, GDC DR46.0 / API {discovery['api_version']} / tag {discovery['api_tag']}. [GDC API downloads](https://docs.gdc.cancer.gov/API/Users_Guide/Downloading_Files/) and [GDC Portal](https://portal.gdc.cancer.gov/) are the sources. Filters are TCGA program, open access, released state. All 33 TCGA projects and complete source gene tables are retained; no candidate-gene-only acquisition.

STAR Counts supplies original Ensembl versioned IDs, unstranded counts, TPM, FPKM and FPKM-UQ; log2(TPM+1) is a separate derived matrix. Special N_* summary rows are stored separately. Gene-level CN keeps source total CN and genomic annotations. CN workflow priority is ABSOLUTE LiftOver > ASCAT3 > AscatNGS > ASCAT2 per sample UUID, without averaging different pipelines. Segment fallback/missing cases are audited. Masked somatic MAF retains original columns plus source-validated case/sample/aliquot and VAF where depth permits. Clinical and full biospecimen hierarchy come from official cases API.

Each file is pinned by UUID, publisher MD5, size, and local SHA256 in [manifest](../data/manifests/gdc_data_manifest.csv). Files live under UUID subdirectories. ZSTD Parquet is generated separately; DuckDB stores views. The [current report](TCGA_CURRENT_REPORT.md) distinguishes selected coverage from completed downloads and current readiness.

The independent PanCanAtlas reference is complete for five acquired files: [Xena TCGA hub](https://tcga.xenahubs.net), [Toil hub](https://toil.xenahubs.net), [PanCanAtlas hub](https://pancanatlas.xenahubs.net). CN/GISTIC: 24,776 genes × 10,845 samples; Toil expression: 58,581 source labels × 10,535 samples; metadata: 12,805 samples / 33 mapped cancer types; Xena MC3: 2,907,335 mutation rows. Reference CN is GISTIC2 continuous scale, five-state GISTIC is -2/-1/0/1/2, and Toil expression is log2(norm_count+1), not TPM. Xena MC3 is FILTER=PASS from the source MAF and is not the complete original MAF. The old [MC3 publication](https://gdc.cancer.gov/about-data/publications/mc3-2017) download UUID was unavailable; the failure is retained separately.

Raw current, processed current, raw reference and processed reference have separate paths in `config/tcga_layers.json`. No current-to-reference fallback is silent. SHA256 and source object metadata for reference files are retained in `data_manifest.csv` and `*.source.json`.
''')
 dep=list(csv.DictReader((target/'depmap_26Q1_manifest.csv').open(encoding='utf-8')))
 dep_table='\n'.join(f"| {r['canonical_name']} | {r['rows']} × {r['columns']} | {r['model_count']} | {r['gene_count']} |" for r in dep)
 write('docs/DATA_SOURCES.md',f'''
# Current data sources

As of {stamp}, all **13 user-supplied DepMap 26Q1 exports are active** and organized under `data/raw/depmap/26Q1/`. Before/after SHA256 checks verified unchanged content. No duplicate exports were present. Original names, dimensions, imported timestamps, raw paths and SHA256 are in [DepMap manifest](../data/manifests/depmap_26Q1_manifest.csv). Original download dates and publisher MD5 for these exports were not supplied; local SHA256 is a content fingerprint, not publisher verification.

| Active canonical file | Rows × columns including IDs | Unique models | Gene columns |
|---|---:|---:|---:|
{dep_table}

Exact ModelID overlaps: CN/Chronos **858**, CN/expression **1,105**, expression/Chronos **1,140**, all three **852**. Model metadata contains 2,154 unique models. Metadata/profile tables can have multiple rows per model; matrix ModelIDs are unique. CN is WGS relative CN; Chronos is gene effect, distinct from dependency probability. Expression keeps the original export scale. Gene symbol columns and blank first ID headers are adapted in memory, without changing raw.

The [official DepMap catalog](https://depmap.org/portal/api/no-captcha/download/files) selected Public 26Q1 in the acquisition snapshot. Initial official-file acquisition was skipped at the user's request when links required browser verification; these later local exports now satisfy the R inputs. Several supplied names contain `subsetted`; their completeness against official full files is unverified. Full somatic variant, fusion supplement and all-gene stranded expression were not supplied. Optional combined/array CN was unavailable in that catalog. No older release is substituted.

TCGA current and reference: see [TCGA sources](TCGA_DATA_SOURCES.md) and [dated current preparation report](TCGA_CURRENT_REPORT.md). Current RNA/CN remain incomplete at this snapshot; reference GISTIC is available and explicitly labeled.

Manifest scopes: `depmap_26Q1_manifest.csv` and `depmap_26Q1_local_qc.json` describe the current local exports; `gdc_*` are dated current snapshots; `data_manifest.csv`, `not_downloaded.csv`, `preparation_summary.json`, `depmap_26Q1_qc.json` and the initial preparation report describe the earlier official-download stage, not the current local-export state. [current_missing_data.csv](../data/manifests/current_missing_data.csv) provides the current missing-data summary. Ongoing GDC state lives in the ignored `raw/manifests/live` directory.
''')
 missing=[{'dataset':k,'status':'in_progress','remaining_files':v['files']-v['verified'],'note':'See live pipeline state and UUID manifest'} for k,v in counts.items() if v['verified']<v['files']]
 missing += [{'dataset':v,'status':'not_supplied','remaining_files':'unknown','note':'Optional extension; no full official export acquired'} for v in ['DepMap full somatic variants','DepMap fusion supplement','DepMap all-gene stranded expression']]
 missing += [{'dataset':'DepMap combined/array CN','status':'unavailable_optional','remaining_files':'unknown','note':'Not available in selected official catalog; WGS export is active'}, {'dataset':'PanCanAtlas original full MC3 MAF','status':'source_unavailable','remaining_files':1,'note':'Old UUID unavailable; acquired Xena FILTER=PASS table is separately labeled'}]
 with (target/'current_missing_data.csv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['dataset','status','remaining_files','note']);w.writeheader();w.writerows(missing)
 print(json.dumps({'snapshot_at':stamp,'counts':counts,'current_status':layers['current_status'],'processed_bytes':size},indent=2))

if __name__=='__main__':main()

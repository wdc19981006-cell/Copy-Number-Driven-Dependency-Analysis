"""Lossless current GDC clinical/biospecimen preparation; no cohort analysis."""
from pathlib import Path
import json
import logging
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,setup,write_json,now,hashes
from scripts.utils.gdc import RAW,PROCESSED,init_dirs
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def json_value(obj):
    return json.dumps(obj,ensure_ascii=False,separators=(',',':'))


def combine(rows,key):
    values=list(dict.fromkeys(str(r[key]) for r in rows if r.get(key) is not None))
    return ';'.join(values) if values else None


def numeric_values(rows,key):
    values=[]
    for row in rows:
        value=row.get(key)
        if value is not None:
            try: values.append(float(value))
            except (ValueError,TypeError): pass
    return values


def save(name,rows):
    if not rows:
        raise ValueError(f'No rows in {name}')
    table=pa.Table.from_pandas(pd.DataFrame(rows),preserve_index=False)
    path=PROCESSED/f'{name}.parquet'
    temporary=path.with_suffix('.parquet.part')
    pq.write_table(table,temporary,compression='zstd',compression_level=6)
    temporary.replace(path)
    return {'rows':len(rows),'columns':table.num_columns,'bytes':path.stat().st_size,'SHA256':hashes(path)[0]}


def main():
    setup('preprocess_gdc_clinical')
    init_dirs()
    source=RAW/'clinical/tcga_cases.json'
    cases=json.loads(source.read_text(encoding='utf-8'))
    clinical,biospecimen=[],[]
    for case in cases:
        project=case['project']['project_id']
        if not project.startswith('TCGA-'):
            raise ValueError('Non-TCGA clinical case')
        demographic=case.get('demographic',{})
        diagnoses=case.get('diagnoses',[])
        treatments=[t for d in diagnoses for t in d.get('treatments',[])]
        followups=[f for d in diagnoses for f in d.get('follow_ups',[])]
        age_days=numeric_values(diagnoses,'age_at_diagnosis')
        age=demographic.get('age_at_index')
        age_source='demographic.age_at_index (years)'
        if age is None:
            age=age_days[0]/365.25 if len(set(age_days))==1 else None
            age_source='diagnosis.age_at_diagnosis / 365.25' if age is not None else 'unavailable or multiple diagnoses'
        observed=numeric_values(diagnoses+followups,'days_to_last_follow_up')
        death=demographic.get('days_to_death')
        if death is None:
            deaths=numeric_values(diagnoses,'days_to_death')
            death=deaths[0] if len(set(deaths))==1 else None
        clinical.append({'project_id':project,'case_id':case['case_id'],'submitter_id':case.get('submitter_id'),
            'primary_diagnosis':combine(diagnoses,'primary_diagnosis'),'disease_type':case.get('disease_type'),
            'primary_site':case.get('primary_site'),'age':age,'age_source':age_source,
            'age_at_diagnosis_days':age_days[0] if len(set(age_days))==1 else None,
            'sex':demographic.get('gender'),'race':demographic.get('race'),
            'stage':combine(diagnoses,'ajcc_pathologic_stage') or combine(diagnoses,'tumor_stage') or combine(diagnoses,'figo_stage'),
            'grade':combine(diagnoses,'tumor_grade'),'vital_status':demographic.get('vital_status'),
            'days_to_death':death,'days_to_last_follow_up':max(observed) if observed else None,
            'treatments':json_value(treatments),'diagnoses_json':json_value(diagnoses),
            'follow_ups_json':json_value(followups),'demographic_json':json_value(demographic)})
        for sample in case.get('samples',[]):
            base={'project_id':project,'case_id':case['case_id'],'case_submitter_id':case.get('submitter_id'),
                  'sample_id':sample['sample_id'],'sample_submitter_id':sample.get('submitter_id'),
                  'sample_type':sample.get('sample_type'),'tissue_type':sample.get('tissue_type'),
                  'tumor_descriptor':sample.get('tumor_descriptor'),'specimen_type':sample.get('specimen_type'),
                  'preservation_method':sample.get('preservation_method'),
                  'sample_json':json_value({k:v for k,v in sample.items() if k!='portions'})}
            leaves=[]
            for portion in sample.get('portions',[]):
                for analyte in portion.get('analytes',[]):
                    for aliquot in analyte.get('aliquots',[]):
                        leaves.append({**base,'portion_id':portion.get('portion_id'),
                            'analyte_id':analyte.get('analyte_id'),'analyte_type':analyte.get('analyte_type'),
                            'aliquot_id':aliquot['aliquot_id'],'aliquot_submitter_id':aliquot.get('submitter_id'),
                            'aliquot_json':json_value(aliquot)})
            biospecimen.extend(leaves or [{**base,'portion_id':None,'analyte_id':None,'analyte_type':None,
                               'aliquot_id':None,'aliquot_submitter_id':None,'aliquot_json':None}])
    result={'TCGA_Clinical':save('TCGA_Clinical',clinical),
            'TCGA_Biospecimen':save('TCGA_Biospecimen',biospecimen)}
    sample_map=[]
    correction={}
    from scripts.utils.gdc import report_path
    audit=report_path('gdc_maf_identity_qc.json')
    by_aliquot={r['aliquot_id']:r for r in biospecimen if r['aliquot_id']}
    if audit.exists():
        for entry in json.loads(audit.read_text(encoding='utf-8'))['file_metadata_discrepancies']:
            correction.setdefault(entry['file_id'],{})[entry['source_tumor_aliquot']]=by_aliquot[entry['source_tumor_aliquot']]
    for kind in ['rna','cn','segments','mutation']:
        selected=RAW/f'metadata/{kind}_selected.json'
        if not selected.exists():
            continue
        for row in json.loads(selected.read_text(encoding='utf-8')):
            sample_map.append({'DataType':kind,'SampleID':row.get('sample_id'),'SampleBarcode':row.get('sample'),
                              'AliquotID':row.get('aliquot_id'),'AliquotBarcode':row.get('aliquot'),
                              'CaseID':row.get('case_id'),'CaseBarcode':row.get('case'),
                              'ProjectID':row.get('project_id'),'SampleType':row.get('sample_type'),
                              'FileID':row['file_id'],'FileName':row['file_name'],'Workflow':row['workflow'],
                              'MappingMethod':row.get('mapping_method')})
            if kind=='mutation' and row['file_id'] in correction:
                api=sample_map.pop()
                for target in correction[row['file_id']].values():
                    sample_map.append({**api,'API_AliquotID':api['AliquotID'],'API_SampleID':api['SampleID'],
                        'SampleID':target['sample_id'],'SampleBarcode':target['sample_submitter_id'],
                        'AliquotID':target['aliquot_id'],'AliquotBarcode':target['aliquot_submitter_id'],
                        'SampleType':target['sample_type'],'MappingMethod':'Source MAF tumor UUID/barcode resolved against full official case biospecimen; original API association retained'})
    if sample_map:
        result['TCGA_Sample_Map']=save('TCGA_Sample_Map',sample_map)
    write_json(ROOT/'data/manifests/gdc_clinical_qc.json',{'created_at':now(),'source_sha256':hashes(source)[0],
               'clinical_cases':len(clinical),'project_count':len({r['project_id'] for r in clinical}),
               'sample_ids':len({r['sample_id'] for r in biospecimen}),
               'aliquot_ids':len({r['aliquot_id'] for r in biospecimen if r['aliquot_id']}),
               'processed':result})
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()

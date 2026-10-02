"""Fetch all current TCGA cases including clinical and complete biospecimen tree."""
from pathlib import Path
import sys
import json
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import setup,write_json,fetch_text,now
from scripts.utils.gdc import API,RAW,init_dirs,paginated,eq


def main():
    setup('download_gdc_clinical')
    init_dirs()
    status=json.loads(fetch_text(API+'/status'))
    if status.get('data_release_version',{}).get('major')!=46:
        raise ValueError('GDC release changed before clinical query')
    rows=paginated('cases',{'filters':eq('project.program.name','TCGA'),
          'expand':'project,demographic,diagnoses,diagnoses.treatments,diagnoses.follow_ups,samples,samples.portions,samples.portions.analytes,samples.portions.analytes.aliquots'},
          'tcga_cases',page_size=250)
    write_json(RAW/'clinical/tcga_cases.json',rows)
    write_json(RAW/'clinical/clinical_query_summary.json',{'queried_at':now(),'release':status['data_release'],
               'cases':len(rows),'projects':sorted({r['project']['project_id'] for r in rows})})
    print(f'Complete TCGA current case metadata: {len(rows)} cases')


if __name__=='__main__':
    main()

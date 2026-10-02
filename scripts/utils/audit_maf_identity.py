from pathlib import Path
import json,gzip,csv,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT
from scripts.utils.gdc import RAW,local_path
rows=json.loads((RAW/'metadata/mutation_selected.json').read_text())
for n,row in enumerate(rows):
 if not 1700<=n<2200:continue
 file=local_path('mutation',row)
 with gzip.open(file,'rt',encoding='utf-8-sig') as f:
  reader=csv.DictReader((line for line in f if not line.startswith('#')),delimiter='\t')
  source={}
  for v in reader:source[(v.get('Tumor_Sample_UUID'),v.get('Tumor_Sample_Barcode'))]=source.get((v.get('Tumor_Sample_UUID'),v.get('Tumor_Sample_Barcode')),0)+1
 if any(key[0]!=row.get('aliquot_id') for key in source):
  report={'number':n,'selected':row,'source_ids':[{ 'uuid':k[0],'barcode':k[1],'rows':v} for k,v in source.items()]}
  (ROOT/'docs/MAF_IDENTITY_AUDIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2));break

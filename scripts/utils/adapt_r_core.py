"""Split the supplied R source at section markers; preserve statistical functions."""
from pathlib import Path
import re, hashlib, json, csv
ROOT=Path(__file__).resolve().parents[2]
SOURCE=Path('C:/Users/22241/.codex/attachments/c293b9c9-85dd-4b79-8456-9adefc59465d/已粘贴的文本.txt')
text=SOURCE.read_text(encoding='utf-8-sig')
target=ROOT/'scripts/R';target.mkdir(parents=True,exist_ok=True)
(ROOT/'docs/USER_SUPPLIED_ANALYSIS.R').write_bytes(SOURCE.read_bytes())
marks=list(re.finditer(r'# =+\n# (\d+)\.',text))
sections={int(m[1]):text[m.start():marks[i+1].start() if i+1<len(marks) else len(text)] for i,m in enumerate(marks)}
for file,nums in [('data_access.R',[3,4,5,6]),('plotting.R',[7,8]),('dependency_screen.R',[9]),('lineage_analysis.R',[10,11])]:
 content='\n'.join(sections[n] for n in nums)
 # R fread assigns V1 to a blank CSV column; only the in-memory identifier changes.
 if file=='data_access.R':content=content.replace('"Unnamed: 0",','"Unnamed: 0",\n    "V1",')
 (target/file).write_text(content,encoding='utf-8')
 provenance={'input_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
 'source':'docs/USER_SUPPLIED_ANALYSIS.R','sections':{f:ns for f,ns in [('data_access.R',[3,4,5,6]),('plotting.R',[7,8]),('dependency_screen.R',[9]),('lineage_analysis.R',[10,11])]},
 'statistical_core':'Function bodies retained; only blank ID column recognition added. Adapters and added modes are separate.'}
 (ROOT/'docs/R_CORE_PROVENANCE.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
 # Fix the coverage manifest for this export format: headers contain gene symbols.
 p=ROOT/'data/manifests/depmap_26Q1_manifest.csv'
 with p.open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f));fields=list(rows[0])
 for row in rows:
  if row['canonical_name'].startswith(('CopyNumber_','CRISPR_Chronos','CRISPR_GeneDependency','Expression_','Mutation_')):
   row['gene_count']=str(int(row['columns'])-1)
   row['notes']='User-supplied 26Q1 portal export; gene-symbol columns. Analyze all supplied columns; official full-release completeness not asserted.'
 with p.open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
 p=ROOT/'docs/DEPMAP_INSPECTION.json';d=json.loads(p.read_text(encoding='utf-8'));d['files']=rows;p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')

if __name__=='__main__':print('Statistical core extracted; source preserved with SHA256.')

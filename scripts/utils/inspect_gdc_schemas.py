"""Print only a few source lines for parser development, never entire matrices."""
from pathlib import Path
import gzip
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.download.download_gdc import selected
from scripts.utils.gdc import local_path

for kind in ['rna','cn','segments','mutation']:
    row=selected(kind)[0]
    path=local_path(kind,row)
    print(kind,row['workflow'],path.name)
    opener=gzip.open if path.name.endswith('.gz') else open
    with opener(path,'rt',encoding='utf-8-sig') as f:
        for number,line in enumerate(f):
            if line.startswith('#'):
                print(line[:150].strip())
                continue
            print(line[:1600].strip())
            if number>7:
                break

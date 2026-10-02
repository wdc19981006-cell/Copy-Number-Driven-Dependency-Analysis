"""Seed ignored local runtime state; keep version-controlled snapshots stable."""
from pathlib import Path
import sys,shutil
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT
from scripts.utils.gdc import runtime_path
for source in (ROOT/'data/manifests').glob('gdc_*.json'):
 target=runtime_path(source.name)
 if not target.exists():shutil.copy2(source,target)
for name in ['gdc_data_manifest.csv']:
 source=ROOT/'data/manifests'/name;target=runtime_path(name)
 if not target.exists():shutil.copy2(source,target)
target=runtime_path('current_layer_status.json')
if not target.exists():shutil.copy2(ROOT/'config/tcga_layers.json',target)
print('Runtime tracking seeded; publish snapshots explicitly after verification.')

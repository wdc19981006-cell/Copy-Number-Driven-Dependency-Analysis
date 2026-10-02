"""Move authorized reference files within the project, verifying contents before/after."""
from pathlib import Path
import json
import os
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT, setup, load_manifest, record, hashes, write_json, now


def checked(path):
    path = path.resolve()
    if not path.is_relative_to(ROOT.resolve()):
        raise ValueError('Move target is outside project workspace')
    return path


def main():
    setup('organize_tcga_reference')
    destination = checked(ROOT / 'data/raw/tcga/pancanatlas_reference')
    destination.mkdir(parents=True,exist_ok=True)
    audit = []
    for row in load_manifest():
        if row['database']!='TCGA' or row['release']!='PANCAN':
            continue
        source = checked(ROOT / row['relative_path']) if row['relative_path'] else None
        target = checked(destination / row['filename'])
        if row['status']=='verified':
            if source and source.exists() and source!=target:
                if hashes(source)[0]!=row['SHA256']:
                    raise ValueError(f'Raw SHA256 failed before moving {source}')
                if target.exists():
                    raise FileExistsError(f'Destination exists; refusing overwrite: {target}')
                source.rename(target)
            if not target.exists() or hashes(target)[0]!=row['SHA256']:
                raise ValueError('Reference relocation content verification failed')
        record({**row,'relative_path':target.relative_to(ROOT).as_posix()})
        audit.append({'file':row['filename'],'old_path':row['relative_path'],
                      'new_path':target.relative_to(ROOT).as_posix(),'SHA256':row['SHA256'],'status':row['status']})
    old = checked(ROOT/'data/processed/tcga')
    new = checked(old/'pancanatlas_reference')
    new.mkdir(parents=True,exist_ok=True)
    for path in list(old.iterdir()):
        if not path.is_file() or path.suffix not in ('.parquet','.json','.duckdb'):
            continue
        target = checked(new/path.name)
        if target.exists():
            raise FileExistsError('Processed destination already exists')
        checked(path).rename(target)
    for provenance in new.glob('*.provenance.json'):
        meta = json.loads(provenance.read_text(encoding='utf-8'))
        meta['source'] = meta['source'].replace('data/raw/tcga/pancan/','data/raw/tcga/pancanatlas_reference/')
        meta['reference_relocated_at']=now()
        write_json(provenance,meta)
    # Rebuild only the tiny DuckDB view file to point at relocated Parquet files.
    from scripts.preprocess.preprocess_tcga import create_duckdb
    create_duckdb(new)
    write_json(ROOT/'config/tcga_layers.json',{'default_layer':'gdc_DR46',
               'gdc_release':'46.0','current_status':'preparing',
               'reference_raw':'data/raw/tcga/pancanatlas_reference',
               'reference_processed':'data/processed/tcga/pancanatlas_reference',
               'current_raw':'data/raw/tcga/gdc_current_DR46',
               'current_processed':'data/processed/tcga/gdc_DR46'})
    write_json(ROOT/'data/manifests/reference_relocation.json',{'relocated_at':now(),'files':audit,'hashes_preserved':True})
    print('Five reference raw files relocated, SHA256 unchanged; reference Parquet and views preserved.')


if __name__=='__main__':
    main()

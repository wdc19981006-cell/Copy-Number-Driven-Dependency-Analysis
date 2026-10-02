"""Record the user's explicit decision to skip DepMap files; performs no network IO."""
from pathlib import Path
import csv
import json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT, CATALOG_URL, latest_release, setup, write_json, record, now


def main():
    setup('depmap_selection')
    with (ROOT / 'data/manifests/depmap_catalog.csv').open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    release = latest_release(rows).removeprefix('DepMap Public ')
    folder = ROOT / 'data/raw/depmap' / release
    folder.mkdir(parents=True, exist_ok=True)
    (ROOT / 'data/processed/depmap' / release).mkdir(parents=True, exist_ok=True)
    official = {r['filename']: r for r in rows if r['release'] == f'DepMap Public {release}'}
    datasets = json.loads((ROOT / 'config/datasets.json').read_text(encoding='utf-8'))['depmap']
    for ds in datasets:
        filename = ds['filename']
        r = official.get(filename)
        record(dict(database='DepMap', release=release, dataset=filename, filename=filename,
                    original_url=r.get('url', '') if r else '',
                    publisher_md5=r['md5_hash'] if r else '',
                    relative_path=(folder / filename).relative_to(ROOT).as_posix(),
                    status='skipped_by_user' if r else ('unavailable_optional' if ds.get('optional') else 'missing_required'),
                    error='User explicitly requested skipping DepMap downloads after official endpoint required browser verification.'
                          if r else 'Not supplied in selected release; no older release substituted.'))
    write_json(ROOT / 'config/current_release.json', {
        'depmap_release': release, 'download_date': None, 'catalog_checked_at': now(),
        'source': 'DepMap Portal', 'catalog_url': CATALOG_URL,
        'installation_status': 'skipped_by_user',
        'note': 'Release discovered, but no DepMap data files have been downloaded. download_date remains null.'})
    print(f'Discovered {release}; 16 supplied requested files skipped, optional OmicsCNGene.csv unavailable.')


if __name__ == '__main__':
    main()

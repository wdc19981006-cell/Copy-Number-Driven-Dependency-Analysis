from pathlib import Path
import argparse
import json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT, setup, download, fetch_text, write_json, now


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--only', nargs='*', help='Optional filenames to download; default is all configured datasets.')
    args = parser.parse_args()
    setup('download_tcga')
    datasets = json.loads((ROOT / 'config/datasets.json').read_text(encoding='utf-8'))['tcga']
    failures = []
    for ds in datasets:
        if args.only and ds['filename'] not in args.only:
            continue
        path = download('TCGA', 'PANCAN', ds['dataset'], ds['filename'], ds['url'], ds.get('md5', ''))
        if path is None:
            failures.append(ds['filename'])
        if 'xenahubs.net/download/' in ds['url']:
            metadata_url = ds['url'].removesuffix('.gz') + '.json'
            meta = ROOT / 'data/manifests' / (ds['filename'] + '.source.json')
            if not meta.exists():
                try:
                    text = fetch_text(metadata_url)
                    value = json.loads(text)
                    write_json(meta, {'url': metadata_url, 'retrieved_at': now(), 'metadata': value})
                except Exception as e:
                    import logging
                    logging.warning('Source dictionary unavailable for %s: %s', ds['filename'], e)
    print(json.dumps({'failed': failures}, indent=2))
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())

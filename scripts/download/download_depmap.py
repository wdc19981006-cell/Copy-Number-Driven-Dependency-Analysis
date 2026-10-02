from pathlib import Path
import argparse
import csv
import json
import logging
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import (ROOT, CATALOG_URL, setup, get_catalog, latest_release,
                                  now, write_json, download, record)


def main():
    parser = argparse.ArgumentParser(description='Download only requested complete analysis-level DepMap files.')
    parser.add_argument('--url-catalog', type=Path, help='Official CSV exported after browser verification (release,filename,url,md5_hash).')
    parser.add_argument('--release', help='Explicitly authorize selection of a release, e.g. 26Q3.')
    args = parser.parse_args()
    setup('download_depmap')
    rows = get_catalog()
    available = latest_release(rows).removeprefix('DepMap Public ')
    config = ROOT / 'config/current_release.json'
    local = json.loads(config.read_text(encoding='utf-8')) if config.exists() else None
    release = args.release or (local['depmap_release'] if local else available)
    if not any(r['release'] == f'DepMap Public {release}' for r in rows):
        raise ValueError(f'Release {release} is not in the official catalog.')
    if local and available != local['depmap_release'] and not args.release:
        logging.info('New release %s available; retaining selected local %s. Use --release only after approval.', available, release)
    if not local or args.release:
        write_json(config, {'depmap_release': release, 'download_date': now(),
                            'source': 'DepMap Portal', 'catalog_url': CATALOG_URL,
                            'installation_status': 'pending'})
    urls = {}
    if args.url_catalog:
        with args.url_catalog.open(encoding='utf-8-sig', newline='') as f:
            for r in csv.DictReader(f):
                if r.get('release') == f'DepMap Public {release}':
                    urls[r['filename']] = r
    official = {r['filename']: r for r in rows if r['release'] == f'DepMap Public {release}'}
    datasets = json.loads((ROOT / 'config/datasets.json').read_text(encoding='utf-8'))['depmap']
    counts = {'verified': 0, 'unavailable_optional': 0, 'blocked_or_failed': 0}
    for ds in datasets:
        filename = ds['filename']
        if filename not in official:
            state = 'unavailable_optional' if ds.get('optional') else 'missing_required'
            record(dict(database='DepMap', release=release, dataset=filename, filename=filename,
                        status=state, error='Not provided in this release; no other release substituted.'))
            counts['unavailable_optional' if ds.get('optional') else 'blocked_or_failed'] += 1
            continue
        r = official[filename]
        supplied = urls.get(filename, {})
        if supplied.get('md5_hash') and supplied['md5_hash'].lower() != r['md5_hash'].lower():
            raise ValueError(f'URL catalog checksum does not match selected release: {filename}')
        path = download('DepMap', release, filename, filename,
                        supplied.get('url') or r.get('url', ''), r['md5_hash'])
        counts['verified' if path else 'blocked_or_failed'] += 1
    selected = json.loads(config.read_text(encoding='utf-8'))
    selected.update(installation_status='complete' if not counts['blocked_or_failed'] else 'incomplete',
                    last_attempt=now(), file_status=counts)
    write_json(config, selected)
    write_json(ROOT / 'data/manifests/depmap_download_summary.json', selected)
    print(json.dumps(selected, indent=2))
    return 1 if counts['blocked_or_failed'] else 0


if __name__ == '__main__':
    raise SystemExit(main())

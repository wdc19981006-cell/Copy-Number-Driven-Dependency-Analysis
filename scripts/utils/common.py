from __future__ import annotations

import csv
import hashlib
import io
import json
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone, timedelta

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / '.runtime'))
CATALOG_URL = 'https://depmap.org/portal/api/no-captcha/download/files'
MANIFEST = ROOT / 'data/manifests/data_manifest.csv'
TCGA_REFERENCE_RAW = ROOT / 'data/raw/tcga/pancanatlas_reference'
TCGA_REFERENCE_PROCESSED = ROOT / 'data/processed/tcga/pancanatlas_reference'
FIELDS = ['database', 'release', 'dataset', 'filename', 'original_url',
          'download_date', 'file_size', 'SHA256', 'status', 'relative_path',
          'publisher_md5', 'MD5', 'etag', 'last_modified', 'error']


def now():
    return datetime.now(timezone(timedelta(hours=8))).isoformat(timespec='seconds')


def setup(name):
    for directory in ['data/raw/depmap', 'data/raw/tcga/pancanatlas_reference',
                      'data/processed/depmap', 'data/processed/tcga',
                      'data/manifests', 'config', 'results', 'logs', 'docs']:
        (ROOT / directory).mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',
                        handlers=[logging.FileHandler(ROOT / f'logs/{name}.log', encoding='utf-8'),
                                  logging.StreamHandler(sys.stdout)], force=True)


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    os.replace(temp, path)


def hashes(path):
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            sha.update(block)
            md5.update(block)
    return sha.hexdigest(), md5.hexdigest()


def curl_binary():
    preferred = Path('C:/Program Files/Git/mingw64/bin/curl.exe')
    return str(preferred) if preferred.exists() else (shutil.which('curl') or 'curl')


def curl(args, **kwargs):
    # Argument array prevents shell expansion and preserves Windows paths.
    return subprocess.run([curl_binary(), '--fail', '--location', '--silent',
                           '--show-error', '--connect-timeout', '30', *map(str, args)],
                          check=True, **kwargs)


def fetch_text(url):
    result = curl(['--max-time', '120', '--retry', '2', url], capture_output=True)
    return result.stdout.decode('utf-8-sig')


def load_manifest():
    if not MANIFEST.exists():
        return []
    with MANIFEST.open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def record(item):
    rows = load_manifest()
    key = tuple(item.get(k, '') for k in ['database', 'release', 'filename'])
    rows = [r for r in rows if tuple(r[k] for k in ['database', 'release', 'filename']) != key]
    rows.append({k: item.get(k, '') for k in FIELDS})
    tmp = MANIFEST.with_suffix('.csv.tmp')
    with tmp.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: (r['database'], r['release'], r['filename'])))
    os.replace(tmp, MANIFEST)


def remote_info(url):
    result = curl(['--head', '--max-time', '90', '--retry', '2', url], capture_output=True)
    text = result.stdout.decode('latin1')
    # Only final response headers matter after redirects.
    block = re.split(r'(?m)^HTTP/\S+ ', text)[-1]
    info = {}
    for line in block.splitlines()[1:]:
        if ':' in line:
            k, v = line.split(':', 1)
            info[k.strip().lower()] = v.strip()
    return info


def download(database, release, dataset, filename, url, publisher_md5=''):
    folder = ROOT / 'data/raw' / ('depmap/' + release if database == 'DepMap' else 'tcga/pancanatlas_reference')
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / filename
    base = dict(database=database, release=release, dataset=dataset, filename=filename,
                original_url=url, relative_path=path.relative_to(ROOT).as_posix(),
                publisher_md5=publisher_md5)
    prior = next((r for r in load_manifest() if r['relative_path'] == base['relative_path']), {})
    try:
        if path.exists():
            sha, md5 = hashes(path)
            if path.stat().st_size == 0:
                raise ValueError('Existing raw file is empty; it will not be overwritten.')
            if publisher_md5 and md5.lower() != publisher_md5.lower():
                raise ValueError('Existing raw file differs from publisher MD5; preserved unchanged.')
            if prior.get('SHA256'):
                if sha != prior['SHA256']:
                    raise ValueError('Raw SHA256 changed; refusing to overwrite.')
                if prior.get('status') == 'verified':
                    logging.info('SHA256 verified; skip %s', filename)
                    return path
            elif not publisher_md5:
                raise ValueError('Untracked raw file: cannot establish download completeness; preserved unchanged.')
            # A user supplied official raw file is acceptable only with publisher verification.
            record({**prior, **base, 'download_date': prior.get('download_date') or now(),
                    'file_size': path.stat().st_size, 'SHA256': sha, 'MD5': md5, 'status': 'verified'})
            return path
        if not url:
            record({**base, 'status': 'blocked_verification',
                    'error': 'Official no-captcha catalog omits download URLs. Supply browser-verified official catalog.'})
            return None
        headers = remote_info(url)
        expected = int(headers['content-length']) if headers.get('content-length') else None
        part = path.with_name(path.name + '.part')
        sidecar = part.with_name(part.name + '.json')
        if part.exists() and part.stat().st_size:
            if not sidecar.exists():
                raise ValueError('Partial file has no remote identity; preserve and inspect before resuming.')
            old = json.loads(sidecar.read_text(encoding='utf-8'))
            if old != {'url': url, 'etag': headers.get('etag', ''), 'length': expected}:
                raise ValueError('Remote identity changed; partial file preserved, refusing unsafe resume.')
        write_json(sidecar, {'url': url, 'etag': headers.get('etag', ''), 'length': expected})
        base.update(etag=headers.get('etag', ''), last_modified=headers.get('last-modified', ''))
        record({**base, 'status': 'downloading'})
        logging.info('Download %s (%s bytes)', filename, expected)
        if expected is None or not part.exists() or part.stat().st_size != expected:
            curl(['--retry', '5', '--retry-delay', '3', '--continue-at', '-',
                  '--output', part, url])
        size = part.stat().st_size
        if size == 0 or (expected is not None and size != expected):
            raise ValueError(f'Download size mismatch: {size} vs {expected}')
        with part.open('rb') as f:
            prefix = f.read(300).lstrip().lower()
        if prefix.startswith((b'<!doctype', b'<html', b'<?xml')):
            raise ValueError('Server returned HTML/XML instead of data.')
        sha, md5 = hashes(part)
        if publisher_md5 and md5.lower() != publisher_md5.lower():
            raise ValueError('Publisher MD5 mismatch; partial retained for inspection.')
        part.rename(path)  # Windows refuses an existing target: raw is never overwritten.
        sidecar.unlink(missing_ok=True)
        record({**base, 'download_date': now(), 'file_size': size, 'SHA256': sha,
                'MD5': md5, 'status': 'verified'})
        logging.info('Verified %s SHA256=%s', filename, sha)
        return path
    except Exception as e:
        logging.error('%s: %s', filename, e)
        record({**prior, **base, 'status': 'failed', 'error': str(e)})
        return None


def latest_release(rows):
    choices = {}
    for r in rows:
        match = re.fullmatch(r'DepMap Public (\d{2})Q([1-4])(?: v(\d+))?', r['release'])
        if match:
            choices[(int(match[1]), int(match[2]), int(match[3] or 1))] = r['release']
    if not choices:
        raise ValueError('Official catalog contains no recognizable DepMap Public release.')
    return choices[max(choices)]


def get_catalog():
    content = fetch_text(CATALOG_URL)
    rows = list(csv.DictReader(io.StringIO(content)))
    if not rows or not {'release', 'filename', 'md5_hash'}.issubset(rows[0]):
        raise ValueError('Official catalog is not a recognized CSV.')
    snapshot = ROOT / 'data/manifests' / ('depmap_catalog_' + datetime.now().strftime('%Y%m%d_%H%M%S') + '.csv')
    snapshot.write_text(content, encoding='utf-8', newline='')
    (ROOT / 'data/manifests/depmap_catalog.csv').write_text(content, encoding='utf-8', newline='')
    return rows


def current_release():
    path = ROOT / 'config/current_release.json'
    if not path.exists():
        raise FileNotFoundError('Run download_depmap.py to discover and select a release first.')
    return json.loads(path.read_text(encoding='utf-8'))['depmap_release']


def verified_raw(database, release, filename):
    row = next((r for r in load_manifest() if (r['database'], r['release'], r['filename']) ==
                (database, release, filename)), None)
    if not row or row['status'] != 'verified':
        return None
    path = ROOT / row['relative_path']
    if not path.exists() or hashes(path)[0] != row['SHA256']:
        raise ValueError(f'Raw SHA256 verification failed: {path}')
    return path, row

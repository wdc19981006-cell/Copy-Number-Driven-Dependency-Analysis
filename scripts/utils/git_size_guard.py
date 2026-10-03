"""Reject data/runtime/secrets and files >50 MB before staging and committing."""
import argparse
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LIMIT = 50_000_000
FORBIDDEN_PREFIXES = ('data/raw/', 'data/processed/', 'data/archive/', '.runtime/',
                      'renv/library/', 'renv/staging/', 'logs/', 'tools/')
FORBIDDEN_SUFFIXES = {'.parquet', '.duckdb', '.rds', '.rdata', '.exe', '.zip',
                      '.feather', '.h5', '.h5ad', '.mtx', '.tmp', '.part', '.pyc'}
SECRET = re.compile(rb'gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|'
                    rb'sk-[A-Za-z0-9]{30,}|-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----')


def git(*args, cwd=ROOT):
    return subprocess.check_output(['git', '-c', f'safe.directory={Path(cwd).as_posix()}',
                                    *args], cwd=cwd)


def forbidden(name):
    path = Path(name)
    return (name.startswith(FORBIDDEN_PREFIXES) or
            bool({'raw', 'processed', 'cache', '__pycache__'} & {p.lower() for p in path.parts}) or
            path.suffix.lower() in FORBIDDEN_SUFFIXES or
            name.endswith(('.duckdb.wal', '.mtx.gz')))


def inspect_worktree(paths, root=ROOT):
    issues = []
    for name in paths:
        target = (root / name).resolve()
        if not target.is_relative_to(root.resolve()):
            issues.append((name, 'outside project'))
            continue
        candidates = list(target.rglob('*')) if target.is_dir() else [target]
        for p in candidates:
            if not p.is_file():
                continue
            rel = p.relative_to(root).as_posix()
            if forbidden(rel):
                issues.append((rel, 'excluded data/runtime/cache'))
            elif p.stat().st_size > LIMIT:
                issues.append((rel, 'larger than 50 MB; upload a top table/manifest instead'))
            elif p.suffix.lower() in {'.py', '.r', '.json', '.yaml', '.yml', '.md', '.txt'} and SECRET.search(p.read_bytes()):
                issues.append((rel, 'credential pattern'))
    return issues


def inspect_index(root=ROOT):
    entries = git('ls-files', '-s', '-z', cwd=root).split(b'\0')
    issues, files = [], []
    parsed = [entry.split(b'\t', 1) for entry in entries if entry]
    objects = [meta.split()[1] for meta, _ in parsed]
    sizes = subprocess.run(['git', '-c', f'safe.directory={Path(root).as_posix()}',
                            'cat-file', '--batch-check=%(objectsize)'], cwd=root,
                           input=b'\n'.join(objects) + b'\n', stdout=subprocess.PIPE,
                           check=True).stdout.splitlines()
    for (meta, name), size_bytes in zip(parsed, sizes):
        _, sha, stage = meta.split()
        name = name.decode('utf-8')
        if stage != b'0':
            issues.append((name, 'unresolved merge'))
            continue
        size = int(size_bytes)
        files.append((size, name))
        if forbidden(name):
            issues.append((name, 'excluded data/runtime/cache'))
        if size > LIMIT:
            issues.append((name, 'larger than 50 MB'))
        if Path(name).suffix.lower() in {'.py', '.r', '.json', '.yaml', '.yml', '.md', '.txt'}:
            if SECRET.search(git('cat-file', 'blob', sha.decode(), cwd=root)):
                issues.append((name, 'credential pattern'))
    return {'indexed_files': len(files), 'largest_files': sorted(files, reverse=True)[:5], 'issues': issues}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=['stage', 'audit'])
    p.add_argument('paths', nargs='*')
    a = p.parse_args()
    if a.action == 'stage':
        if not a.paths:
            p.error('stage requires explicit paths')
        # Ignore already ignored children of directories (e.g. __pycache__),
        # while still rejecting an explicitly named forbidden file.
        eligible = git('ls-files', '--cached', '--others', '--exclude-standard', '-z',
                       '--', *a.paths).decode('utf-8').split('\0')
        explicit_files = [name for name in a.paths if not (ROOT / name).is_dir()]
        issues = inspect_worktree(sorted(set(filter(None, eligible + explicit_files))))
        if issues:
            print(json.dumps({'issues': issues}, ensure_ascii=False, indent=2))
            raise SystemExit('Stage refused; no files added')
        git('add', '--', *a.paths)
    report = inspect_index()
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if report['issues']:
        raise SystemExit('Index guard failed; commit/push forbidden')


if __name__ == '__main__':
    main()

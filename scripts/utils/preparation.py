"""Provenance checks for derived files. Raw inputs are read only."""
from pathlib import Path
import json
import os
import pyarrow.parquet as pq
from scripts.utils.common import ROOT, hashes, now, write_json


def output_ready(path, source_sha, version):
    path = Path(path)
    meta = path.with_suffix(path.suffix + '.provenance.json')
    if not path.exists() or not meta.exists():
        return False
    m = json.loads(meta.read_text(encoding='utf-8'))
    if m.get('source_sha256') != source_sha or m.get('processor_version') != version:
        return False
    return hashes(path)[0] == m.get('output_sha256')


def finish(temp, path, source, source_sha, version, extra=None):
    path = Path(path)
    # Only derived artifacts are replaceable, never raw files.
    os.replace(temp, path)
    write_json(path.with_suffix(path.suffix + '.provenance.json'), {
        'source': Path(source).relative_to(ROOT).as_posix(), 'source_sha256': source_sha,
        'output_sha256': hashes(path)[0], 'created_at': now(),
        'processor_version': version, **(extra or {})})


def parquet_qc(path):
    file = pq.ParquetFile(path)
    return {'rows': file.metadata.num_rows, 'columns': len(file.schema_arrow.names),
            'row_groups': file.metadata.num_row_groups,
            'compression': sorted({file.metadata.row_group(g).column(c).compression
                                   for g in range(min(2, file.metadata.num_row_groups))
                                   for c in range(file.metadata.num_columns)})}

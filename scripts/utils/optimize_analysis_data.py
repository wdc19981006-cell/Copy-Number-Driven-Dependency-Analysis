"""Offline, gated preparation/cleanup for the two public DR46/26Q1 workflows.

No download or statistical code lives here. Existing matrices stay byte-identical.
Cleanup requires independently produced regression and isolated-workflow evidence.
"""
import argparse
import csv
import hashlib
import json
import os
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.utils.common import hashes, now, write_json
from scripts.utils.tcga_current import (BASELINE_METHOD, sample_baselines, select_rna,
                                        prepared_table)
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

TCGA = ROOT / 'data/processed/tcga/gdc_DR46'
DEPMAP = ROOT / 'data/processed/depmap/26Q1'
WORK = ROOT / '.runtime/database_optimization'
MANIFESTS = ROOT / 'data/manifests'
REPORT = MANIFESTS / 'analysis_ready_validation.json'
CORE = [TCGA / (n + '.parquet') for n in (
    'TCGA_GeneLevel_CN', 'TCGA_GeneLevel_CN.genes', 'TCGA_STAR_TPM',
    'TCGA_STAR_TPM.genes', 'TCGA_Sample_Baselines', 'TCGA_Sample_Metadata',
    'TCGA_RNA_Selection')] + [DEPMAP / (n + '.parquet') for n in (
        'OmicsCNGeneWGS', 'CRISPRGeneEffect', 'Model')]
CORE_RESOLVED = {p.resolve() for p in CORE}
TUMORS = ['Primary Tumor', 'Recurrent Tumor', 'Metastatic', 'Additional Metastatic',
          'Additional - New Primary', 'Primary Blood Derived Cancer - Peripheral Blood',
          'Primary Blood Derived Cancer - Bone Marrow']


def signature(path):
    stat = path.stat()
    return dict(bytes=stat.st_size, mtime_ns=stat.st_mtime_ns, SHA256=hashes(path)[0])


def persist(name, frame, inputs, **extra):
    path = TCGA / f'{name}.parquet'
    if path.exists():
        raise FileExistsError(f'Refusing to overwrite prepared database: {path}')
    frame.to_parquet(path, compression='zstd', compression_level=6, index=False)
    write_json(path.with_suffix('.provenance.json'), dict(
        release='46.0', algorithm_version='analysis-ready-v1', created_at=now(),
        SHA256=hashes(path)[0], inputs=inputs, rows=len(frame), **extra))
    actual, _ = prepared_table(name)
    pd.testing.assert_frame_equal(actual, frame, check_exact=True)
    return path


def prepare():
    if REPORT.exists():
        raise FileExistsError('Preparation report already exists; inspect before resuming')
    WORK.mkdir(parents=True, exist_ok=True)
    before = json.loads((WORK / 'before.json').read_text())
    inputs = {}
    coverage = {}
    for path in CORE[:4] + CORE[7:]:
        identity = signature(path)
        side = path.with_suffix('.provenance.json')
        if not side.exists():
            side = Path(str(path) + '.provenance.json')
        if side.exists():
            meta = json.loads(side.read_text())
            if meta.get('SHA256') and meta['SHA256'] != identity['SHA256']:
                raise ValueError(f'Original processed SHA256 mismatch: {path.name}')
        pf = pq.ParquetFile(path)
        codecs = Counter(pf.metadata.row_group(r).column(c).compression
                         for r in range(pf.num_row_groups)
                         for c in range(pf.metadata.num_columns))
        if set(codecs) != {'ZSTD'}:
            raise ValueError(f'Unexpected codec: {path.name}: {codecs}')
        inputs[path.name] = identity
        id_col = 'SampleID' if path.parent == TCGA else 'ModelID'
        item = dict(**identity, rows=pf.metadata.num_rows, columns=len(pf.schema_arrow),
                    codecs=dict(codecs), numerical_types=dict(Counter(str(f.type) for f in pf.schema_arrow)),
                    optimization='retain existing ZSTD bytes; no re-encoding or precision change')
        if id_col in pf.schema_arrow.names:
            ids = pf.read(columns=[id_col])[id_col].to_pylist()
            item['unique_ids'] = len(set(ids))
            item['ordered_ids_sha256'] = hashlib.sha256(json.dumps(ids).encode()).hexdigest()
        coverage[path.name] = item
        print('Verified core', path.name, flush=True)
    # Check every gene and missing cell against the soon-to-be-retired log matrix.
    tpm = pq.ParquetFile(TCGA / 'TCGA_STAR_TPM.parquet')
    logged = pq.ParquetFile(TCGA / 'TCGA_STAR_log2TPMplus1.parquet')
    genes = pq.read_table(TCGA / 'TCGA_STAR_TPM.genes.parquet').to_pandas()
    oldgenes = pq.read_table(TCGA / 'TCGA_STAR_log2TPMplus1.genes.parquet').to_pandas()
    pd.testing.assert_frame_equal(genes, oldgenes, check_exact=True)
    gene_ids = set(genes.gene_id)
    fields = [c for c in tpm.schema_arrow.names if c not in gene_ids]
    assert tpm.read(columns=fields).equals(logged.read(columns=fields))
    max_error = 0.0
    nulls = 0
    for start in range(0, len(genes), 256):
        columns = genes.gene_id.iloc[start:start+256].tolist()
        a = tpm.read(columns=columns, use_threads=False)
        b = logged.read(columns=columns, use_threads=False)
        for col_a, col_b in zip(a.columns, b.columns):
            assert col_a.is_null().equals(col_b.is_null())
        av = a.to_pandas().to_numpy(); bv = b.to_pandas().to_numpy()
        calculated = np.log2(av + 1)
        np.testing.assert_allclose(calculated, bv, rtol=2e-15, atol=0, equal_nan=True)
        finite = np.isfinite(bv)
        if finite.any():
            max_error = max(max_error, float(np.max(np.abs(calculated[finite]-bv[finite]))))
        nulls += int(np.isnan(av).sum())
        if start % 5120 == 0:
            print('TPM/log validation', start, '/', len(genes), flush=True)
    cn_file = TCGA / 'TCGA_GeneLevel_CN.parquet'
    cn_annotation = TCGA / 'TCGA_GeneLevel_CN.genes.parquet'
    # Force the unchanged complete-autosome algorithm, then compare old saved baselines.
    baselines, baseline_report = sample_baselines(cache_dir=WORK / 'baseline_recheck', recompute=True)
    old = pd.read_csv(ROOT / 'results/MTAP_Analysis/Provenance/Data_Audit/TCGA_Sample_Baselines.csv')
    pd.testing.assert_frame_equal(baselines, old[baselines.columns], check_dtype=False, check_exact=True)
    persist('TCGA_Sample_Baselines', baselines,
            {p.name: inputs[p.name] for p in [cn_file, cn_annotation]},
            method=BASELINE_METHOD, report=baseline_report)
    cn_gene_ids = set(pq.read_table(cn_annotation, columns=['gene_id'])['gene_id'].to_pylist())
    cn_fields = [c for c in pq.ParquetFile(cn_file).schema_arrow.names if c not in cn_gene_ids]
    cn = pq.read_table(cn_file, columns=cn_fields).to_pandas()
    rna = tpm.read(columns=fields).to_pandas()
    _, selection = select_rna(cn, rna)
    selection = selection.reset_index(drop=True)
    selection['CN_AliquotID'] = selection.SampleID.map(cn.set_index('SampleID').AliquotID)
    selection['MatchStatus'] = np.where(selection.CN_AliquotID.isna(), 'no_CN_sample',
                                       np.where(selection.Matches_CN_Aliquot, 'exact_CN_aliquot', 'sample_only'))
    selection['SelectionBasis'] = 'matching CN aliquot first; lexical FileID; one representative per SampleID'
    persist('TCGA_RNA_Selection', selection, {p.name: inputs[p.name] for p in [cn_file, TCGA/'TCGA_STAR_TPM.parquet']})
    source_map = TCGA / 'TCGA_Sample_Map.parquet'
    metadata = pq.read_table(source_map).to_pandas()
    metadata['CancerType'] = metadata.ProjectID.str.replace('TCGA-', '', regex=False)
    metadata['TumorNormal'] = np.where(metadata.SampleType.isin(TUMORS), 'Tumor',
                                     np.where(metadata.SampleType.str.contains('Normal', na=False), 'Normal', 'Other'))
    # The retired sample map is recorded as lineage, not a runtime dependency.
    persist('TCGA_Sample_Metadata', metadata,
            {p.name: inputs[p.name] for p in [cn_file, TCGA/'TCGA_STAR_TPM.parquet']},
            retired_source=dict(path=source_map.relative_to(ROOT).as_posix(), **signature(source_map)))
    for frame in [cn, rna]:
        assert set(frame.SampleID) <= set(metadata.SampleID)
        assert set(frame.CaseID) <= set(metadata.CaseID)
    cancers = sorted(set(cn.loc[cn.SampleType.isin(TUMORS), 'ProjectID'].str.removeprefix('TCGA-')))
    expected = sorted(set(json.loads((ROOT/'config/tcga_cancer_types.json').read_text())['mapping'].values()))
    assert cancers == expected and len(cancers) == 33
    report = dict(status='PREPARED_NOT_AUTHORIZED_FOR_DELETION', created_at=now(),
                  initial_project=before, coverage=coverage, cancers=cancers,
                  tcga=dict(CN_genes=len(cn_gene_ids), TPM_genes=len(genes), CN_samples=len(cn),
                            RNA_rows=len(rna), RNA_samples=rna.SampleID.nunique(),
                            metadata_rows=len(metadata), baseline_samples=len(baselines),
                            baseline_exact_match=True, TPM_log_all_genes_checked=True,
                            TPM_missing_cells=nulls, log_max_absolute_error=max_error,
                            log_tolerance='rtol=2e-15, atol=0; NaN and Arrow-null masks identical',
                            patient_deduplication=False),
                  core={p.relative_to(ROOT).as_posix(): signature(p) for p in CORE},
                  maintenance_code_sha256={p.relative_to(ROOT).as_posix():hashes(p)[0] for p in [
                      ROOT/'scripts/data_access.py', ROOT/'scripts/utils/tcga_current.py',
                      ROOT/'scripts/utils/export_workflow_tcga.py', ROOT/'scripts/R/local_data.R']})
    write_json(REPORT, report)
    print('Prepared core tables; deletion remains gated.', flush=True)


def safe_target(relative):
    path = ROOT / relative
    resolved = path.resolve(strict=True)
    data = (ROOT / 'data').resolve()
    if path.is_symlink() or not resolved.is_relative_to(data) or not path.is_file():
        raise ValueError(f'Unsafe deletion target: {relative}')
    if not relative.startswith(('data/raw/', 'data/processed/')):
        raise ValueError(f'Out of cleanup scope: {relative}')
    if resolved in CORE_RESOLVED:
        raise ValueError(f'Core deletion forbidden: {relative}')
    return path


def plan():
    report = json.loads(REPORT.read_text())
    if report['status'] != 'PREPARED_NOT_AUTHORIZED_FOR_DELETION':
        raise ValueError('Invalid preparation state')
    # Retain source lists/status records outside raw. Large API responses are not
    # needed: gdc_data_manifest.csv already preserves all source file identities.
    preserved = []
    folder = ROOT / 'data/raw/tcga/gdc_current_DR46/manifests'
    for path in sorted(folder.rglob('*')):
        if path.is_file():
            dest = MANIFESTS/'historical_DR46'/path.relative_to(folder)
            payload = path.read_bytes()
            canonical = MANIFESTS/path.name
            if canonical.exists() and canonical.read_bytes() == payload:
                dest = canonical
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists() and dest.read_bytes() != payload:
                raise ValueError('Historical provenance collision')
            dest.write_bytes(payload)
            preserved.append(dict(source=path.relative_to(ROOT).as_posix(),
                                  retained=dest.relative_to(ROOT).as_posix(), SHA256=hashes(dest)[0]))
    candidates = [p for p in (ROOT/'data/raw').rglob('*') if p.is_file()]
    candidates += [p for p in (ROOT/'data/processed').rglob('*.parquet') if p not in CORE]
    expected = {}
    for name, path_field, sha_field, md5_field in [
            ('gdc_data_manifest.csv','relative_path','SHA256','publisher_md5'),
            ('depmap_26Q1_manifest.csv','path','sha256',None),
            ('data_manifest.csv','relative_path','SHA256','publisher_md5')]:
        with (MANIFESTS/name).open(encoding='utf-8-sig',newline='') as f:
            for row in csv.DictReader(f):
                if row.get(path_field) and row.get(sha_field):
                    expected[row[path_field]]=(row[sha_field],row.get(md5_field,'') if md5_field else '')
    checkpoint=WORK/'cleanup_inventory_checkpoint.jsonl'
    entries=[json.loads(line) for line in checkpoint.read_text(encoding='utf-8').splitlines()] if checkpoint.exists() else []
    completed={e['path']:e for e in entries}
    def inspect(path):
            relative = path.relative_to(ROOT).as_posix()
            safe_target(relative)
            if relative in completed:
                old=completed[relative];stat=path.stat()
                if stat.st_size!=old['bytes'] or stat.st_mtime_ns!=old['mtime_ns']:
                    raise ValueError(f'Inventory checkpoint input changed: {relative}')
                return None
            initial=path.stat()
            sha, md5 = hashes(path)
            if path.stat().st_size!=initial.st_size or path.stat().st_mtime_ns!=initial.st_mtime_ns:
                raise ValueError(f'Inventory input changed during hashing: {relative}')
            if relative in expected:
                expected_sha, expected_md5 = expected[relative]
                if sha != expected_sha or (expected_md5 and md5.lower() != expected_md5.lower()):
                    raise ValueError(f'Original raw source hash mismatch: {relative}')
            category = '/'.join(Path(relative).parts[:4])
            entry=dict(path=relative, bytes=path.stat().st_size,
                       mtime_ns=path.stat().st_mtime_ns, SHA256=sha, MD5=md5, category=category)
            return entry
    ordered=sorted(candidates)
    with checkpoint.open('a',encoding='utf-8',newline='\n') as progress, ThreadPoolExecutor(max_workers=4) as pool:
        for start in range(0,len(ordered),256):
            for entry in pool.map(inspect,ordered[start:start+256]):
                if entry is not None:
                    entries.append(entry);progress.write(json.dumps(entry,ensure_ascii=False)+'\n');progress.flush()
            if start % 2048 == 0:
                print('Cleanup inventory SHA256/MD5',len(entries),'/',len(candidates),flush=True)
    if set(completed)-{p.relative_to(ROOT).as_posix() for p in candidates}:
        raise ValueError('Inventory checkpoint contains removed targets')
    plan_path = MANIFESTS / 'data_cleanup_files.jsonl'
    with plan_path.open('w', encoding='utf-8', newline='\n') as f:
        for entry in entries:
            f.write(json.dumps(entry, ensure_ascii=False) + '\n')
    write_json(MANIFESTS/'data_cleanup_report.json',dict(
        status='PLANNED_NOT_DELETED', created_at=now(),
        exact_file_inventory=plan_path.relative_to(ROOT).as_posix(), inventory_SHA256=hashes(plan_path)[0],
        planned_files=len(entries), planned_bytes=sum(e['bytes'] for e in entries),
        historical_provenance=preserved, retained_core=report['core'],
        deletion_prerequisites=['core hashes', 'synthetic tests', 'isolated workflows',
                                'real statistics equality', 'no active project workers']))


def no_workers():
    import psutil
    for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
        try:
            if proc.pid == os.getpid():
                continue
            name = (proc.info['name'] or '').lower()
            if name not in ['python.exe', 'rscript.exe', 'rterm.exe', 'r.exe', 'curl.exe', 'gdc-client.exe']:
                continue
            args = ' '.join(proc.info['cmdline'] or []).replace('\\', '/').lower()
            cwd = proc.cwd().replace('\\', '/').lower()
            root = ROOT.as_posix().lower()
            # Codex's external GEO MCP service can inherit this cwd. It is not
            # a project worker: allow only its exact external server entry point
            # with no project path arguments or open project data handles.
            external_geo_server = any(
                str(a).replace('\\','/').lower().endswith('/geo-mcp-server/server.py')
                and not str(a).replace('\\','/').lower().startswith(root+'/')
                for a in (proc.info['cmdline'] or []))
            if external_geo_server and root not in args:
                handles=proc.open_files()
                if not any(str(h.path).replace('\\','/').lower().startswith(root+'/data/') for h in handles):
                    continue
            if root in args or cwd.startswith(root):
                raise RuntimeError(f'Active project worker: {proc.pid} {name}')
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue


def cleanup():
    report = json.loads(REPORT.read_text())
    evidence = json.loads((MANIFESTS/'analysis_ready_regression.json').read_text())
    if evidence.get('status') != 'PASS' or not all(evidence.get(k) for k in [
            'synthetic_R', 'synthetic_Python', 'isolated_geneA_screen', 'isolated_geneA_geneB',
            'all_saved_CSVs_equal', 'no_raw_or_network', 'independent_source_check', 'PDF_visual_review']):
        raise ValueError('Deletion blocked: incomplete validation evidence')
    for relative, sig in report['core'].items():
        if signature(ROOT/relative) != sig:
            raise ValueError(f'Core changed after validation: {relative}')
    for relative, sha in evidence['code_sha256'].items():
        if hashes(ROOT/relative)[0] != sha:
            raise ValueError(f'Code changed after validation: {relative}')
    no_workers()
    cleanup_report = json.loads((MANIFESTS/'data_cleanup_report.json').read_text())
    if cleanup_report['status'] != 'PLANNED_NOT_DELETED':
        raise ValueError('Cleanup already started; inspect journal before proceeding')
    plan_path = ROOT / cleanup_report['exact_file_inventory']
    if hashes(plan_path)[0] != cleanup_report['inventory_SHA256']:
        raise ValueError('Deletion inventory changed')
    entries = [json.loads(line) for line in plan_path.read_text(encoding='utf-8').splitlines()]
    by_path={e['path']:e for e in entries}
    for entry in entries:
        if not entry['path'].endswith('.parquet'):
            continue
        p=ROOT/entry['path']
        side=p.with_suffix('.provenance.json')
        if not side.exists():side=Path(str(p)+'.provenance.json')
        if side.exists():
            meta=json.loads(side.read_text(encoding='utf-8'))
            expected_sha=meta.get('SHA256') or meta.get('output_sha256')
            if expected_sha and expected_sha!=entry['SHA256']:
                raise ValueError(f'Retired processed source SHA256 mismatch: {p.name}')
            if meta.get('raw_source') in by_path and meta.get('raw_sha256')!=by_path[meta['raw_source']]['SHA256']:
                raise ValueError(f'Retired processed/raw mapping mismatch: {p.name}')
    # Validate *all* paths and content before the first unlink. The deletion loop
    # checks stat again and journals each success; any failure stops immediately.
    def revalidate(entry):
        p = safe_target(entry['path'])
        stat = p.stat()
        if stat.st_size != entry['bytes'] or stat.st_mtime_ns != entry['mtime_ns'] or hashes(p)[0] != entry['SHA256']:
            raise ValueError(f'Deletion input changed: {entry["path"]}')
    with ThreadPoolExecutor(max_workers=4) as pool:
        for start in range(0,len(entries),256):
            list(pool.map(revalidate,entries[start:start+256]))
            if start % 4096 == 0:
                print('Pre-delete revalidation',min(start+256,len(entries)),'/',len(entries),flush=True)
    no_workers()
    cleanup_report['status'] = 'DELETING'
    write_json(MANIFESTS/'data_cleanup_report.json',cleanup_report)
    deleted = Counter(); sizes = Counter()
    journal = MANIFESTS/'data_cleanup_deleted.jsonl'
    try:
        with journal.open('x', encoding='utf-8', newline='\n') as log:
            for i, entry in enumerate(entries):
                p = safe_target(entry['path'])
                stat = p.stat()
                if stat.st_size != entry['bytes'] or stat.st_mtime_ns != entry['mtime_ns']:
                    raise ValueError(f'Deletion input changed: {entry["path"]}')
                p.unlink()
                log.write(json.dumps(dict(path=entry['path'], SHA256=entry['SHA256']))+'\n')
                log.flush()
                deleted[entry['category']] += 1
                sizes[entry['category']] += entry['bytes']
                if i % 4000 == 0:
                    print('Deleted', i, '/', len(entries), flush=True)
    except Exception:
        cleanup_report.update(status='STOPPED_ON_ERROR', deleted_files=sum(deleted.values()),
                              deleted_bytes=sum(sizes.values()), deleted_by_category=dict(deleted))
        write_json(MANIFESTS/'data_cleanup_report.json',cleanup_report)
        raise
    # Remove only explicitly enumerated empty parents, deepest first. rmdir
    # cannot remove a nonempty directory and never traverses into retained data.
    parents=set()
    for entry in entries:
        for parent in (ROOT/entry['path']).parents:
            if parent in [ROOT/'data/raw',ROOT/'data/processed']:
                parents.add(parent);break
            if parent==ROOT/'data':break
            parents.add(parent)
    empty_removed=0
    for directory in sorted(parents,key=lambda p:len(p.parts),reverse=True):
        if not directory.is_dir() or directory.is_symlink():continue
        if not directory.resolve().is_relative_to((ROOT/'data').resolve()):
            raise ValueError(f'Empty directory outside data: {directory}')
        try:directory.rmdir();empty_removed+=1
        except OSError:pass
    cleanup_report.update(status='PASS', finished_at=now(), deleted_files=sum(deleted.values()),
                          deleted_bytes=sum(sizes.values()), deleted_by_category=dict(deleted),
                          deleted_bytes_by_category=dict(sizes), empty_directories_removed=empty_removed,
                          validation_evidence='data/manifests/analysis_ready_regression.json',
                          historical_raw_integrity_note='Prior raw QC reports are historical records; deleted files are no longer present.')
    write_json(MANIFESTS/'data_cleanup_report.json',cleanup_report)
    report['status'] = 'PASS'
    report['validated_for_cleanup_at'] = now()
    write_json(REPORT,report)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['prepare', 'plan', 'cleanup'])
    {'prepare': prepare, 'plan': plan, 'cleanup': cleanup}[p.parse_args().action]()

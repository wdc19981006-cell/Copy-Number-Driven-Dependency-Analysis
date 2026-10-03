"""Processed-only DR46 baseline and exact sample-level RNA selection.

The baseline is an analysis estimate, not a published ploidy measurement.
No raw files, reference matrices, network calls or whole CN matrix loads.
"""
import hashlib
import json
from pathlib import Path

from scripts.data_access import ROOT, required
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

FOLDER = ROOT / 'data/processed/tcga/gdc_DR46'
BASELINE_METHOD = 'autosomal_gene_CN_integer_mode_v1; ties choose smallest mode'
STATES = ['Deep Deletion', 'Shallow Deletion', 'Diploid', 'Gain', 'Amplification']


def source_identity(path):
    path = required(path)
    side = required(path.with_suffix('.provenance.json'))
    provenance = json.loads(side.read_text(encoding='utf-8'))
    if provenance.get('release') != '46.0' or not provenance.get('parsed_all_inputs'):
        raise ValueError(f'Unverified DR46 processed source: {path.name}')
    return dict(path=path.name, bytes=path.stat().st_size,
                mtime_ns=path.stat().st_mtime_ns, provenance=provenance)


def classify_cn(cn, baseline):
    cn, baseline = np.asarray(cn, dtype=float), np.asarray(baseline, dtype=float)
    valid = np.isfinite(cn) & (cn >= 0) & np.isfinite(baseline) & (baseline > 0)
    state = np.full(cn.shape, np.nan)
    state[valid & (cn == 0)] = -2
    state[valid & (cn > 0) & (cn < baseline)] = -1
    state[valid & (cn == baseline)] = 0
    state[valid & (cn > baseline) & (cn < 2 * baseline)] = 1
    state[valid & (cn >= 2 * baseline)] = 2
    return state


def select_rna(cn, rna):
    if cn.SampleID.duplicated().any() or rna.FileID.duplicated().any():
        raise ValueError('CN sample or RNA file identities are not unique')
    rna = rna.copy()
    rna['Matches_CN_Aliquot'] = rna.AliquotID.eq(cn.set_index('SampleID').AliquotID.reindex(rna.SampleID).to_numpy())
    rna = rna.sort_values(['SampleID', 'Matches_CN_Aliquot', 'FileID'], ascending=[True, False, True])
    rna['Selected'] = ~rna.SampleID.duplicated()
    return rna.loc[rna.Selected].copy(), rna


def sample_baselines(folder=FOLDER, cache_dir=None, chunk_genes=512):
    """Count integer CNs in bounded gene-column blocks; cache only sample summaries."""
    path = required(folder / 'TCGA_GeneLevel_CN.parquet')
    annotation_path = required(folder / 'TCGA_GeneLevel_CN.genes.parquet')
    identity = dict(source=source_identity(path), method=BASELINE_METHOD,
                    annotation_sha256=hashlib.sha256(annotation_path.read_bytes()).hexdigest())
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    cache_dir = Path(cache_dir) if cache_dir else ROOT / '.runtime/tcga_baselines'
    cache = cache_dir / f'{key}.parquet'
    stamp = cache.with_suffix('.json')
    if cache.exists() and stamp.exists():
        saved = json.loads(stamp.read_text())
        if saved.get('sha256') == hashlib.sha256(cache.read_bytes()).hexdigest():
            return pq.read_table(cache).to_pandas(), saved['report']
    annotation = pq.read_table(annotation_path, columns=['gene_id', 'chromosome']).to_pandas()
    autosomes = {str(i) for i in range(1, 23)}
    columns = annotation.loc[annotation.chromosome.str.removeprefix('chr').isin(autosomes), 'gene_id'].tolist()
    if not columns or len(columns) != len(set(columns)):
        raise ValueError('Missing or duplicate autosomal gene annotation')
    file = pq.ParquetFile(path)
    ids = file.read(columns=['SampleID']).to_pandas()
    if ids.SampleID.duplicated().any():
        raise ValueError('CN SampleID must be unique')
    histogram = np.zeros((len(ids), 16), dtype=np.int32)
    for start in range(0, len(columns), chunk_genes):
        values = file.read(columns=columns[start:start + chunk_genes], use_threads=False).to_pandas().to_numpy()
        finite = np.isfinite(values)
        observed = values[finite]
        if np.any(observed < 0) or np.any(observed != np.rint(observed)):
            raise ValueError('Autosomal CN must contain nonnegative integers; no rounding is permitted')
        if observed.size:
            width = int(observed.max()) + 1
            if width > 10001:
                raise ValueError('Extreme CN > 10000 requires source review')
            if width > histogram.shape[1]:
                histogram = np.pad(histogram, ((0, 0), (0, width - histogram.shape[1])))
            rows = np.nonzero(finite)[0]
            np.add.at(histogram, (rows, observed.astype(np.int64)), 1)
        if start % (chunk_genes * 20) == 0:
            print(f'Baseline: {min(start + chunk_genes, len(columns))}/{len(columns)} autosomal genes', flush=True)
    total = histogram.sum(axis=1)
    mode = histogram.argmax(axis=1)
    support = histogram.max(axis=1)
    if np.any(total == 0) or np.any(mode <= 0):
        raise ValueError('Missing or nonpositive sample baseline; cannot define CNA states')
    ids['BaselineCN'] = mode
    ids['AutosomalFiniteGenes'] = total
    ids['ModalGeneCount'] = support
    ids['BaselineTies'] = (histogram == support[:, None]).sum(axis=1)
    ids['BaselineMethod'] = BASELINE_METHOD
    report = dict(method=BASELINE_METHOD, source_identity=identity, samples=len(ids),
                  autosomal_genes=len(columns), distribution={str(k): int(v) for k, v in ids.BaselineCN.value_counts().sort_index().items()},
                  tied_samples=int((ids.BaselineTies > 1).sum()),
                  extreme_samples=ids.loc[(ids.BaselineCN < 1) | (ids.BaselineCN > 8), ['SampleID', 'BaselineCN']].to_dict('records'),
                  interpretation='Estimated sample-specific baseline; not measured ploidy or official GISTIC',
                  published_ploidy='No available ploidy in wide CN/clinical/biospecimen schema; processed segment ploidy is null by its preparer')
    cache_dir.mkdir(parents=True, exist_ok=True)
    ids.to_parquet(cache, index=False)
    stamp.write_text(json.dumps(dict(sha256=hashlib.sha256(cache.read_bytes()).hexdigest(), report=report), indent=2), encoding='utf-8')
    return ids, report

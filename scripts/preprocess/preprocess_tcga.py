"""Stream gene-by-sample gzip matrices into one Parquet row group per gene.

Each gene retains a sample vector in original order. Gene statistics let DuckDB
prune unrelated row groups. No complete matrix is loaded into RAM.
"""
from pathlib import Path
import csv
import gzip
import json
import logging
import sys
import re
import os
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT, setup, verified_raw, write_json, now
from scripts.utils.preparation import output_ready, finish, parquet_qc
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import duckdb

VERSION = 'tcga-gene-vectors-v1'
MATRICES = {
 'Gistic2_CopyNumber_Gistic2_all_data_by_genes.tsv.gz': 'CN',
 'Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.tsv.gz': 'GISTIC',
 'tcga_RSEM_Hugo_norm_count.tsv.gz': 'Expression',
}


def preprocess_matrix(source, row, kind, folder):
    target = folder / f'{kind}.parquet'
    if (output_ready(target, row['SHA256'], VERSION)
            and (folder / f'{kind}.samples.parquet').exists()
            and (folder / f'{kind}.genes.parquet').exists()):
        logging.info('Verified derived file; skip %s', target.name)
        return json.loads(target.with_suffix('.parquet.provenance.json').read_text(encoding='utf-8'))['qc']
    temp = target.with_suffix('.parquet.part')
    schema = pa.schema([('Gene', pa.string()), ('GeneSymbol', pa.string()),
                        ('Values', pa.list_(pa.float64()))], metadata={
        b'orientation': b'one gene per row group; Values aligned to <dataset>.samples.parquet',
        b'source_sha256': row['SHA256'].encode(), b'compression': b'ZSTD',
        b'dataset': kind.encode()})
    genes, missing_cells = [], 0
    minimum, maximum = None, None
    with gzip.open(source, 'rt', encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        samples = header[1:]
        if not samples or any(not re.match(r'^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}', s) for s in samples):
            raise ValueError(f'{kind}: invalid TCGA sample headers')
        if len(set(samples)) != len(samples):
            raise ValueError(f'{kind}: duplicate sample IDs')
        with pq.ParquetWriter(temp, schema, compression='zstd', compression_level=6) as writer:
            batch = []
            for number, fields in enumerate(reader, 1):
                if len(fields) != len(header):
                    raise ValueError(f'{kind}: row {number} width {len(fields)} != {len(header)}')
                gene = fields[0]
                if not gene:
                    raise ValueError(f'{kind}: empty gene at row {number}')
                vector = np.asarray([v if v not in ('', 'NA', 'NaN', 'nan', 'null') else 'nan'
                                     for v in fields[1:]], dtype=np.float64)
                if np.isinf(vector).any():
                    raise ValueError(f'{kind}: infinite value at row {number}')
                nonnull = vector[~np.isnan(vector)]
                missing_cells += len(vector) - len(nonnull)
                if kind == 'GISTIC' and not np.isin(nonnull, [-2,-1,0,1,2]).all():
                    raise ValueError(f'{kind}: unexpected categorical value')
                if len(nonnull):
                    lo, hi = float(nonnull.min()), float(nonnull.max())
                    minimum = lo if minimum is None else min(minimum, lo)
                    maximum = hi if maximum is None else max(maximum, hi)
                symbol = gene.split('|')[0]
                genes.append({'Gene': gene, 'GeneSymbol': symbol, 'RowGroup': number - 1,
                              'NonMissingSamples': len(nonnull)})
                batch.append((gene, symbol, vector))
                if len(batch) == 64:
                    write_gene_batch(writer, schema, batch, len(samples))
                    batch.clear()
                if number % 2000 == 0:
                    logging.info('%s: %s genes processed', kind, number)
            if batch:
                write_gene_batch(writer, schema, batch, len(samples))
    if not genes:
        raise ValueError(f'{kind}: matrix contains no genes')
    pq.write_table(pa.table({'SampleIndex': range(len(samples)), 'SampleID': samples}),
                   folder / f'{kind}.samples.parquet', compression='zstd')
    pq.write_table(pa.Table.from_pylist(genes), folder / f'{kind}.genes.parquet', compression='zstd')
    qc = {**parquet_qc(temp), 'samples': len(samples), 'genes': len(genes),
          'duplicate_gene_labels': len(genes) - len({g['Gene'] for g in genes}),
          'missing_cells': missing_cells, 'min': minimum, 'max': maximum,
          'sample_id_format': 'passed', 'parsed_all_rows': True}
    finish(temp, target, source, row['SHA256'], VERSION, {'qc': qc})
    logging.info('%s complete: %s genes x %s samples', kind, len(genes), len(samples))
    return qc


def write_gene_batch(writer, schema, batch, sample_count):
    values = pa.array(np.concatenate([b[2] for b in batch]), from_pandas=True)
    vectors = pa.ListArray.from_arrays(pa.array(np.arange(len(batch) + 1, dtype=np.int32) * sample_count), values)
    table = pa.Table.from_arrays([pa.array([b[0] for b in batch]),
                                 pa.array([b[1] for b in batch]), vectors], schema=schema)
    writer.write_table(table, row_group_size=1)


SAMPLE_TYPES = {'01':'Primary Solid Tumor', '02':'Recurrent Solid Tumor',
                '03':'Primary Blood Derived Cancer - Peripheral Blood',
                '04':'Recurrent Blood Derived Cancer - Bone Marrow',
                '05':'Additional - New Primary', '06':'Metastatic',
                '07':'Additional Metastatic', '08':'Human Tumor Original Cells',
                '09':'Primary Blood Derived Cancer - Bone Marrow', '10':'Blood Derived Normal',
                '11':'Solid Tissue Normal', '12':'Buccal Cell Normal', '13':'EBV Immortalized Normal',
                '14':'Bone Marrow Normal', '20':'Control Analyte', '40':'Recurrent Blood Derived Cancer - Peripheral Blood'}


def preprocess_metadata(source, row, folder):
    frame = pd.read_csv(source, sep='\t', dtype=str)
    sample_col = next((c for c in ['sample', 'Sample', 'sampleID', 'sample_id', 'SampleID'] if c in frame), frame.columns[0])
    frame = frame.rename(columns={sample_col: 'SampleID'})
    if frame['SampleID'].duplicated().any():
        raise ValueError('Phenotype has duplicate sample IDs')
    all_samples = set(frame['SampleID'])
    for kind in MATRICES.values():
        sample_path = folder / f'{kind}.samples.parquet'
        if sample_path.exists():
            all_samples.update(pq.read_table(sample_path, columns=['SampleID'])['SampleID'].to_pylist())
    frame = frame.set_index('SampleID').reindex(sorted(all_samples)).reset_index()
    frame['PatientID'] = frame['SampleID'].str[:12]
    frame['SampleTypeCode'] = frame['SampleID'].str[13:15]
    frame['SampleType'] = frame['SampleTypeCode'].map(SAMPLE_TYPES)
    frame['TumorNormal'] = frame['SampleTypeCode'].map(
        lambda c: 'Tumor' if c in {'01','02','03','04','05','06','07','08','09','40'}
        else ('Normal' if c in {'10','11','12','13','14'} else 'Other'))
    cancer_col = next((c for c in ['_primary_disease', 'cancer type abbreviation', 'cancer_type', 'acronym', '_study'] if c in frame), None)
    # Do not guess cancer from tissue-site codes or use disease names as abbreviations.
    if 'cancer type abbreviation' in frame:
        frame['CancerType'] = frame['cancer type abbreviation']
    elif 'acronym' in frame:
        frame['CancerType'] = frame['acronym']
    elif '_study' in frame:
        frame['CancerType'] = frame['_study'].str.replace('TCGA-', '', regex=False)
    elif 'cancer_type' in frame:
        frame['CancerType'] = frame['cancer_type']
    elif '_primary_disease' in frame:
        code_map = json.loads((ROOT / 'config/tcga_cancer_types.json').read_text(encoding='utf-8'))['mapping']
        unknown = set(frame['_primary_disease'].dropna()) - set(code_map)
        if unknown:
            raise ValueError(f'Unmapped Xena primary-disease labels: {sorted(unknown)}')
        frame['CancerType'] = frame['_primary_disease'].map(code_map)
    else:
        raise ValueError(f'Cannot determine cancer abbreviation from phenotype columns {list(frame)}')
    # A dataset may contain a sample absent from phenotype while another sample
    # from the same patient is annotated. Only fill a cancer when that patient's
    # observed samples all agree, and retain the provenance flag.
    frame['CancerTypeSource'] = frame['CancerType'].map(lambda x: 'phenotype' if pd.notna(x) else 'unmapped')
    patient_groups = frame.dropna(subset=['CancerType']).groupby('PatientID')['CancerType']
    patient_types = {p: g.iloc[0] for p, g in patient_groups if g.nunique() == 1}
    absent = frame['CancerType'].isna()
    frame.loc[absent, 'CancerType'] = frame.loc[absent, 'PatientID'].map(patient_types)
    frame.loc[absent & frame['CancerType'].notna(), 'CancerTypeSource'] = 'same_patient_unambiguous'
    path = folder / 'Metadata.parquet'
    temp = path.with_suffix('.parquet.part')
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), temp, compression='zstd')
    finish(temp, path, source, row['SHA256'], VERSION)
    overlap = {}
    for kind in MATRICES.values():
        sample_path = folder / f'{kind}.samples.parquet'
        if sample_path.exists():
            ids = set(pq.read_table(sample_path)['SampleID'].to_pylist())
            source_ids = set(pd.read_csv(source, sep='\t', usecols=[sample_col], dtype=str)[sample_col])
            overlap[kind] = {'matrix_samples': len(ids), 'phenotype_overlap': len(ids & source_ids),
                             'missing_cancer_type': int(frame.loc[frame.SampleID.isin(ids), 'CancerType'].isna().sum())}
    return {'rows': len(frame), 'cancer_types': sorted(frame['CancerType'].dropna().unique().tolist()),
            'sample_overlap': overlap, 'tumor_normal': frame['TumorNormal'].value_counts().to_dict()}


def create_duckdb(folder):
    path = folder / 'tcga.duckdb'
    with duckdb.connect(str(path)) as con:
        con.execute("SET memory_limit='2GB'")
        for kind in MATRICES.values():
            parquet = folder / f'{kind}.parquet'
            if parquet.exists():
                escaped = parquet.as_posix().replace("'", "''")
                con.execute(f'CREATE OR REPLACE VIEW "{kind}" AS SELECT * FROM read_parquet(\'{escaped}\')')
        metadata = folder / 'Metadata.parquet'
        if metadata.exists():
            escaped = metadata.as_posix().replace("'", "''")
            con.execute(f"CREATE OR REPLACE VIEW Metadata AS SELECT * FROM read_parquet('{escaped}')")
    # The DB contains views only, not another copy of the matrix.


def main():
    setup('preprocess_tcga')
    folder = ROOT / 'data/processed/tcga/pancanatlas_reference'
    folder.mkdir(parents=True, exist_ok=True)
    qc, failures = {}, {}
    for filename, kind in MATRICES.items():
        try:
            raw = verified_raw('TCGA', 'PANCAN', filename)
            if not raw:
                failures[kind] = 'Verified raw file unavailable'
                continue
            qc[kind] = preprocess_matrix(*raw, kind, folder)
        except Exception as e:
            logging.exception('%s preprocessing failed', kind)
            failures[kind] = str(e)
    try:
        raw = verified_raw('TCGA', 'PANCAN', 'TCGA_phenotype_denseDataOnlyDownload.tsv.gz')
        if raw:
            qc['metadata'] = preprocess_metadata(*raw, folder)
        else:
            failures['metadata'] = 'Verified raw file unavailable'
    except Exception as e:
        logging.exception('Metadata preprocessing failed')
        failures['metadata'] = str(e)
    create_duckdb(folder)
    write_json(ROOT / 'data/manifests/tcga_qc.json', {'created_at': now(), 'qc': qc, 'failures': failures})
    print(json.dumps({'qc': qc, 'failures': failures}, indent=2))
    return bool(failures)


if __name__ == '__main__':
    raise SystemExit(main())

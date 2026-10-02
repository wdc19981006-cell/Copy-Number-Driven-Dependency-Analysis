"""Integrity and query smoke checks, not biological dependency analysis."""
from pathlib import Path
import csv
import gzip
import json
import sys
import time
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT, setup, load_manifest, hashes, now, write_json, current_release
from scripts.data_access import get_tcga_cn, get_tcga_gistic, get_tcga_expression, get_tcga_metadata, get_depmap_dependency, DataUnavailableError
import pyarrow.parquet as pq
import numpy as np


def raw_gene(filename, gene):
    path = ROOT / 'data/raw/tcga/pancanatlas_reference' / filename
    with gzip.open(path, 'rt', encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        samples = next(reader)[1:]
        for fields in reader:
            if fields[0] == gene:
                values = [float(v) if v not in ('', 'NA', 'NaN', 'nan', 'null') else np.nan for v in fields[1:]]
                return samples, np.array(values)
    raise KeyError(gene)


def main():
    setup('validate_database')
    checks, errors = [], []
    for row in load_manifest():
        if row['status'] != 'verified':
            continue
        path = ROOT / row['relative_path']
        ok = path.exists() and path.stat().st_size == int(row['file_size']) and hashes(path)[0] == row['SHA256']
        checks.append({'check': 'raw_integrity', 'file': row['filename'], 'passed': ok})
        if not ok:
            errors.append(f'Raw integrity failed: {path}')
    for provenance in (ROOT / 'data/processed').rglob('*.provenance.json'):
        info = json.loads(provenance.read_text(encoding='utf-8'))
        path = provenance.with_name(provenance.name.removesuffix('.provenance.json'))
        ok = path.exists() and hashes(path)[0] == info['output_sha256']
        checks.append({'check': 'processed_integrity', 'file': path.name, 'passed': ok})
        if not ok:
            errors.append(f'Processed integrity failed: {path}')
    # Three different genes for each dataset, selected to exercise exact matching
    # near different parts of the original alphabetic gene order.
    for dataset, filename, query in [
        ('CN', 'Gistic2_CopyNumber_Gistic2_all_data_by_genes.tsv.gz', get_tcga_cn),
        ('GISTIC', 'Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.tsv.gz', get_tcga_gistic),
        ('Expression', 'tcga_RSEM_Hugo_norm_count.tsv.gz', get_tcga_expression),
    ]:
        gene_index = ROOT / f'data/processed/tcga/pancanatlas_reference/{dataset}.genes.parquet'
        if not gene_index.exists():
            errors.append(f'Missing prepared dataset: {dataset}')
            continue
        genes = pq.read_table(gene_index, columns=['Gene'])['Gene'].to_pylist()
        for gene in [genes[0], genes[len(genes)//2], genes[-1]]:
            try:
                begin = time.perf_counter()
                result = query(gene, layer='reference')
                duration = time.perf_counter() - begin
                samples, values = raw_gene(filename, gene)
                np.testing.assert_array_equal(result.SampleID.to_numpy(), np.array(samples))
                np.testing.assert_allclose(result.iloc[:, 1].to_numpy(), values, rtol=0, atol=0, equal_nan=True)
                if dataset == 'GISTIC':
                    assert set(result.iloc[:, 1].dropna()).issubset({-2,-1,0,1,2})
                checks.append({'check': 'single_gene_matches_raw', 'dataset': dataset, 'gene': gene,
                               'samples': len(result), 'query_seconds': round(duration, 3), 'passed': True})
            except Exception as e:
                errors.append(f'{dataset} query {gene}: {e}')
    try:
        meta = get_tcga_metadata(layer='reference')
        assert not meta.SampleID.duplicated().any()
        assert meta.SampleID.str.match(r'^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}').all()
        assert (meta.PatientID == meta.SampleID.str[:12]).all()
        assert meta.CancerType.dropna().nunique() == 33
        for dataset in ['CN', 'GISTIC', 'Expression']:
            ids = pq.read_table(ROOT / f'data/processed/tcga/pancanatlas_reference/{dataset}.samples.parquet')['SampleID'].to_pylist()
            aligned = meta.set_index('SampleID').reindex(ids)
            missing = int(aligned.CancerType.isna().sum())
            unclassified_biological = int((aligned.CancerType.isna() & aligned.TumorNormal.isin(['Tumor', 'Normal'])).sum())
            expected_controls = int((aligned.CancerType.isna() & (aligned.SampleTypeCode == '20')).sum())
            checks.append({'check': 'matrix_metadata_coverage', 'dataset': dataset,
                           'samples': len(ids), 'missing_cancer_type': missing,
                           'unclassified_control_analytes': expected_controls,
                           'unclassified_tumor_normal': unclassified_biological,
                           'passed': missing == expected_controls})
            if missing != expected_controls:
                errors.append(f'{dataset}: {missing - expected_controls} biological/unknown samples have no cancer classification')
    except Exception as e:
        errors.append(f'Metadata validation: {e}')
    # If DepMap was skipped, queries must fail clearly, never substitute TCGA or old-release data.
    config = json.loads((ROOT / 'config/current_release.json').read_text(encoding='utf-8'))
    if config.get('installation_status') == 'skipped_by_user':
        try:
            get_depmap_dependency('TP53')
            errors.append('Skipped DepMap incorrectly returned a dataset')
        except DataUnavailableError:
            checks.append({'check': 'skipped_depmap_fails_explicitly', 'passed': True})
    result = {'created_at': now(), 'checks': checks, 'errors': errors,
              'passed': not errors, 'depmap_model_overlap': 'Not measured: user skipped DepMap downloads'}
    write_json(ROOT / 'data/manifests/validation_report.json', result)
    print(json.dumps(result, indent=2))
    return bool(errors)


if __name__ == '__main__':
    raise SystemExit(main())

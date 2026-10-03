"""Validate exact processed-source cohorts and PDF integrity without rerunning modules."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.data_access import (depmap_gene, get_depmap_dependency, get_tcga_cn,
                                 get_tcga_expression, gdc_gene)
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from collections import Counter


def validate_current(case, gene):
    """Independent selection/classification and processed-column re-reads."""
    folder = ROOT / 'data/processed/tcga/gdc_DR46'
    saved = pd.read_csv(case / 'Tables/TCGA_Current_Samples.csv')
    baseline = pd.read_csv(case / 'Tables/TCGA_Sample_Baselines.csv')
    cn = get_tcga_cn(gene, layer='current')
    rna = gdc_gene('TCGA_STAR_TPM', gene, 'TPM')
    assert not saved.SampleID.duplicated().any() and not baseline.SampleID.duplicated().any()
    assert saved.DataLayer.eq('gdc_current_DR46').all()
    direct = cn.set_index('SampleID').loc[saved.SampleID]
    assert direct.FileID.tolist() == saved.FileID.tolist()
    np.testing.assert_allclose(direct.CopyNumber, saved.CopyNumber, rtol=1e-13, equal_nan=True)
    for column in ['CaseID', 'ProjectID', 'SampleType', 'AliquotID']:
        assert direct[column].tolist() == saved[column].tolist()
    preferred = cn.set_index('SampleID').AliquotID.to_dict()
    rna['preference'] = [a == preferred.get(s) for s, a in zip(rna.SampleID, rna.AliquotID)]
    rna = rna.sort_values(['SampleID', 'preference', 'FileID'], ascending=[True, False, True]).drop_duplicates('SampleID').set_index('SampleID')
    exact = rna.reindex(saved.SampleID)
    assert exact.FileID.fillna('').tolist() == saved.RNA_FileID.fillna('').tolist()
    np.testing.assert_allclose(exact.TPM, saved.RNA_TPM, rtol=1e-13, equal_nan=True)
    np.testing.assert_allclose(np.log2(exact.TPM+1), saved.Expression, rtol=1e-13, equal_nan=True)
    # A separately prepared log2TPM matrix must agree for the chosen RNA files.
    logged = get_tcga_expression(gene, layer='current').set_index('FileID').Expression
    use = saved.RNA_FileID.notna()
    np.testing.assert_allclose(logged.loc[saved.loc[use, 'RNA_FileID']], saved.loc[use, 'Expression'], rtol=1e-13)
    states = []
    for value, base in zip(saved.CopyNumber, saved.BaselineCN):
        assert np.isfinite(base) and base > 0 and base == int(base)
        states.append(np.nan if not np.isfinite(value) else -2 if value == 0 else
                      -1 if value < base else 0 if value == base else 1 if value < 2*base else 2)
    np.testing.assert_equal(np.asarray(states), saved.CNAState.to_numpy())
    # Random sample baselines: different accumulator implementation (Counter),
    # rereading source columns without using the baseline cache or exporter.
    selected = set(baseline.sample(n=min(20, len(baseline)), random_state=4106).SampleID)
    selected.update(baseline.loc[(baseline.BaselineCN > 8) | (baseline.BaselineTies > 1), 'SampleID'])
    pf = pq.ParquetFile(folder / 'TCGA_GeneLevel_CN.parquet')
    ids = pf.read(columns=['SampleID']).column('SampleID').to_pylist()
    indices = [i for i, s in enumerate(ids) if s in selected]
    genes = pq.read_table(folder/'TCGA_GeneLevel_CN.genes.parquet').to_pandas()
    columns = genes.loc[genes.chromosome.str.replace('chr', '', regex=False).isin([str(i) for i in range(1, 23)]), 'gene_id'].tolist()
    counters = [Counter() for _ in indices]
    for start in range(0, len(columns), 384):
        values = pf.read(columns=columns[start:start+384], use_threads=False).take(indices).to_pandas().to_numpy()
        for row, counter in zip(values, counters):
            counter.update(row[np.isfinite(row)].tolist())
    indexed = baseline.set_index('SampleID')
    for i, counter in zip(indices, counters):
        mode = sorted(counter, key=lambda value: (-counter[value], value))[0]
        assert mode == indexed.loc[ids[i], 'BaselineCN']
        assert sum(counter.values()) == indexed.loc[ids[i], 'AutosomalFiniteGenes']
    return dict(current_source_exact_sample_join=True, STAR_TPM_and_log2_processed_agree=True,
                analysis_defined_CN_states=True, unique_positive_integer_baselines=True,
                independently_reread_baseline_samples=len(indices), random_sample_seed=4106,
                baseline_distribution={str(k): int(v) for k, v in baseline.BaselineCN.value_counts().sort_index().items()},
                tumor_samples=len(saved), finite_gene_CN=int(np.isfinite(saved.CopyNumber).sum()))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', default='VPS4B_VPS4A_Analysis')
    a = p.parse_args()
    case = (ROOT / 'results' / a.case).resolve()
    if not case.is_relative_to((ROOT / 'results').resolve()):
        raise ValueError('Case outside results')
    provenance = case / 'Provenance'
    info = json.loads((provenance / 'Run_Metadata.json').read_text(encoding='utf-8'))
    gene_a, gene_b = info['geneA'], info['geneB']
    checks = {}
    checks.update(validate_current(case, gene_a))
    assert info['TCGA_main_layer'] == 'gdc_current_DR46'
    versions = json.loads((provenance/'Input_Versions.json').read_text())['tcga']
    assert all('data/processed/tcga/gdc_DR46/' in entry['path'] for entry in versions)
    assert not any(term in json.dumps(versions).lower() for term in ['pancanatlas', 'xena', 'toil'])
    checks['all_main_TCGA_actual_inputs_current_DR46'] = True
    cn = depmap_gene('OmicsCNGeneWGS', gene_a, 'CN_relative')
    checks['DepMap_processed_readable'] = len(cn) > 0
    if gene_b:
        direct = cn.merge(get_depmap_dependency(gene_b).rename(columns={'GeneEffect': 'Chronos'}),
                          on='ModelID', validate='one_to_one')
        direct['CN_log'] = np.log2(direct.CN_relative + 1)
        direct = direct.loc[np.isfinite(direct.CN_log) & np.isfinite(direct.Chronos)]
        pair = pd.read_csv(case / 'Tables' / f'{gene_a}_{gene_b}_CellLines.csv')
        columns = ['ModelID', 'CN_relative', 'CN_log', 'Chronos']
        pd.testing.assert_frame_equal(direct[columns].sort_values('ModelID').reset_index(drop=True),
                                      pair[columns].sort_values('ModelID').reset_index(drop=True),
                                      check_dtype=False, rtol=1e-12, atol=1e-12)
        checks['targeted_processed_source_exact_model_join'] = True
    index = set(pd.read_csv(provenance / 'File_Index.csv').Path)
    actual = {p.relative_to(case).as_posix() for p in case.rglob('*') if p.is_file()}
    assert index == actual
    checks['file_index_exact'] = True
    summary = (case / '00_Analysis_Summary.txt').read_text(encoding='utf-8')
    assert all(name in summary for name in index)
    checks['summary_lists_every_file'] = True
    # Available local PDF renderer checks that every generated PDF opens cleanly.
    from pypdf import PdfReader
    pdfs = sorted((case / 'Main_Results').rglob('*.pdf'))
    for path in pdfs:
        doc = PdfReader(path)
        assert len(doc.pages) == 1
        assert len(doc.pages[0].extract_text()) > 20
    checks['all_main_PDFs_valid_single_page'] = True
    # Validation.json is the existing indexed report; updating it does not add a
    # file that would require analysis or index regeneration.
    report = json.loads((provenance / 'Validation.json').read_text(encoding='utf-8'))
    report['processed_source_checks'] = checks
    report['source_validation_status'] = 'PASS'
    report['main_pdf_count'] = len(pdfs)
    (provenance / 'Validation.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'main_PDFs': len(pdfs), 'checks': checks}, indent=2))


if __name__ == '__main__':
    main()

"""Validate exact processed-source cohorts and PDF integrity without rerunning modules."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.utils.export_workflow_tcga import reference_samples
from scripts.data_access import (depmap_gene, get_depmap_dependency, get_tcga_cn,
                                 get_tcga_expression)
import numpy as np
import pandas as pd


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
    ref = reference_samples(gene_a)
    saved = pd.read_csv(case / 'Tables/TCGA_Reference_Samples.csv')
    pd.testing.assert_frame_equal(ref.sort_values('SampleID').reset_index(drop=True),
                                  saved.sort_values('SampleID').reset_index(drop=True),
                                  check_dtype=False, rtol=1e-12, atol=1e-12)
    checks['reference_source_exact_sample_join'] = True
    current_cn = get_tcga_cn(gene_a, layer='current')
    current_rna = get_tcga_expression(gene_a, layer='current')
    assert len(current_cn) and len(current_rna)
    checks['DR46_processed_CN_and_expression_readable'] = True
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

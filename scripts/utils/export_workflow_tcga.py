"""Local-only single-gene DR46 extraction for default TCGA figures 01-03."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.data_access import get_tcga_cn, gdc_gene
from scripts.utils.tcga_current import (FOLDER, STATES, classify_cn, sample_baselines,
                                       select_rna, source_identity)
import numpy as np


def current_samples(gene):
    source_identity(FOLDER / 'TCGA_STAR_TPM.parquet')
    cn = get_tcga_cn(gene, layer='current')
    rna = gdc_gene('TCGA_STAR_TPM', gene, 'TPM')
    selected, audit = select_rna(cn, rna)
    baseline, report = sample_baselines()
    fields = ['SampleID', 'CaseID', 'ProjectID', 'AliquotID', 'FileID', 'TPM']
    selected = selected[fields].rename(columns={c: 'RNA_' + c for c in fields if c != 'SampleID'})
    frame = cn.merge(baseline, on='SampleID', validate='one_to_one')
    frame = frame.merge(selected, on='SampleID', how='left', validate='one_to_one')
    matched = frame.RNA_FileID.notna()
    if not ((frame.loc[matched, 'CaseID'] == frame.loc[matched, 'RNA_CaseID']).all()
            and (frame.loc[matched, 'ProjectID'] == frame.loc[matched, 'RNA_ProjectID']).all()):
        raise ValueError('Exact SampleID match has inconsistent case/project metadata')
    if (frame.RNA_TPM.dropna() < 0).any():
        raise ValueError('Negative STAR TPM')
    tumors = ['Primary Tumor', 'Recurrent Tumor', 'Metastatic', 'Additional Metastatic',
              'Additional - New Primary', 'Primary Blood Derived Cancer - Peripheral Blood',
              'Primary Blood Derived Cancer - Bone Marrow']
    frame = frame.loc[frame.SampleType.isin(tumors)].copy()
    frame['TumorNormal'] = 'Tumor'
    frame['CancerType'] = frame.ProjectID.str.replace('TCGA-', '', regex=False)
    frame['Expression'] = np.log2(frame.RNA_TPM + 1)
    frame['CNAState'] = classify_cn(frame.CopyNumber, frame.BaselineCN)
    frame['CNA'] = frame.CNAState.map(dict(zip(range(-2, 3), STATES)))
    frame['DataLayer'] = 'gdc_current_DR46'
    frame['CNUnit'] = 'absolute gene-level copy number'
    frame['ExpressionUnit'] = 'STAR log2(TPM + 1)'
    audit['DataLayer'] = 'gdc_current_DR46'
    return frame, baseline, report, audit


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gene', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    target = Path(a.output).resolve()
    if not target.is_relative_to((ROOT / 'results').resolve()):
        raise ValueError('Output must be inside project results')
    frame, baseline, report, audit = current_samples(a.gene)
    target.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(target, index=False)
    baseline.to_csv(target.parent / 'TCGA_Sample_Baselines.csv', index=False)
    audit.to_csv(target.parent / 'TCGA_RNA_Representative_Selection.csv', index=False)
    (ROOT / 'results' / target.relative_to((ROOT / 'results').resolve()).parts[0] / 'Provenance/TCGA_Baseline_Method.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({'samples': len(frame), 'finite_CN': int(np.isfinite(frame.CopyNumber).sum()),
                      'matched_RNA': int(np.isfinite(frame.Expression).sum()),
                      'cancers': int(frame.CancerType.nunique()), 'layer': 'gdc_current_DR46'}))


if __name__ == '__main__':
    main()

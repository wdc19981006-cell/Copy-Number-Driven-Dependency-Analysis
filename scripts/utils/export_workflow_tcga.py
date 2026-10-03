"""Local-only, single-gene reference extraction; never inspect raw or use a network."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from scripts.data_access import (get_tcga_cn, get_tcga_expression,
                                 get_tcga_gistic, get_tcga_metadata)


def reference_samples(gene):
    meta = get_tcga_metadata(layer='reference')
    fields = ['SampleID', 'CancerType', 'SampleType', 'TumorNormal']
    frame = get_tcga_cn(gene, layer='reference').merge(meta[fields], on='SampleID',
                                                     validate='one_to_one')
    frame = frame.merge(get_tcga_gistic(gene, layer='reference'), on='SampleID',
                        how='left', validate='one_to_one')
    frame = frame.merge(get_tcga_expression(gene, layer='reference'), on='SampleID',
                        how='left', validate='one_to_one')
    frame = frame.loc[(frame.TumorNormal == 'Tumor') & frame.CancerType.notna()].copy()
    finite = frame.loc[np.isfinite(frame.GISTIC), 'GISTIC']
    if not set(finite).issubset({-2, -1, 0, 1, 2}):
        raise ValueError('Reference GISTIC must have exactly the five source codes')
    frame['DataLayer'] = 'PanCanAtlas_Xena_reference'
    frame['CNUnit'] = 'source continuous GISTIC2 scale'
    frame['ExpressionUnit'] = 'log2(norm_count+1); not TPM'
    return frame


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gene', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    target = Path(a.output).resolve()
    if not target.is_relative_to((ROOT / 'results').resolve()):
        raise ValueError('Output must be inside project results')
    frame = reference_samples(a.gene)
    target.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(target, index=False)
    print(json.dumps({'samples': len(frame), 'cancers': int(frame.CancerType.nunique()),
                      'layer': 'PanCanAtlas_Xena_reference'}))


if __name__ == '__main__':
    main()

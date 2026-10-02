"""Uniform gene query API. Paths and release selection are centralized here."""
from pathlib import Path
import json
import re
import sys
import os
from functools import lru_cache
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.utils.common import ROOT, current_release
import pandas as pd
import pyarrow.parquet as pq
import duckdb


class DataUnavailableError(FileNotFoundError):
    pass


def required(path):
    if not path.exists():
        raise DataUnavailableError(f'Prepared dataset unavailable: {path}. See data/manifests/data_manifest.csv and run the matching preparation script.')
    return path


def depmap_path(dataset):
    return ROOT / 'data/processed/depmap' / current_release() / f'{dataset}.parquet'


def resolve_column(path, gene):
    names = pq.ParquetFile(path).schema_arrow.names
    if gene in names and gene != 'ModelID':
        return gene
    matches = [c for c in names if c != 'ModelID' and
               (re.sub(r'\s+\([^)]*\)$', '', c) == str(gene) or
                (re.search(r'\((\d+)\)$', c) and re.search(r'\((\d+)\)$', c)[1] == str(gene)))]
    if not matches:
        raise KeyError(f'Gene {gene!r} is absent from {path.name}')
    if len(matches) != 1:
        raise ValueError(f'Ambiguous gene {gene!r}: {matches}. Supply exact original label or Entrez ID.')
    return matches[0]


def depmap_gene(dataset, gene, value_name):
    path = required(depmap_path(dataset))
    column = resolve_column(path, gene)
    frame = pq.read_table(path, columns=['ModelID', column]).to_pandas()
    frame = frame.rename(columns={column: value_name})
    frame.attrs.update(gene=gene, original_gene_column=column, release=current_release(), dataset=dataset)
    return frame


def get_depmap_cn(gene):
    primary = depmap_gene('OmicsCNGeneWGS', gene, 'CopyNumber')
    primary['CNSource'] = 'OmicsCNGeneWGS'
    supplemental = depmap_path('OmicsCNGene')
    if supplemental.exists():
        try:
            fallback = depmap_gene('OmicsCNGene', gene, 'CopyNumber')
        except KeyError:
            return primary
        # Fallback is only for models absent from WGS, not a WGS value missing for this gene.
        fallback = fallback.loc[~fallback.ModelID.isin(primary.ModelID)].copy()
        fallback['CNSource'] = 'OmicsCNGene'
        primary = pd.concat([primary, fallback], ignore_index=True)
    return primary


def get_depmap_expression(gene):
    return depmap_gene('ExpressionProteinCoding', gene, 'Expression')


def get_depmap_dependency(gene, probability=False):
    return depmap_gene('CRISPRGeneDependency' if probability else 'CRISPRGeneEffect', gene,
                       'DependencyProbability' if probability else 'GeneEffect')


def get_depmap_model_metadata():
    return pq.read_table(required(depmap_path('Model'))).to_pandas()


def get_depmap_mutation(gene, kind='damaging'):
    if kind not in ('damaging', 'hotspot'):
        raise ValueError('kind must be damaging or hotspot')
    return depmap_gene('MutationDamaging' if kind == 'damaging' else 'MutationHotspot', gene, 'Mutation')


@lru_cache(maxsize=3)
def tcga_index(path, modified_ns):
    index = pq.read_table(path, columns=['Gene', 'GeneSymbol', 'RowGroup']).to_pydict()
    exact, symbols = {}, {}
    for label, symbol, row_group in zip(index['Gene'], index['GeneSymbol'], index['RowGroup']):
        exact.setdefault(label, []).append(row_group)
        symbols.setdefault(symbol, []).append(row_group)
    return exact, symbols


@lru_cache(maxsize=3)
def tcga_parquet(path, modified_ns):
    # Cache the footer only, never the full matrix.
    return pq.ParquetFile(path, memory_map=True)


def tcga_gene(dataset, gene):
    folder = ROOT / 'data/processed/tcga/pancanatlas_reference'
    path = required(folder / f'{dataset}.parquet')
    index_path = required(folder / f'{dataset}.genes.parquet')
    exact, symbols = tcga_index(str(index_path), index_path.stat().st_mtime_ns)
    groups = exact.get(str(gene)) or symbols.get(str(gene))
    if not groups:
        raise KeyError(f'Gene {gene!r} absent from TCGA {dataset}')
    if len(groups) != 1:
        raise ValueError(f'Ambiguous gene {gene!r}: row groups {groups}; supply the exact source label.')
    file = tcga_parquet(str(path), path.stat().st_mtime_ns)
    rows = file.read_row_group(groups[0], columns=['Gene', 'Values']).to_pylist()
    if len(rows) != 1 or (rows[0]['Gene'] != str(gene) and rows[0]['Gene'].split('|')[0] != str(gene)):
        raise ValueError('Gene index and data file do not agree. Regenerate processed files.')
    samples = pq.read_table(required(folder / f'{dataset}.samples.parquet'), columns=['SampleID'])['SampleID'].to_pylist()
    label = {'CN': 'CopyNumber', 'GISTIC': 'GISTIC', 'Expression': 'Expression'}[dataset]
    frame = pd.DataFrame({'SampleID': samples, label: rows[0]['Values']})
    frame.attrs.update(gene=gene, source_gene=rows[0]['Gene'], dataset=dataset,
                       unit='log2(norm_count+1)' if dataset == 'Expression' else 'source GISTIC2 scale')
    return frame


def get_tcga_cn(gene, layer='current'):
    if layer == 'current':
        return gdc_gene('TCGA_GeneLevel_CN', gene, 'CopyNumber')
    if layer != 'reference':
        raise ValueError('layer must be current or reference')
    return tcga_gene('CN', gene)


def get_tcga_gistic(gene, layer='reference'):
    if layer != 'reference':
        raise ValueError('Thresholded GISTIC is available only in the reference layer')
    return tcga_gene('GISTIC', gene)


def get_tcga_expression(gene, layer='current'):
    if layer == 'current':
        return gdc_gene('TCGA_STAR_log2TPMplus1', gene, 'Expression')
    if layer != 'reference':
        raise ValueError('layer must be current or reference')
    return tcga_gene('Expression', gene)


def get_tcga_metadata(layer='current'):
    if layer=='current':
        return pq.read_table(required(ROOT/'data/processed/tcga/gdc_DR46/TCGA_Sample_Map.parquet')).to_pandas()
    if layer!='reference':
        raise ValueError('layer must be current or reference')
    return pq.read_table(required(ROOT / 'data/processed/tcga/pancanatlas_reference/Metadata.parquet')).to_pandas()


def gdc_gene(dataset, gene, value_name):
    folder = ROOT/'data/processed/tcga/gdc_DR46'
    path = required(folder/f'{dataset}.parquet')
    annotation = pq.read_table(required(folder/f'{dataset}.genes.parquet'),columns=['gene_id','gene_name']).to_pandas()
    matches = annotation.loc[(annotation.gene_id==str(gene)) | (annotation.gene_name==str(gene)), 'gene_id'].tolist()
    if str(gene) in matches:
        matches=[str(gene)]
    if not matches:
        # Accept a versionless Ensembl ID only if it uniquely resolves.
        matches=annotation.loc[annotation.gene_id.str.split('.').str[0]==str(gene),'gene_id'].tolist()
    if len(matches)!=1:
        raise KeyError(f'Absent or ambiguous gene {gene!r}; resolved IDs: {matches}')
    columns=['SampleID','SampleBarcode','AliquotID','AliquotBarcode','CaseID','CaseBarcode','ProjectID','SampleType','FileID','Workflow',matches[0]]
    schema=pq.ParquetFile(path).schema_arrow.names
    frame=pq.read_table(path,columns=[c for c in columns if c in schema]).to_pandas().rename(columns={matches[0]:value_name})
    frame.attrs.update(gene=gene, gene_id=matches[0], layer='gdc_DR46',release='46.0',
                       unit='log2(TPM+1)' if 'log2' in dataset else ('TPM' if 'TPM' in dataset else 'source total copy number'))
    return frame


def request_gene(gene, database='tcga', dataset='cn'):
    functions = {('tcga', 'cn'): get_tcga_cn, ('tcga', 'gistic'): get_tcga_gistic,
                 ('tcga', 'expression'): get_tcga_expression, ('depmap', 'cn'): get_depmap_cn,
                 ('depmap', 'expression'): get_depmap_expression,
                 ('depmap', 'dependency'): get_depmap_dependency, ('depmap', 'mutation'): get_depmap_mutation}
    if (database, dataset) not in functions:
        raise ValueError('Unknown database/dataset combination')
    return functions[(database, dataset)](gene)

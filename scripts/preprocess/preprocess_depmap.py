"""Convert verified selected-release CSVs; fails safely on ambiguous model mapping."""
from pathlib import Path
import csv
import json
import logging
import re
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT, current_release, setup, verified_raw, write_json, now
from scripts.utils.preparation import output_ready, finish, parquet_qc
import pyarrow as pa
import pyarrow.csv as pcsv
import pyarrow.parquet as pq
import pandas as pd

VERSION = 'depmap-model-columns-v1'
WIDE = {
 'CRISPRGeneEffect.csv': 'CRISPRGeneEffect',
 'CRISPRGeneDependency.csv': 'CRISPRGeneDependency',
 'OmicsCNGeneWGS.csv': 'OmicsCNGeneWGS',
 'OmicsCNGene.csv': 'OmicsCNGene',
 'OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv': 'ExpressionProteinCoding',
 'OmicsExpressionTPMLogp1HumanAllGenesStranded.csv': 'ExpressionAllGenesStranded',
 'OmicsSomaticMutationsMatrixDamaging.csv': 'MutationDamaging',
 'OmicsSomaticMutationsMatrixHotspot.csv': 'MutationHotspot',
}
TABLES = {
 'Model.csv': 'Model', 'ModelCondition.csv': 'ModelCondition',
 'OmicsProfiles.csv': 'OmicsProfiles', 'CRISPRScreenMap.csv': 'CRISPRScreenMap',
 'ScreenSequenceMap.csv': 'ScreenSequenceMap',
 'OmicsSomaticMutations.csv': 'SomaticMutations',
 'OmicsGlobalSignatures.csv': 'GlobalSignatures',
 'OmicsInferredMolecularSubtypes.csv': 'MolecularSubtypes',
 'OmicsFusionFilteredSupplementary.csv': 'FusionSupplementary',
}


def stream_table(source, row, target):
    if output_ready(target, row['SHA256'], VERSION):
        return parquet_qc(target)
    with source.open(encoding='utf-8-sig', newline='') as f:
        columns = next(csv.reader(f))
    # Keep heterogeneous metadata and variant annotations as strings without lossy inference.
    types = {c: pa.string() for c in columns}
    reader = pcsv.open_csv(source, read_options=pcsv.ReadOptions(block_size=32 * 1024 * 1024),
                          convert_options=pcsv.ConvertOptions(column_types=types, strings_can_be_null=True))
    temp = target.with_suffix('.parquet.part')
    with pq.ParquetWriter(temp, reader.schema, compression='zstd', compression_level=6) as writer:
        for batch in reader:
            writer.write_batch(batch)
    if pq.ParquetFile(temp).metadata.num_rows == 0:
        raise ValueError('Empty CSV table')
    qc = parquet_qc(temp)
    finish(temp, target, source, row['SHA256'], VERSION, {'qc': qc})
    return qc


def mapping(folder):
    profiles = folder / 'OmicsProfiles.parquet'
    conditions = folder / 'ModelCondition.parquet'
    pmap, defaults, cmap = {}, set(), {}
    if profiles.exists():
        frame = pq.read_table(profiles).to_pandas()
        profile_col = next((c for c in ['ProfileID', 'OmicsProfileID'] if c in frame), None)
        if profile_col and 'ModelID' in frame:
            for profile, group in frame.groupby(profile_col):
                ids = group['ModelID'].dropna().unique()
                if len(ids) != 1:
                    raise ValueError(f'Profile {profile} has ambiguous ModelID')
                pmap[profile] = ids[0]
            default_col = next((c for c in ['IsDefaultEntryForModel', 'IsDefaultProfile', 'IsDefault'] if c in frame), None)
            if default_col:
                defaults = set(frame.loc[frame[default_col].astype(str).str.lower().isin(['true','1','yes']), profile_col])
    if conditions.exists():
        frame = pq.read_table(conditions).to_pandas()
        if {'ModelConditionID', 'ModelID'}.issubset(frame):
            for condition, group in frame.groupby('ModelConditionID'):
                ids = group['ModelID'].dropna().unique()
                if len(ids) != 1:
                    raise ValueError(f'Condition {condition} has ambiguous ModelID')
                cmap[condition] = ids[0]
    return pmap, defaults, cmap


def stream_matrix(source, row, target, maps, model_ids):
    if output_ready(target, row['SHA256'], VERSION):
        ids = pq.read_table(target, columns=['ModelID'])['ModelID'].to_pylist()
        return {**parquet_qc(target), 'unique_models': len(set(ids)),
                'model_metadata_overlap': len(set(ids) & model_ids),
                'model_metadata_unmatched': len(set(ids) - model_ids), 'ModelID_format': 'passed'}
    with source.open(encoding='utf-8-sig', newline='') as f:
        header = next(csv.reader(f))
    if len(header) != len(set(header)):
        raise ValueError('Duplicate matrix column labels are ambiguous')
    ids_col = 'ModelID' if 'ModelID' in header else header[0]
    other_ids = {c for c in header if c in ['ModelID','ModelConditionID','ProfileID','OmicsProfileID','SequencingID']}
    other_ids.add(ids_col)
    genes = [c for c in header if c not in other_ids]
    if not genes:
        raise ValueError('Matrix has no gene columns')
    types = {c: (pa.string() if c in other_ids else pa.float64()) for c in header}
    reader = pcsv.open_csv(source, read_options=pcsv.ReadOptions(block_size=32 * 1024 * 1024),
                          convert_options=pcsv.ConvertOptions(column_types=types,
                                                             null_values=['', 'NA', 'NaN', 'nan']))
    schema = pa.schema([('ModelID', pa.string())] + [(g, pa.float64()) for g in genes], metadata={
        b'source_sha256': row['SHA256'].encode(), b'orientation': b'models x genes',
        b'numeric_precision': b'float64, source scale preserved'})
    temp = target.with_suffix('.parquet.part')
    pmap, defaults, cmap = maps
    seen, skipped_defaults, invalid = set(), 0, []
    with pq.ParquetWriter(temp, schema, compression='zstd', compression_level=6) as writer:
        for batch in reader:
            rawids = batch.column(batch.schema.get_field_index(ids_col)).to_pylist()
            keep, models = [], []
            for index, value in enumerate(rawids):
                if value is None:
                    raise ValueError('Matrix has a missing row ID')
                if re.fullmatch(r'ACH-\d{6}', value):
                    model = value
                elif value in pmap:
                    if defaults and value not in defaults:
                        skipped_defaults += 1
                        continue
                    model = pmap[value]
                elif value in cmap:
                    model = cmap[value]
                else:
                    raise ValueError(f'Unrecognized matrix ID {value}; explicit mapping required')
                if not re.fullmatch(r'ACH-\d{6}', model):
                    raise ValueError(f'Invalid ModelID: {model}')
                if model in seen:
                    raise ValueError(f'Multiple rows map to {model}. Default profile/condition selection is required; no averaging performed.')
                seen.add(model)
                models.append(model)
                keep.append(index)
            if keep:
                arrays = [pa.array(models, type=pa.string())] + [batch.column(batch.schema.get_field_index(g)).take(pa.array(keep)) for g in genes]
                writer.write_table(pa.Table.from_arrays(arrays, schema=schema))
    if not seen:
        raise ValueError('No model rows survived mapping')
    qc = {**parquet_qc(temp), 'unique_models': len(seen), 'gene_columns': len(genes),
          'model_metadata_overlap': len(seen & model_ids), 'model_metadata_unmatched': len(seen - model_ids),
          'ModelID_format': 'passed', 'nondefault_profiles_excluded': skipped_defaults}
    if not (seen & model_ids):
        raise ValueError('Zero overlap with Model.csv')
    finish(temp, target, source, row['SHA256'], VERSION, {'qc': qc})
    gene_index = [{'column': g, 'symbol': re.sub(r'\s+\([^)]*\)$', '', g),
                   'entrez_id': (re.search(r'\((\d+)\)$', g)[1] if re.search(r'\((\d+)\)$', g) else None)} for g in genes]
    write_json(target.with_suffix('.genes.json'), gene_index)
    return qc


def main():
    setup('preprocess_depmap')
    release = current_release()
    folder = ROOT / 'data/processed/depmap' / release
    folder.mkdir(parents=True, exist_ok=True)
    qc, failures = {}, {}
    for filename, name in TABLES.items():
        try:
            raw = verified_raw('DepMap', release, filename)
            if raw:
                qc[name] = stream_table(*raw, folder / f'{name}.parquet')
            else:
                failures[name] = 'Raw unavailable: DepMap downloads skipped by user or require official browser verification'
        except Exception as e:
            logging.exception('Failed %s', filename)
            failures[name] = str(e)
    model_path = folder / 'Model.parquet'
    if model_path.exists():
        frame = pq.read_table(model_path).to_pandas()
        if 'ModelID' not in frame or frame.ModelID.duplicated().any() or not frame.ModelID.str.fullmatch(r'ACH-\d{6}').all():
            raise ValueError('Model.csv has missing, duplicate, or invalid ModelID')
        model_ids = set(frame.ModelID)
        maps = mapping(folder)
        for filename, name in WIDE.items():
            try:
                raw = verified_raw('DepMap', release, filename)
                if raw:
                    qc[name] = stream_matrix(*raw, folder / f'{name}.parquet', maps, model_ids)
                else:
                    failures[name] = 'Verified raw file unavailable (optional datasets may be absent in this release)'
            except Exception as e:
                logging.exception('Failed %s', filename)
                failures[name] = str(e)
    else:
        failures['model_overlap'] = 'Not measurable: Model.csv and all DepMap raw matrices are unavailable'
    write_json(ROOT / f'data/manifests/depmap_{release}_qc.json', {'created_at': now(), 'release': release,
               'qc': qc, 'failures': failures, 'status': 'complete' if not failures else 'incomplete'})
    print(json.dumps({'release': release, 'qc': qc, 'failures': failures}, indent=2))
    return bool(failures)


if __name__ == '__main__':
    raise SystemExit(main())

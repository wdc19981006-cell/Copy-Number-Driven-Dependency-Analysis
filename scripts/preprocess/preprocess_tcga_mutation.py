"""Prepare and validate the complete Xena MC3 table, without claiming it is MAF."""
from pathlib import Path
import gzip
import csv
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT, setup, verified_raw, write_json, now
from scripts.utils.preparation import finish, output_ready, parquet_qc
import pyarrow as pa
import pyarrow.csv as pcsv
import pyarrow.parquet as pq
import pyarrow.compute as pc

VERSION = 'mc3-xena-v1'


def main():
    setup('preprocess_tcga_mutation')
    raw = verified_raw('TCGA', 'PANCAN', 'mc3.v0.2.8.PUBLIC.xena.tsv.gz')
    if not raw:
        raise FileNotFoundError('Verified Xena MC3 table is unavailable')
    source, row = raw
    target = ROOT / 'data/processed/tcga/pancanatlas_reference/MC3.parquet'
    if output_ready(target, row['SHA256'], VERSION):
        print('MC3 derived file is verified; no conversion needed.')
        return
    with gzip.open(source, 'rt', encoding='utf-8-sig', newline='') as f:
        columns = next(csv.reader(f, delimiter='\t'))
    if not {'sample', 'gene', 'chr', 'start', 'end', 'effect'}.issubset(columns):
        raise ValueError('Unexpected MC3 Xena schema')
    reader = pcsv.open_csv(source, read_options=pcsv.ReadOptions(block_size=32 * 1024 * 1024),
                          parse_options=pcsv.ParseOptions(delimiter='\t'),
                          convert_options=pcsv.ConvertOptions(column_types={c: pa.string() for c in columns},
                                                             strings_can_be_null=True))
    temp = target.with_suffix('.parquet.part')
    sample_ids, gene_ids = set(), set()
    with pq.ParquetWriter(temp, reader.schema, compression='zstd', compression_level=6) as writer:
        for batch in reader:
            samples = batch.column(batch.schema.get_field_index('sample'))
            if samples.null_count or not pc.all(pc.match_substring_regex(samples, r'^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}')).as_py():
                raise ValueError('Invalid TCGA MC3 sample IDs')
            sample_ids.update(pc.unique(samples).to_pylist())
            gene_ids.update(pc.unique(batch.column(batch.schema.get_field_index('gene'))).to_pylist())
            writer.write_batch(batch)
    gene_ids.discard(None)
    qc = {**parquet_qc(temp), 'unique_samples': len(sample_ids), 'unique_genes': len(gene_ids),
          'parsed_all_rows': True, 'sample_id_format': 'passed',
          'source_filter': 'Xena retained source MAF FILTER=PASS only',
          'is_original_maf': False}
    finish(temp, target, source, row['SHA256'], VERSION, {'qc': qc})
    write_json(ROOT / 'data/manifests/tcga_mc3_qc.json', {'created_at': now(), 'qc': qc})
    print(qc)


if __name__ == '__main__':
    main()

"""Synthetic safety tests. They do not count as validation of real DepMap data."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.utils import common
from scripts.preprocess.preprocess_depmap import stream_matrix
from scripts.data_access import resolve_column
import pyarrow.parquet as pq
import hashlib


class InfrastructureTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='database_test_', dir=common.ROOT / 'results')
        self.folder = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_latest_release_uses_numeric_quarter_and_revision(self):
        rows = [{'release': s} for s in ['DepMap Public 25Q3', 'DepMap Public 26Q1',
                'DepMap Public 26Q3', 'DepMap Public 26Q3 v2', 'Harmonized Proteomics 27Q1']]
        self.assertEqual(common.latest_release(rows), 'DepMap Public 26Q3 v2')

    def test_model_mapping_and_projection_keep_all_genes(self):
        source = self.folder / 'matrix.csv'
        source.write_text(',GENEA (1),GENEB (2)\nPR-1,0.5,-0.9\nPR-2,1.1,NA\n', encoding='utf-8')
        target = self.folder / 'matrix.parquet'
        qc = stream_matrix(source, {'SHA256': common.hashes(source)[0]}, target,
                           ({'PR-1':'ACH-000001','PR-2':'ACH-000002'}, set(), {}),
                           {'ACH-000001','ACH-000002','ACH-000003'})
        self.assertEqual(qc['model_metadata_overlap'], 2)
        self.assertEqual(qc['gene_columns'], 2)
        self.assertEqual(resolve_column(target, 'GENEB'), 'GENEB (2)')
        self.assertEqual(resolve_column(target, '1'), 'GENEA (1)')
        table = pq.read_table(target, columns=['ModelID', 'GENEB (2)'])
        self.assertEqual(table.column_names, ['ModelID', 'GENEB (2)'])
        self.assertEqual(table['GENEB (2)'].null_count, 1)

    def test_ambiguous_model_mapping_fails_without_averaging(self):
        source = self.folder / 'ambiguous.csv'
        source.write_text(',G (1)\nPR-1,0\nPR-2,1\n', encoding='utf-8')
        target = self.folder / 'ambiguous.parquet'
        with self.assertRaisesRegex(ValueError, 'Multiple rows'):
            stream_matrix(source, {'SHA256': common.hashes(source)[0]}, target,
                          ({'PR-1':'ACH-000001','PR-2':'ACH-000001'}, set(), {}), {'ACH-000001'})
        self.assertFalse(target.exists())

    def test_default_profile_selection_is_explicit(self):
        source = self.folder / 'profiles.csv'
        source.write_text(',G (1)\nPR-1,0\nPR-2,1\n', encoding='utf-8')
        target = self.folder / 'profiles.parquet'
        qc = stream_matrix(source, {'SHA256': common.hashes(source)[0]}, target,
                          ({'PR-1':'ACH-000001','PR-2':'ACH-000001'}, {'PR-2'}, {}), {'ACH-000001'})
        self.assertEqual(qc['nondefault_profiles_excluded'], 1)
        self.assertEqual(pq.read_table(target)['G (1)'].to_pylist(), [1.0])

    def test_download_resumes_verifies_skips_and_preserves_corruption(self):
        payload = b'ModelID,G\nACH-000001,1\n'
        digest = hashlib.md5(payload).hexdigest()
        (self.folder / 'data/manifests').mkdir(parents=True)
        def fake_curl(args, **kwargs):
            output = Path(args[args.index('--output') + 1])
            # Simulate an actual resumed body starting after the local byte offset.
            offset = output.stat().st_size if output.exists() else 0
            with output.open('ab') as f:
                f.write(payload[offset:])
        with patch.object(common, 'ROOT', self.folder), \
             patch.object(common, 'MANIFEST', self.folder / 'data/manifests/data_manifest.csv'), \
             patch.object(common, 'remote_info', return_value={'content-length': str(len(payload)), 'etag':'test'}), \
             patch.object(common, 'curl', side_effect=fake_curl) as transfer:
            folder = self.folder / 'data/raw/depmap/26Q1'
            folder.mkdir(parents=True)
            part = folder / 'Model.csv.part'
            part.write_bytes(payload[:5])
            common.write_json(part.with_name(part.name + '.json'),
                              {'url':'https://example.org/data','etag':'test','length':len(payload)})
            path = common.download('DepMap','26Q1','Model.csv','Model.csv','https://example.org/data',digest)
            self.assertEqual(path.read_bytes(), payload)
            self.assertEqual(transfer.call_count, 1)
            common.download('DepMap','26Q1','Model.csv','Model.csv','https://example.org/data',digest)
            self.assertEqual(transfer.call_count, 1)
            path.write_bytes(b'changed')
            self.assertIsNone(common.download('DepMap','26Q1','Model.csv','Model.csv','https://example.org/data',digest))
            self.assertEqual(path.read_bytes(), b'changed')
            self.assertEqual(transfer.call_count, 1)


if __name__ == '__main__':
    unittest.main()

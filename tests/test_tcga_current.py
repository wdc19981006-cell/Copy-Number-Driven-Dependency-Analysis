"""Synthetic tests: per-sample modes, exact thresholds, RNA representatives."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.utils.tcga_current import classify_cn, sample_baselines, select_rna
import numpy as np
import pandas as pd


class CurrentTCGA(unittest.TestCase):
    def test_thresholds_and_missing(self):
        actual = classify_cn([0, 1, 3, 5, 6, np.nan, -1, 2], [3, 3, 3, 3, 3, 3, 3, 0])
        np.testing.assert_equal(actual, [-2, -1, 0, 1, 2, np.nan, np.nan, np.nan])

    def test_modes_exclude_sex_chromosomes_and_keep_sample_baselines(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            pd.DataFrame({'SampleID': ['s1', 's2', 's3'], 'g1': [2, 4, 1],
                          'g2': [2, 4, 2], 'g3': [3, 2, np.nan], 'gx': [9, 9, 9]}).to_parquet(folder/'TCGA_GeneLevel_CN.parquet')
            pd.DataFrame({'gene_id': ['g1', 'g2', 'g3', 'gx'], 'chromosome': ['chr1', '2', 'chr22', 'chrX']}).to_parquet(folder/'TCGA_GeneLevel_CN.genes.parquet')
            (folder/'TCGA_GeneLevel_CN.provenance.json').write_text(json.dumps({'release': '46.0', 'parsed_all_inputs': True}))
            result, report = sample_baselines(folder, folder/'cache', chunk_genes=1)
            self.assertEqual(result.BaselineCN.tolist(), [2, 4, 1])
            self.assertEqual(result.BaselineTies.tolist(), [1, 1, 2])
            self.assertEqual(result.AutosomalFiniteGenes.tolist(), [3, 3, 2])
            self.assertEqual(report['tied_samples'], 1)
            cached, _ = sample_baselines(folder, folder/'cache')
            pd.testing.assert_frame_equal(cached, result)

    def test_representative_is_independent_of_input_order(self):
        cn = pd.DataFrame({'SampleID': ['s1', 's2'], 'AliquotID': ['a2', 'b']})
        rna = pd.DataFrame({'SampleID': ['s1', 's1', 's2', 's2'],
                            'AliquotID': ['a1', 'a2', 'c', 'd'], 'FileID': ['f1', 'f9', 'f3', 'f2']})
        for frame in (rna, rna.iloc[::-1]):
            selected, audit = select_rna(cn, frame)
            self.assertEqual(selected.FileID.tolist(), ['f9', 'f2'])
            self.assertEqual(audit.groupby('SampleID').Selected.sum().tolist(), [1, 1])


if __name__ == '__main__':
    unittest.main()

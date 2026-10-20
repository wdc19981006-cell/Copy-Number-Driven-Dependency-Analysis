"""Offline fixtures for persisted identity, precision, selection, and cleanup gates."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.data_access import get_tcga_expression, DataUnavailableError, depmap_gene
from scripts.utils import tcga_current as tc
from scripts.utils import optimize_analysis_data as opt
from scripts.utils.common import hashes
from scripts.utils.validate_analysis_ready import compare_csvs
import numpy as np
import pandas as pd


class AnalysisReady(unittest.TestCase):
    def test_expression_preserves_double_missing_and_identifiers(self):
        values = pd.DataFrame({'SampleID':['a','b','c','d'],
                               'Expression':[0., np.nan, 1.23456789012345, 1e-16]})
        with patch('scripts.data_access.gdc_gene', return_value=values.copy()) as read:
            result = get_tcga_expression('ARBITRARY_GENE')
        read.assert_called_once_with('TCGA_STAR_TPM', 'ARBITRARY_GENE', 'Expression')
        self.assertEqual(result.SampleID.tolist(), values.SampleID.tolist())
        np.testing.assert_array_equal(result.Expression, np.log2(values.Expression+1))
        self.assertEqual(str(result.Expression.dtype), 'float64')

    def test_prepared_table_tamper_and_stale_input_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            source = folder/'input.parquet';source.write_bytes(b'original input')
            path = folder/'TCGA_Sample_Baselines.parquet'
            pd.DataFrame({'SampleID':['s'], 'BaselineCN':[4]}).to_parquet(path)
            meta = dict(release='46.0', algorithm_version='analysis-ready-v1',
                        SHA256=hashes(path)[0], inputs={'input.parquet':opt.signature(source)})
            path.with_suffix('.provenance.json').write_text(json.dumps(meta))
            result, _ = tc.prepared_table('TCGA_Sample_Baselines', folder)
            self.assertEqual(result.BaselineCN.tolist(), [4])
            source.write_bytes(b'changed input')
            with self.assertRaisesRegex(ValueError, 'Stale'):
                tc.prepared_table('TCGA_Sample_Baselines', folder)
            path.write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
                tc.prepared_table('TCGA_Sample_Baselines', folder)

    def test_persisted_selection_retains_all_files_and_missing_expression(self):
        cn = pd.DataFrame({'SampleID':['s1'], 'AliquotID':['a2']})
        rna = pd.DataFrame({'SampleID':['s1','s1','s2'], 'AliquotID':['a1','a2','b'],
                            'FileID':['f1','f9','f2'], 'CaseID':['c1','c1','c2'],
                            'ProjectID':['TCGA-ACC']*3, 'TPM':[7.,np.nan,0.]})
        direct, audit = tc.select_rna(cn, rna)
        saved = audit.drop(columns='TPM').reset_index(drop=True)
        with patch.object(tc, 'prepared_table', return_value=(saved, {})):
            selected, actual = tc.selected_rna(cn, rna.iloc[::-1])
            pd.testing.assert_frame_equal(actual, audit.reset_index(drop=True))
            self.assertEqual(selected.FileID.tolist(), direct.FileID.tolist())
            self.assertTrue(selected.TPM.isna().any())
            with self.assertRaisesRegex(ValueError, 'coverage'):
                tc.selected_rna(cn, rna.iloc[:1])

    def test_optional_data_absence_never_substitutes(self):
        with tempfile.TemporaryDirectory() as tmp, patch('scripts.data_access.depmap_path', return_value=Path(tmp)/'missing.parquet'):
            with self.assertRaises(DataUnavailableError):
                depmap_gene('ExpressionProteinCoding', 'GENE', 'Expression')

    def test_cleanup_cannot_delete_code_or_traverse_outside_data(self):
        with self.assertRaises(ValueError):
            opt.safe_target('scripts/data_access.py')
        with self.assertRaises(ValueError):
            opt.safe_target('data/../scripts/data_access.py')

    def test_real_csv_comparison_handles_boolean_nulls_and_rejects_tiny_p_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp)/'before';actual=Path(tmp)/'after'
            for folder in [base,actual]:(folder/'Main_Results').mkdir(parents=True)
            frame=pd.DataFrame({'SampleID':['s1','s2'],'Selected':[True,False],
                                'P':[1e-180,np.nan],'Eligible_Rank':[1.,np.nan]})
            for i in range(20):
                for folder in [base,actual]:frame.to_csv(folder/f'Main_Results/{i}.csv',index=False)
            compared,error=compare_csvs(base,actual)
            self.assertEqual(len(compared),20);self.assertEqual(error,0)
            changed=frame.copy();changed.loc[0,'P']=0
            changed.to_csv(actual/'Main_Results/0.csv',index=False)
            with self.assertRaises(AssertionError):compare_csvs(base,actual)

    def test_worker_guard_distinguishes_idle_external_mcp_from_project_jobs(self):
        proc=SimpleNamespace(pid=999999,info={'name':'python.exe','cmdline':[
            'C:/Python312/python.exe','D:/CodexProjects/geo-mcp-server/server.py']},
            cwd=lambda:str(ROOT),open_files=lambda:[])
        with patch('psutil.process_iter',return_value=[proc]):opt.no_workers()
        proc.open_files=lambda:[SimpleNamespace(path=str(ROOT/'data/raw/active.tsv'))]
        with patch('psutil.process_iter',return_value=[proc]):
            with self.assertRaisesRegex(RuntimeError,'Active project worker'):opt.no_workers()
        proc.open_files=lambda:[]
        proc.info['cmdline']=['C:/Python312/python.exe',str(ROOT/'scripts/R/run_analysis.R')]
        with patch('psutil.process_iter',return_value=[proc]):
            with self.assertRaisesRegex(RuntimeError,'Active project worker'):opt.no_workers()


if __name__ == '__main__':
    unittest.main()

import importlib.util
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

spec = importlib.util.spec_from_file_location('guard', Path(__file__).resolve().parents[1] / 'scripts/utils/git_size_guard.py')
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)


class GuardTests(unittest.TestCase):
    def test_all_cases_and_exclusions(self):
        for case in ['VPS4B_VPS4A_Analysis', 'SMARCA4_Analysis', 'ENO1_Analysis']:
            self.assertFalse(guard.forbidden(f'results/{case}/Main_Results/a.pdf'))
        for name in ['data/raw/a.csv', 'data/processed/a.csv', '.runtime/a.txt',
                     'logs/a.txt', 'results/CASE/Tables/cache/a.csv', 'a.parquet', 'a.duckdb']:
            self.assertTrue(guard.forbidden(name), name)

    def test_refuse_oversize_before_stage(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            path = root / 'large.csv'
            with path.open('wb') as out:
                out.truncate(guard.LIMIT + 1)
            self.assertEqual(len(guard.inspect_worktree(['large.csv'], root)), 1)
            self.assertFalse((root / '.git').exists())

    def test_actual_index_blob_even_if_worktree_changed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            subprocess.run(['git', 'init', '-q', str(root)], check=True)
            (root / 'a.txt').write_text('small', encoding='utf-8')
            guard.git('add', 'a.txt', cwd=root)
            (root / 'a.txt').write_text('different', encoding='utf-8')
            report = guard.inspect_index(root)
            self.assertEqual(report['largest_files'][0][0], 5)
            self.assertEqual(report['issues'], [])

    def test_root_alias_keeps_size_and_containment_checks(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / 'project'
            (root / 'step').mkdir(parents=True)
            alias = root / 'step' / '..'
            (root / 'large.csv').write_bytes(b'12345')
            (root.parent / 'outside.csv').write_bytes(b'outside')
            with patch.object(guard, 'LIMIT', 4):
                issues = guard.inspect_worktree(['large.csv'], alias)
            self.assertEqual(len(issues), 1)
            self.assertIn('larger than', issues[0][1])
            self.assertEqual(guard.inspect_worktree(['../outside.csv'], alias),
                             [('../outside.csv', 'outside project')])


if __name__ == '__main__':
    unittest.main()

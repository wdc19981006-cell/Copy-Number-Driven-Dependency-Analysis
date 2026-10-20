"""Run public workflows in a core-only offline fixture and compare saved CSVs.

The fixture contains hard links to immutable core Parquets, copies of code, and
the installed R/Python libraries. Neither raw nor retired processed files exist.
No source result is modified. Evidence is consumed by the separate cleanup gate.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.utils.optimize_analysis_data import CORE, WORK, MANIFESTS, REPORT
from scripts.utils.common import hashes, write_json, now
import numpy as np
import pandas as pd


def hardlink_tree(source, target):
    shutil.copytree(source, target, copy_function=os.link,
                    ignore=shutil.ignore_patterns('__pycache__', '*.pyc', 'staging'))


def create_fixture(resume=False):
    target = WORK/'core_only'
    if target.exists() and not resume:
        raise FileExistsError('Core-only fixture exists; inspect before resuming')
    if target.exists():
        if (target/'data/raw').exists():
            raise ValueError('Resume fixture contains forbidden raw data')
        expected={p.relative_to(ROOT).as_posix() for p in CORE}
        actual={p.relative_to(target).as_posix() for p in (target/'data/processed').rglob('*.parquet')}
        if expected!=actual:
            raise ValueError('Resume fixture contains missing/extra processed data')
        for name in ['scripts','tests','config','docs']:
            shutil.copytree(ROOT/name,target/name,dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__'))
        for name in ['duckdb-1.5.6.dist-info','pyarrow-25.0.1.dist-info']:
            if not (target/'.runtime'/name).exists():
                hardlink_tree(ROOT/'.runtime'/name,target/'.runtime'/name)
        return target
    target.mkdir()
    for name in ['scripts', 'tests', 'config', 'docs']:
        shutil.copytree(ROOT/name, target/name, ignore=shutil.ignore_patterns('__pycache__'))
    for name in ['.Rprofile', 'renv.lock', 'DESCRIPTION', '.R-version']:
        shutil.copy2(ROOT/name, target/name)
    hardlink_tree(ROOT/'renv', target/'renv')
    runtime = target/'.runtime';runtime.mkdir()
    for name in ['pyarrow', 'pyarrow.libs', 'duckdb', '_duckdb.cp312-win_amd64.pyd',
                 'duckdb-1.5.6.dist-info','pyarrow-25.0.1.dist-info']:
        src = ROOT/'.runtime'/name;dest=runtime/name
        if src.is_dir():
            hardlink_tree(src,dest)
        else:
            os.link(src,dest)
    (target/'results').mkdir()
    (target/'data/manifests').mkdir(parents=True)
    shutil.copy2(ROOT/'data/manifests/depmap_26Q1_manifest.csv',target/'data/manifests/depmap_26Q1_manifest.csv')
    for path in CORE:
        dest=target/path.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True)
        os.link(path,dest)
        for side in [path.with_suffix('.provenance.json'), Path(str(path)+'.provenance.json')]:
            if side.exists():
                shutil.copy2(side,target/side.relative_to(ROOT))
    guard = runtime/'offline_guard';guard.mkdir()
    (guard/'sitecustomize.py').write_text('''import sys, os
from pathlib import Path
log = Path(os.environ['ANALYSIS_ACCESS_LOG'])
with log.open('a') as f: f.write('PYTHON_OFFLINE_GUARD_STARTED\\n')
def audit(event, args):
    if event in ('socket.connect', 'socket.getaddrinfo', 'socket.bind'):
        with log.open('a') as f: f.write('DENIED_NETWORK '+event+'\\n')
        raise RuntimeError('Offline regression forbids network')
    if event == 'open' and args and isinstance(args[0], (str, bytes)):
        path = os.fsdecode(args[0]).replace('\\\\','/').lower()
        if '/data/raw/' in '/'+path or path.startswith('data/raw/'):
            with log.open('a') as f: f.write('DENIED_RAW '+path+'\\n')
            raise RuntimeError('Offline regression forbids raw data')
sys.addaudithook(audit)
''',encoding='utf-8')
    (runtime/'offline_run.R').write_text('''deny_network <- function(...)stop("Offline regression forbids network")
for(n in c("url","socketConnection"))trace(n,where=baseenv(),tracer=quote(deny_network()),print=FALSE)
trace("download.file",where=asNamespace("utils"),tracer=quote(deny_network()),print=FALSE)
stopifnot(!dir.exists("data/raw"))
source("scripts/R/run_analysis.R",encoding="UTF-8")
''',encoding='utf-8')
    return target


def run_logged(command, cwd, log, env):
    with log.open('w',encoding='utf-8') as out:
        result=subprocess.run(command,cwd=cwd,stdout=out,stderr=subprocess.STDOUT,env=env)
    if result.returncode:
        raise RuntimeError(f'Validation failed ({result.returncode}): {command}; see {log}')


def compare_csvs(reference, actual):
    compared=[];max_error=0.0
    paths=[p for p in reference.rglob('*.csv') if any(
        p.relative_to(reference).as_posix().startswith(prefix) for prefix in [
            'Main_Results/', 'Supplementary/', 'Provenance/Data_Audit/'])]
    for path in sorted(paths):
        rel=path.relative_to(reference)
        before=pd.read_csv(path, low_memory=False)
        after=pd.read_csv(actual/rel, low_memory=False)
        # Preserve row order and identifiers; only CSV double serialization admits
        # relative machine-precision tolerance (never absolute tolerance for P).
        if list(before.columns)!=list(after.columns) or before.shape!=after.shape:
            raise AssertionError(f'CSV shape/columns changed: {rel}')
        for column in before.columns:
            a=before[column];b=after[column]
            if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
                np.testing.assert_allclose(a,b,rtol=5e-13,atol=0,equal_nan=True,err_msg=f'{rel}: {column}')
                use=np.isfinite(a)&np.isfinite(b)
                if use.any():
                    max_error=max(max_error,float(np.max(np.abs(a[use].to_numpy(dtype=float)-b[use].to_numpy(dtype=float)))))
            else:
                pd.testing.assert_series_equal(a,b,check_dtype=False,check_names=False)
        compared.append(dict(path=rel.as_posix(),rows=len(before),columns=len(before.columns),
                             reference_SHA256=hashes(path)[0],actual_SHA256=hashes(actual/rel)[0]))
    if len(compared)<20:
        raise AssertionError('Incomplete real-case reference')
    return compared,max_error


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geneA',default='MTAP');parser.add_argument('--geneB',default='WDR77')
    parser.add_argument('--resume',action='store_true',help='Resume an inspected core-only fixture after an environment repair')
    parser.add_argument('--verify-saved',action='store_true',help='Validate completed isolated runs without executing workflows again')
    args=parser.parse_args()
    prepared=json.loads(REPORT.read_text())
    for name,sig in prepared['core'].items():
        if hashes(ROOT/name)[0]!=sig['SHA256']:raise ValueError('Prepared core hash changed')
    for log in ['synthetic_R.log','synthetic_Python.log']:
        content=(WORK/log).read_text(encoding='utf-8',errors='replace')
        if 'FAIL' in content or ('PASS:' not in content if log.startswith('synthetic_R') else '\nOK' not in content):
            raise ValueError('Synthetic tests not passed: '+log)
    fixture=WORK/'core_only' if args.verify_saved else create_fixture(args.resume)
    if args.verify_saved and (not fixture.is_dir() or (fixture/'data/raw').exists()):
        raise ValueError('Invalid completed core-only fixture')
    env=os.environ.copy()
    env.update(DATA_MODE='local',PYTHONPATH=str(fixture/'.runtime/offline_guard'),
               ANALYSIS_ACCESS_LOG=str(WORK/'offline_access.log'),RENV_CONFIG_AUTOLOADER_ENABLED='TRUE')
    rscript=json.loads((ROOT/'config/R_runtime.json').read_text())['rscript']
    entry=str(fixture/'.runtime/offline_run.R')
    for workflow in ['geneA_screen','geneA_geneB']:
        cmd=[rscript,'--vanilla',entry,'--workflow',workflow,'--geneA',args.geneA,'--project',str(fixture)]
        if workflow=='geneA_geneB':cmd+=['--geneB',args.geneB]
        if args.verify_saved:
            if 'Analysis completed:' not in (WORK/(workflow+'.log')).read_text(encoding='utf-8',errors='replace'):
                raise ValueError('Isolated workflow has not completed: '+workflow)
        else:
            print('Running isolated',workflow,flush=True)
            run_logged(cmd,fixture,WORK/(workflow+'.log'),env)
    reference=ROOT/'results'/f'{args.geneA}_Analysis'
    actual=fixture/'results'/reference.name
    compared,error=compare_csvs(reference,actual)
    loaded_R=['run_analysis.R','data_access.R','plotting.R','dependency_screen.R','lineage_analysis.R',
              'adapters.R','additional_modes.R','tcga_analysis.R','extensions_statistics.R',
              'extension_modules.R','extension_summary.R','local_data.R','workflow_statistics.R',
              'workflow_plots.R','final_workflow.R','geneA_project.R','validate_final_workflow.R']
    workflow_sources=[ROOT/'scripts/R'/n for n in loaded_R]+[ROOT/p for p in [
        'scripts/data_access.py','scripts/utils/tcga_current.py','scripts/utils/export_workflow_tcga.py']]
    for source in workflow_sources:
        if hashes(source)[0]!=hashes(fixture/source.relative_to(ROOT))[0]:
            raise ValueError(f'Workflow source changed during regression: {source.name}')
    # The second workflow must reuse every A module, proving validated caching.
    runs=pd.read_csv(actual/f'Provenance/GeneB/{args.geneB}/Module_Runs.csv')
    assert runs.loc[runs.Module.isin(['tcga','genomewide_dependency','cn_covariation']),'Execution'].eq('cached').all()
    run_logged([sys.executable,str(fixture/'scripts/utils/validate_final_sources.py'),'--case',reference.name],
               fixture,WORK/'independent_sources.log',env)
    access=(WORK/'offline_access.log').read_text()
    assert 'PYTHON_OFFLINE_GUARD_STARTED' in access and 'DENIED_' not in access
    assert not (fixture/'data/raw').exists()
    for folder in ['Main_Results/03_TCGA_CN_mRNA']:
        assert len(list((actual/folder).glob('*.pdf')))==33
    evidence=dict(status='PASS',validated_at=now(),geneA=args.geneA,geneB=args.geneB,
        synthetic_R=True,synthetic_Python=True,isolated_geneA_screen=True,isolated_geneA_geneB=True,
        all_saved_CSVs_equal=True,no_raw_or_network=True,independent_source_check=True,
        base_modules_reused=True,core_only_fixture='.runtime/database_optimization/core_only',
        comparison_tolerance='rtol=5e-13, atol=0; missing masks/order/identifiers must match',
        maximum_absolute_CSV_difference=error,compared_CSVs=compared,
        code_sha256={p.relative_to(ROOT).as_posix():hashes(p)[0] for p in sorted(
            list((ROOT/'scripts').rglob('*.py'))+list((ROOT/'scripts/R').glob('*.R')))},
        logs={p.name:hashes(p)[0] for p in WORK.glob('*.log')})
    write_json(MANIFESTS/'analysis_ready_regression.json',evidence)
    print('PASS: isolated workflows and',len(compared),'real-case CSVs',flush=True)


if __name__=='__main__':main()

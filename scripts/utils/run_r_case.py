"""Run modes sequentially with exact R and monitor process-tree RSS."""
from pathlib import Path
import argparse,json,subprocess,time,sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'.runtime'))
import psutil
RSCRIPT=Path('D:/R/R-4.5.0/bin/Rscript.exe')
MODES=['qc','depmap_cn_expression','genomewide_dependency','targeted_dependency','lineage_dependency','adjusted_dependency','reverse_dependency','cn_covariation','mutation_dependency','tcga_cn_landscape','tcga_cna_prevalence','tcga_cn_expression','full']

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--geneA',default='VPS4B');parser.add_argument('--geneB',default='VPS4A');parser.add_argument('--modes',nargs='+',default=MODES);parser.add_argument('--output_case');args=parser.parse_args()
 case=args.output_case or f'{args.geneA}_{args.geneB}'
 if any(v in case for v in ['/','\\']) or case in {'.','..'}:raise ValueError('Invalid case folder')
 report=ROOT/f'results/{case}/Summary/Resource_Monitor.json'
 reports=json.loads(report.read_text(encoding='utf-8')) if report.exists() else []
 for mode in args.modes:
  log=ROOT/f'logs/{args.geneA}_{args.geneB}_{mode}.log';start=time.monotonic();peak=0
  with log.open('w',encoding='utf-8') as out:
   command=[str(RSCRIPT),'--vanilla',str(ROOT/'scripts/R/run_analysis.R'),'--project',str(ROOT),'--mode',mode,'--geneA',args.geneA,'--geneB',args.geneB]
   if args.output_case:command.extend(['--output_case',args.output_case])
   proc=subprocess.Popen(command,cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW if sys.platform=='win32' else 0)
   ps=psutil.Process(proc.pid)
   while proc.poll() is None:
    try:peak=max(peak,sum(p.memory_info().rss for p in [ps,*ps.children(recursive=True)] if p.is_running()))
    except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    time.sleep(.5)
  row=dict(mode=mode,exit_code=proc.returncode,wall_seconds=round(time.monotonic()-start,3),approx_peak_RSS_MiB=round(peak/2**20,2),log=str(log.relative_to(ROOT)))
  reports.append(row);print(json.dumps(row),flush=True)
  report.parent.mkdir(parents=True,exist_ok=True);report.write_text(json.dumps(reports,indent=2),encoding='utf-8')
  if proc.returncode:raise SystemExit(f'Failed {mode}; inspect {log}')

if __name__=='__main__':main()

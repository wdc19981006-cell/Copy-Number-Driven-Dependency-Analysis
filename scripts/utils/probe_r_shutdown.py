"""Isolate shutdown failure; report real native process exit codes."""
from pathlib import Path
import subprocess,json
ROOT=Path(__file__).resolve().parents[2]
script=ROOT/'scripts/R/probe_package.R'
script.write_text('Sys.setlocale("LC_CTYPE","English_United States.utf8");source(".Rprofile");args<-commandArgs(TRUE);for(p in args)library(p,character.only=TRUE);cat("PROBE FINISHED\\n")\n',encoding='utf-8')
results=[]
for pkg in ['rlang','vctrs','S7','cli','tibble','tidyselect','glue','magrittr','pillar']:
 result=subprocess.run(['D:/R/R-4.5.0/bin/Rscript.exe','--vanilla',str(script),pkg],cwd=ROOT,capture_output=True)
 row=dict(package=pkg,exit_code=result.returncode,stderr_tail=result.stderr.decode('utf-8',errors='replace')[-500:]);results.append(row);print(json.dumps(row),flush=True)
(ROOT/'docs/R_SHUTDOWN_DIAGNOSTICS.json').write_text(json.dumps(results,indent=2),encoding='utf-8')

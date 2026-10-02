from pathlib import Path
import subprocess,json,os
ROOT=Path(__file__).resolve().parents[2]
for utf in [False,True]:
 for local in [False,True]:
  code=('Sys.setlocale("LC_CTYPE","English_United States.utf8");' if utf else 'Sys.setlocale("LC_CTYPE","C");')
  if local:code+=' .libPaths(c(file.path(getwd(),"renv/library/windows/R-4.5/x86_64-w64-mingw32"),.libPaths()));'
  code+='library(cli);cat("Finished\\n")'
  args=['D:/R/R-4.5.0/bin/Rscript.exe','--vanilla','-e',code]
  result=subprocess.run(args,cwd=ROOT,capture_output=True)
  print(json.dumps({'utf8':utf,'project_lib':local,'exit':result.returncode,'tail':result.stderr.decode('utf-8',errors='replace')[-150:]}),flush=True)

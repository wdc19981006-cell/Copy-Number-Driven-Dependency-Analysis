from pathlib import Path
import os,subprocess,json
ROOT=Path(__file__).resolve().parents[2]
env=os.environ.copy()
env.update(LANG='English_United States.utf8',LC_ALL='English_United States.utf8',LC_CTYPE='English_United States.utf8')
env.pop('R_HOME',None)
exe=ROOT/'tools/R450/runtime/app/bin/Rscript.exe'
for label,executable in [('original',Path('D:/R/R-4.5.0/bin/Rscript.exe')),('clean',exe)]:
 code='.libPaths(c(file.path(getwd(),"renv/library/windows/R-4.5/x86_64-w64-mingw32"),.libPaths()));library(rlang);library(cli);print(R.version.string);print(R.home());print(warnings())'
 p=subprocess.run([str(executable),'--vanilla','-e',code],env=env,cwd=ROOT,capture_output=True)
 print(json.dumps({'runtime':label,'exit':p.returncode,'output':p.stdout.decode('utf-8',errors='replace'),'stderr':p.stderr.decode('utf-8',errors='replace')}),flush=True)

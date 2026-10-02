from pathlib import Path
import subprocess,tarfile
ROOT=Path(__file__).resolve().parents[2]
dest=ROOT/'.runtime/cli-source';dest.mkdir(exist_ok=True)
with tarfile.open(ROOT/'.runtime/compatible-src/cli_3.6.4.tar.gz') as t:
 t.extractall(dest,filter='data')
print('CLI native source inspected in project runtime folder')

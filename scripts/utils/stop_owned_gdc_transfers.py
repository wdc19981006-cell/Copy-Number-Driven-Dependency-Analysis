"""Stop only this project's transfer processes so verified files can be resumed."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,curl_binary
import psutil

targets={(ROOT/'scripts/download/download_gdc.py').resolve(),(ROOT/'scripts/download/download_gdc_api.py').resolve(),(ROOT/'scripts/download/download_gdc_bundles.py').resolve()}
toolroot=(ROOT/'tools/gdc-client').resolve()
curlpath=Path(curl_binary()).resolve()
console_host=Path('C:/Windows/System32/conhost.exe').resolve()
for process in psutil.process_iter(['pid','cmdline']):
    try:
        args=process.info['cmdline'] or []
        if not any(Path(arg).resolve() in targets for arg in args if 'download_gdc' in arg):
            continue
        children=process.children(recursive=True)
        safe=[]
        for child in children:
            if Path(child.exe()).resolve().is_relative_to(toolroot) or Path(child.exe()).resolve() in {curlpath,console_host}:
                safe.append(child)
            else:
                raise ValueError(f'Unexpected child process {child.pid}; stopping requires explicit review')
        print('Stopping owned transfer coordinator',process.pid,'and',len(safe),'verified transfer/console children')
        process.terminate()
        for child in safe:
            try: child.terminate()
            except psutil.NoSuchProcess: pass
        psutil.wait_procs([process,*safe],timeout=10)
    except (psutil.NoSuchProcess,psutil.AccessDenied):
        continue

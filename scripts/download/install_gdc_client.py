"""Install official verified GDC executable locally inside this project only."""
from pathlib import Path
import zipfile
import sys
import subprocess
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT,setup,curl,hashes,write_json,now

URL='https://gdc.cancer.gov/system/files/public/file/gdc-client_2.3_Windows_x64-py3.8-windows-2019.zip'
MD5='525ce44bb5f3f0624066b906c7dbdaf4'


def main():
    setup('install_gdc_client')
    folder=ROOT/'tools/gdc-client'
    folder.mkdir(parents=True,exist_ok=True)
    package=folder/'gdc-client_2.3_Windows_x64-py3.8-windows-2019.zip'
    if not package.exists():
        curl(['--retry','3','--output',package,URL])
    sha,md5=hashes(package)
    if md5!=MD5:
        raise ValueError('Official installer MD5 mismatch; executable will not be run')
    distribution=folder/'distribution'
    distribution.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(package) as z:
        for member in z.infolist():
            target=(distribution/member.filename).resolve()
            if not target.is_relative_to(distribution.resolve()):
                raise ValueError('Unsafe zip path')
        z.extractall(distribution)
    for nested in distribution.rglob('*.zip'):
        nested_target=nested.parent/'contents'
        with zipfile.ZipFile(nested) as z:
            if any(not (nested_target/m.filename).resolve().is_relative_to(nested_target.resolve()) for m in z.infolist()):
                raise ValueError('Unsafe nested zip path')
            z.extractall(nested_target)
    executable=next(folder.rglob('gdc-client.exe'))
    write_json(ROOT/'config/gdc_client.json',{'executable':executable.relative_to(ROOT).as_posix(),
               'download_url':URL,'publisher_md5':MD5,'SHA256':sha,'installed_at':now()})
    subprocess.run([str(executable),'download','--help'],check=True)


if __name__=='__main__':
    main()

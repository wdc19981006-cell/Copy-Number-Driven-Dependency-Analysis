"""Explicit database maintenance. Never imported or called by ordinary analyses."""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMMANDS = {
    'update_tcga': ['scripts/download/discover_gdc.py', 'scripts/download/download_gdc.py',
                    'scripts/preprocess/preprocess_gdc.py'],
    'update_depmap': ['scripts/download/download_depmap.py', 'scripts/preprocess/preprocess_depmap.py'],
    'prepare_tcga_reference': ['scripts/preprocess/preprocess_tcga.py'],
    'prepare_depmap_local': ['scripts/preprocess/preprocess_local_depmap.py'],
    'validate_data': ['scripts/utils/validate_database.py'],
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=COMMANDS)
    args = parser.parse_args()
    for script in COMMANDS[args.command]:
        subprocess.run([sys.executable, str(ROOT / script)], cwd=ROOT, check=True)


if __name__ == '__main__':
    main()

from pathlib import Path
import json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT, get_catalog, latest_release, setup


def main():
    setup('check_updates')
    available = latest_release(get_catalog()).removeprefix('DepMap Public ')
    path = ROOT / 'config/current_release.json'
    local = json.loads(path.read_text(encoding='utf-8'))['depmap_release'] if path.exists() else None
    if local == available:
        print('DepMap database is up to date.')
        print('Release selection matches the catalog; this does not assert download completeness.')
    else:
        print(f'Current local: {local or "none"}')
        print(f'New DepMap release available: {available}')
        print('Existing databases were not changed. Explicit selection is required to install another release.')


if __name__ == '__main__':
    main()

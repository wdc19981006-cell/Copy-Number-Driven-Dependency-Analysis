"""Render final PDFs to ignored QA contact sheets; never recompute analyses."""
import argparse
import subprocess
from pathlib import Path
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[2]
POPPLER = Path.home() / '.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', default='VPS4B_VPS4A_Analysis')
    a = p.parse_args()
    base = ROOT / 'results' / a.case / 'Main_Results'
    out = ROOT / '.runtime/qa/final_workflow'
    out.mkdir(parents=True, exist_ok=True)
    paths = sorted(base.rglob('*.pdf'))
    thumbnails = []
    for i, path in enumerate(paths):
        prefix = out / f'page_{i:02d}'
        subprocess.run([str(POPPLER), '-singlefile', '-scale-to', '1050', '-png',
                        str(path), str(prefix)], check=True, stdout=subprocess.DEVNULL)
        image = Image.open(prefix.with_suffix('.png')).convert('RGB')
        thumbnails.append((path.name, image.copy()))
    for start in range(0, len(thumbnails), 6):
        canvas = Image.new('RGB', (1600, 1860), 'white')
        draw = ImageDraw.Draw(canvas)
        for j, (name, image) in enumerate(thumbnails[start:start+6]):
            image.thumbnail((780, 565))
            x, y = (j % 2) * 800, (j // 2) * 620
            draw.text((x + 12, y + 8), name, fill='black')
            canvas.paste(image, (x + (800-image.width)//2, y + 38))
        target = out / f'contact_{start//6:02d}.png'
        canvas.save(target)
        print(target.as_posix())


if __name__ == '__main__':
    main()

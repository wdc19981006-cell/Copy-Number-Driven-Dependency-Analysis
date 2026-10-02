"""Measure indexed gene retrieval only; no statistical or biological analysis."""
from pathlib import Path
import sys
import time
import json
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.utils.common import ROOT, write_json, setup
from scripts.data_access import get_tcga_cn, get_tcga_expression


def main():
    setup('benchmark_queries')
    results = []
    for query in [get_tcga_expression, get_tcga_cn]:
        for gene in ['TP53', 'EGFR']:
            begin = time.perf_counter()
            frame = query(gene, layer='reference')
            results.append({'dataset': query.__name__, 'gene': gene, 'samples': len(frame),
                            'seconds': round(time.perf_counter() - begin, 3)})
    write_json(ROOT / 'data/manifests/query_timing.json', results)
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()

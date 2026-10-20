# Manifest scopes

- `depmap_26Q1_manifest.csv` / `depmap_26Q1_local_qc.json`: active user-provided exports, unchanged raw SHA256, dimensions and full supplied-column preprocessing.
- `gdc_data_manifest.csv`, `gdc_download_summary.json`, `gdc_*_qc.json`: dated DR46 acquisition/preparation snapshots. Former live source manifests are now archived under `historical_DR46/`; raw has been deleted. `gdc_progress.py` / `publish_gdc_snapshot.py` belong to separate acquisition maintenance, not current analysis readiness.
- `current_missing_data.csv`: outstanding sources/coverage at the snapshot, separate from historical initial skips.
- `data_manifest.csv`, `not_downloaded.csv`, `preparation_summary.json`, `depmap_26Q1_qc.json`, `validation_report.json`: historical initial official-download/preparation stage. Official full DepMap downloads were skipped then; later local exports are active and are described in their own manifest. These historical files must not be used to claim the current local exports are absent.
- `*.source.json`, `tcga_qc.json`, `tcga_mc3_qc.json`: reference Xena objects, scales and preparation.
- `R_package_versions.csv`: actual project package versions, with `renv.lock` and `docs/RENV_STATUS.txt`.

Current interpretation and source limitations: [DATA_SOURCES](../../docs/DATA_SOURCES.md). Processed databases remain local and are excluded from Git.

## Analysis-ready storage

`analysis_ready_validation.json`, `analysis_ready_regression.json`, `data_cleanup_report.json` and the exact JSONL deletion inventories describe the current minimal DR46/26Q1 database and actual cleanup. See [storage contract](../../docs/ANALYSIS_READY_STORAGE.md).

After cleanup, prior raw integrity/download/QC records are **historical verification records**, not evidence that those raw files still exist. Source file IDs, original MD5/SHA256, preprocessing versions and coverage are retained; the `historical_DR46/` snapshot preserves the former live manifests. Only the analysis-ready report defines current processed availability. Do not run historical raw validators as default analysis checks.

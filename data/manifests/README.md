# Manifest scopes

- `depmap_26Q1_manifest.csv` / `depmap_26Q1_local_qc.json`: active user-provided exports, unchanged raw SHA256, dimensions and full supplied-column preprocessing.
- `gdc_data_manifest.csv`, `gdc_download_summary.json`, `gdc_*_qc.json`: explicitly dated current DR46 acquisition/preparation snapshot. Actual local progress is in ignored `data/raw/tcga/gdc_current_DR46/manifests/live/`; use `scripts/utils/gdc_progress.py` to inspect it and `publish_gdc_snapshot.py` to publish a new snapshot deliberately.
- `current_missing_data.csv`: outstanding sources/coverage at the snapshot, separate from historical initial skips.
- `data_manifest.csv`, `not_downloaded.csv`, `preparation_summary.json`, `depmap_26Q1_qc.json`, `validation_report.json`: historical initial official-download/preparation stage. Official full DepMap downloads were skipped then; later local exports are active and are described in their own manifest. These historical files must not be used to claim the current local exports are absent.
- `*.source.json`, `tcga_qc.json`, `tcga_mc3_qc.json`: reference Xena objects, scales and preparation.
- `R_package_versions.csv`: actual project package versions, with `renv.lock` and `docs/RENV_STATUS.txt`.

Current interpretation and source limitations: [DATA_SOURCES](../../docs/DATA_SOURCES.md). Raw matrices, metadata query pages and processed databases remain local and are excluded from Git.

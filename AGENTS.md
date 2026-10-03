# Repository workflow rules

- Windows / Git Bash. Use `/c/Python312/python.exe`; never the Windows Store `python`/`python3`. Pass Python script paths as Windows-native absolute paths. R must be 4.5.0; runtime location is in `config/R_runtime.json`.
- Start tasks with `git pull --ff-only origin main`, preserving existing local work.
- Ordinary analyses use `DATA_MODE=local`, existing processed Parquet/DuckDB only. Do not download, update, read raw or silently substitute data layers. Missing processed data must fail with the two-line message in `docs/FINAL_WORKFLOWS.md`.
- User workflows: `geneA_screen` and `geneA_geneB`. Do not automatically add modules. Hidden legacy modes remain explicit developer calls. Avoid repeated full/screen/sensitivity runs; use validated per-module caches, and only rerun affected modules after a fix.
- TCGA main results use the explicitly documented reference layer. Keep NCI GDC DR46 and PanCanAtlas/Xena reference provenance separate; user-facing names say TCGA. Keep Chronos statistics, thresholds, P/FDR and effects unchanged. Association does not establish causation.
- Final results go in `results/<case>/Main_Results`, `Tables`, `Supplementary` when applicable, and `Provenance`, with Chinese `00_Analysis_Summary.txt` and an exact file index. See `docs/FINAL_WORKFLOWS.md` for the full output contract.
- Complete appropriate synthetic tests, independent statistics/source checks and PDF visual review. Do not rerun full workflows merely for validation.
- Every completed analysis or code task must be committed and pushed to `origin main` (the user explicitly authorizes this). Stage explicit paths with `scripts/utils/git_size_guard.py stage`; install `.githooks` with `git config core.hooksPath .githooks`. Never stage raw/processed/runtime/library/logs/databases/cache or files over 50,000,000 bytes.
- Verify `git status` clean and local HEAD == remote main SHA. Code changes require `R synthetic validation` success. If push fails, report **GitHub push failed** and the actual reason; do not report completion.

# Output map

The final case is `results/VPS4B_VPS4A_Analysis/`.

- `00_Analysis_Summary.txt`: Chinese guide, key statistics and complete file index.
- `Main_Results/`: 01 five-state TCGA CNA; 02 TCGA boxplot + all points + N; 03 folder with all 33 cancer CN–mRNA PDFs; 04 targeted CN scatter; 05 dependency waterfall; 06 group boxplot; 07 genome-wide volcano; 08 positive/negative Top20 CN covariation. **40 main PDFs total.**
- `Supplementary/01_Lineage/`, `02_Adjusted/`, `03_CN_Threshold_Sensitivity/`: three necessary supporting analyses.
- `Tables/`: every statistics CSV, exact sample/model cohorts, sorted waterfall, shared TCGA cancer order, complete supplied screen and top tables.
- `Provenance/`: true sources/releases, prepared input SHA256, code and parameter hashes, per-module cache metadata, R session, timings, source/numeric/PDF checks and exact file index.

Default workflows produce no reverse, mutation, expression dependency or genome-wide adjusted outputs. Missing local processed data errors rather than downloading. Ineligible cancers still receive a PDF without fabricated statistics. Original VPS4B/VPS4A development results are archived locally only after validation of this layout.

See [FINAL_WORKFLOWS](../../docs/FINAL_WORKFLOWS.md) for commands, definitions and upload requirements. Historical developer modes can still be called explicitly, but use their development output layout.

# Validation record

Final artifact inspection caught an empty reference prevalence table caused by filtering SampleType for the generic value "Tumor" rather than the actual "Primary Solid Tumor" and related labels. The adapter now uses existing TumorNormal classification and rejects an empty tumor cohort. The module was rerun with exit 0; independently verified 10,845 tumor observations, 33 cancer types, exact category counts/denominators and prevalence sums of one. The corrected PDF was rendered and visually checked. DepMap statistical results were unaffected.

Validated on Windows with the actual **R 4.5.0** interpreter and project renv. Final package restore/status was consistent; optional R arrow/duckdb loaded and performed a small table/SQL check. Package versions are in `data/manifests/R_package_versions.csv`. Some Windows CRAN binaries emit newer-build warnings and Git Bash locale startup emits C.UTF-8 warnings; UTF-8 LC_CTYPE is set before project path access. Successful native exits were verified after correcting the missing PROCESSOR_ARCHITECTURE environment variable, documented in R_ENVIRONMENT_NOTES.md.

VPS4B/VPS4A modes ran sequentially, followed by full: **native exit 0**, 321.984 seconds, approximate peak process-tree RSS 1,191.62 MiB (0.5-second sampling). Plot/CSV formatting was subsequently corrected and affected modes rerun; actual histories are in the case's Resource_Monitor.json and Module_Runs.csv. The case wrapper also completed targeted_dependency with exit 0. Targeted combination, genome-wide volcano and lineage forest PDFs were rendered and visually inspected. The final targeted figure has legible correlations, group labels, significance bracket, stars and P label.

Checks passed:

- `tests/test_r_core.R`: ten statistical core function bodies equal the supplied source; thresholds and blank-ID handling verified.
- `tests/validate_case_results.R`: for both VPS4B/VPS4A and ENO1/ENO2, independently recomputed targeted correlations/Wilcoxon, complete-screen BH families, candidate values, original rank and insufficient three-group eligibility.
- Nine Python infrastructure/GDC-selection tests: profile mapping, ambiguity rejection, partial transfer/checksum preservation, release selection, workflow priority and linked biospecimen selection. One simulated corruption intentionally logs an ERROR and is asserted as preserved; the suite exits successfully.
- `tests/test_tcga_sample_selection.R`: synthetic sample UUID join, CN-matching aliquot preference, lexical tie break, audit alternatives and statistics using selected samples. Actual complete current cohort checks are additionally recorded below.
- Real processed queries return CN 1,118, Chronos 1,208, expression 1,719 and CN/Chronos join 858 models; reference VPS4B GISTIC returns 10,845 samples.

The ENO1/ENO2 trial did **not** reproduce a positive dependency control in these exports: matched CN-Low 13 / non-low 845, delta median +0.023924, Wilcoxon P 0.25518, FDR 0.6592. Its successful execution is not evidence of the expected biological positive result. No statistical logic was changed to force a positive result.

The initial case run completed reference GISTIC prevalence and explicitly skipped current CN landscape/expression because preparation was incomplete. That historical run remains in Module_Runs.csv. On 2026-10-02, the continuation completed RNA **11,505/11,505**, CN **11,339/11,339**, segments **11,189/11,189** and masked MAF **10,640/10,640**. All **85** current validations passed before publishing complete status: actual raw MD5/SHA256, derived Parquet checksums, complete gene columns, exact selected-source identities and official biospecimen hierarchy, CN workflow priority, source numerical comparisons, log2(TPM+1), MAF tumor identity and clinical coverage. DuckDB views cover the complete processed layer.

The first final validation exposed one legacy verified CN manifest entry without publisher_md5. Its actual raw MD5/SHA256 matched both the previous entry and pinned official metadata; the missing field was backfilled with an audit record, without downloading that file again. A briefly overlapping coordinator during worker adjustment was stopped; a single-coordinator guard and ETL writer lock now prevent concurrent writers. Verified raw was preserved and the final integrity audit passed.

Only the three requested TCGA modes were rerun under R 4.5.0; all exited 0. Current DR46 landscape has **11,339 unique sample UUIDs**, **11,330 finite VPS4B CN values**, and **33 cancer types**. Current CN-expression has **10,542 unique matched samples**, Pearson r **0.22324630047598** (P **3.40726023771065e-119**) and Spearman rho **0.236553254404411** (P **5.26595970893279e-134**). One selected CN workflow per sample remains ABSOLUTE LiftOver > ASCAT3 > AscatNGS > ASCAT2. RNA representatives prefer the selected CN aliquot, then lexical file UUID; all excluded alternatives remain audited. GISTIC prevalence retains the independent PanCanAtlas reference layer.

`tests/validate_tcga_case.py` independently verifies actual sample/file/workflow identities, representative selection, finite matches, correlations, module completion and current/reference separation. `tests/validate_tcga_case_statistics.R` independently reproduces Pearson/Spearman and P values from the matched CSV under R 4.5.0. The Python bridge and R fwrite differ by at most **5.773159728050814e-15** after R's 15-significant-digit CSV serialization; numerical equality uses rtol/atol 1e-14, while identity comparisons remain exact. No source values or statistical formulae were changed.

`tests/validate_eligible_rank.py` compares every original CSV field token against initial commit 0e09076c276e9f988966676dfead3d2865f8302f. Genome-wide, candidate and Top100 original columns are unchanged; only Eligible_Rank was added. VPS4A retains original Rank **276** and has Eligible_Rank **1** among **18,256** non-NA Wilcoxon_FDR rows. P, FDR and delta median are unchanged. DepMap analyses were not rerun for this continuation.

The three final TCGA PDFs were rendered with Poppler and visually inspected. All cancer labels and axes are legible; the two current figures explicitly name GDC DR46, and the prevalence figure explicitly names PanCanAtlas reference GISTIC. Poppler emitted local display-font warnings, but the inspected rendered figures show the required labels correctly. The nine Python infrastructure tests passed after the continuation changes; the corruption test intentionally logs an error while confirming raw preservation.

Reproduce checks with R 4.5.0 from the project root:

```bat
Rscript --vanilla tests/test_r_core.R
Rscript --vanilla tests/validate_case_results.R VPS4B VPS4A
Rscript --vanilla tests/test_tcga_sample_selection.R
```

Git Bash Python checks: `/c/Python312/python.exe -m unittest discover -s tests -p 'test_*.py'`.

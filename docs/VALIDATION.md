# Validation record

Final artifact inspection caught an empty reference prevalence table caused by filtering SampleType for the generic value "Tumor" rather than the actual "Primary Solid Tumor" and related labels. The adapter now uses existing TumorNormal classification and rejects an empty tumor cohort. The module was rerun with exit 0; independently verified 10,845 tumor observations, 33 cancer types, exact category counts/denominators and prevalence sums of one. The corrected PDF was rendered and visually checked. DepMap statistical results were unaffected.

Validated on Windows with the actual **R 4.5.0** interpreter and project renv. Final package restore/status was consistent; optional R arrow/duckdb loaded and performed a small table/SQL check. Package versions are in `data/manifests/R_package_versions.csv`. Some Windows CRAN binaries emit newer-build warnings and Git Bash locale startup emits C.UTF-8 warnings; UTF-8 LC_CTYPE is set before project path access. Successful native exits were verified after correcting the missing PROCESSOR_ARCHITECTURE environment variable, documented in R_ENVIRONMENT_NOTES.md.

VPS4B/VPS4A modes ran sequentially, followed by full: **native exit 0**, 321.984 seconds, approximate peak process-tree RSS 1,191.62 MiB (0.5-second sampling). Plot/CSV formatting was subsequently corrected and affected modes rerun; actual histories are in the case's Resource_Monitor.json and Module_Runs.csv. The case wrapper also completed targeted_dependency with exit 0. Targeted combination, genome-wide volcano and lineage forest PDFs were rendered and visually inspected. The final targeted figure has legible correlations, group labels, significance bracket, stars and P label.

Checks passed:

- `tests/test_r_core.R`: ten statistical core function bodies equal the supplied source; thresholds and blank-ID handling verified.
- `tests/validate_case_results.R`: for both VPS4B/VPS4A and ENO1/ENO2, independently recomputed targeted correlations/Wilcoxon, complete-screen BH families, candidate values, original rank and insufficient three-group eligibility.
- Nine Python infrastructure/GDC-selection tests: profile mapping, ambiguity rejection, partial transfer/checksum preservation, release selection, workflow priority and linked biospecimen selection. One simulated corruption intentionally logs an ERROR and is asserted as preserved; the suite exits successfully.
- `tests/test_tcga_sample_selection.R`: synthetic sample UUID join, CN-matching aliquot preference, lexical tie break, audit alternatives and statistics using selected samples. This does not validate the future full current RNA/CN cohort.
- Real processed queries return CN 1,118, Chronos 1,208, expression 1,719 and CN/Chronos join 858 models; reference VPS4B GISTIC returns 10,845 samples.

The ENO1/ENO2 trial did **not** reproduce a positive dependency control in these exports: matched CN-Low 13 / non-low 845, delta median +0.023924, Wilcoxon P 0.25518, FDR 0.6592. Its successful execution is not evidence of the expected biological positive result. No statistical logic was changed to force a positive result.

TCGA reference GISTIC prevalence completed. Current CN landscape and CN-expression were explicitly skipped at case run time because full GDC RNA/CN preparation was incomplete. Raw segments and masked MAF completed with publisher MD5 and local SHA256; current clinical/biospecimen mapping and MAF row identity were validated during ETL. Full current numerical/coverage validation remains pending and will be performed by `complete_gdc_pipeline.py`; readiness is not marked complete before that validation.

Reproduce checks with R 4.5.0 from the project root:

```bat
Rscript --vanilla tests/test_r_core.R
Rscript --vanilla tests/validate_case_results.R VPS4B VPS4A
Rscript --vanilla tests/test_tcga_sample_selection.R
```

Git Bash Python checks: `/c/Python312/python.exe -m unittest discover -s tests -p 'test_*.py'`.

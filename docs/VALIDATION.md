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

## Framework extensions (2026-10-02)

R 4.5.0 was confirmed; `scripts/R/verify_environment.R` restored the existing renv and verified synchronization (exit 0). No TCGA/DepMap download, raw write or re-ETL was performed. Before the requested full rerun, all prior VPS4B/VPS4A outputs were copied to `.runtime/prior_case/d07ad561ed91064182949ca99f185d8e156818b3`; Git history also preserves those outputs.

All synthetic/core checks passed: test_r_core.R, test_cn_threshold_sensitivity.R, test_expression_dependency.R, test_tcga_by_cancer.R, test_genomewide_adjusted.R and test_tcga_sample_selection.R. They check threshold equality/group sizes/insufficient-N skips, type-7 10/20% quantiles/ties/NA exclusion, correlations/Wilcoxon and separate BH families, per-cancer N>=20/Fisher-z CI, raw and standardized adjusted regression, exact shuffled ModelID alignment and aliased CN/lineage rejection. Statuses are recorded in synthetic_test_results.csv. GitHub Actions runs these same six tests with R 4.5.0 and a scoped restore of locked test dependencies; it does not download biological data. Remote run status is available in the repository's Actions tab.

The prior VPS4B/VPS4A validation passed before full. The expanded full run exited 0 in 174.579 seconds (approximate peak RSS 1,416.39 MiB); the explicit adjusted genome-wide run exited 0 in 67.187 seconds (1,456.13 MiB). After adding explicit rank-deficient-model protection and plot-label wrapping, only affected extension modes were rerun. Histories retain all actual executions. ENO1/ENO2 threshold sensitivity and expression dependency both exited 0 in the separately named Method_Test folder; no positive result was forced.

`tests/validate_preserved_results.py` verified all original field tokens in **32 prior CSVs** against completed commit d07ad561ed91064182949ca99f185d8e156818b3. `tests/validate_eligible_rank.py` additionally verifies initial screen tokens, including P/FDR/delta/original Rank. User-facing VPS4A Eligible_Rank remains **1 / 18,256**, with historical Rank **276** retained. New Eligible_N and reordered candidate columns are metadata additions.

`tests/validate_extension_results.R` independently recomputed both cases' threshold/expression tests from existing processed Parquet, all current DR46 per-cancer correlations/BH/Fisher intervals and adjusted regression, both complete adjusted-screen BH families, adjusted rank, and processed-source regressions for VPS4A, ENO2 and RPL5. All passed; tiny positive P values are compared on the log scale, while matching floating-point-underflow zeroes are explicitly distinguished. The original TCGA mapping/overall validation also passed unchanged: matched N 10,542, Pearson 0.22324630047598, Spearman 0.236553254404411.

Actual new VPS4B/VPS4A results: cancer-adjusted CN beta **0.201073386062131**, P **1.16268549540614e-170**; standardized beta **0.231535006800212**. All 33 cancers are eligible; 23 Pearson and 22 Spearman cancer tests have their respective BH-FDR <0.05. Pan-cancer overall correlation may be influenced by between-cancer differences. The lineage-adjusted complete screen has **18,435** eligible continuous tests; VPS4A adjusted rank **1**, Beta_CN **0.875849024824707**, FDR_CN **8.78814561637107e-43**, Beta_CNLow **-0.683392634066967**, FDR_CNLow **2.1560852112041e-49**. This adjusted-model rank is distinct from original Wilcoxon Eligible_Rank.

VPS4B threshold 0.585/0.50/0.40 has N_low 76/39/3, delta median -0.560742/-0.689454/-1.695244 and eligible-threshold BH-FDR 3.46013e-22/1.53670e-12/0.00335167. Threshold 0.35 has zero low models and is not tested. Expression bottom 10%/20% has N_low 114/228, delta -0.165163/-0.171665, FDR 2.17803e-10/4.89806e-18. Its continuous matched N is 1,140; Pearson 0.318281, Spearman 0.318021.

ENO1/ENO2 has no significant signal in these defaults: threshold 0.585 N_low=13, delta +0.0239242, P/FDR=0.255183; threshold 0.50 has only two low models and 0.40/0.35 have none, so those tests are skipped. Expression bottom 10% has delta +0.0121205, P=0.343428/FDR=0.686855; bottom 20% delta +0.000348954, P/FDR=0.952566. Continuous expression correlation is Pearson -0.0285117 (P=0.336147), Spearman -0.0288874 (P=0.329816), N=1,140. These results do not explain the earlier non-replication by a significant CN threshold or expression-low state in the supplied exports.

All six newly generated result PDFs were rendered and visually inspected. The per-cancer forest retains Fisher-z intervals; the adjusted volcano highlights VPS4A; expression plots separate continuous, 10% and 20% panels. Sensitivity captions explicitly identify skipped thresholds and their group sizes. Display-font warnings from local Poppler do not affect the inspected labels.

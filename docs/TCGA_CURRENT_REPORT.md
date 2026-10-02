# GDC current preparation snapshot

Snapshot: **2026-10-02T21:49:41+08:00**. Current status: **complete**. All selected raw files, complete RNA/CN ETL, DuckDB views and identity/numerical validation have completed. The current layer is available; case output status is determined by its latest R run.

1. Official release: **46.0**, API version 1, tag 8.5.0; query snapshot 2026-10-02T17:01:45+08:00.
2. TCGA projects: **33**.
3–4, 6–7. File coverage and bytes (MD5 verified / selected):

| Kind | Files verified / selected | Bytes verified / selected |
|---|---:|---:|
| rna | 11,505 / 11,505 | 48,743,334,626 / 48,743,334,626 |
| cn | 11,339 / 11,339 | 38,946,368,973 / 38,946,368,973 |
| segments | 11,189 / 11,189 | 72,026,450 / 72,026,450 |
| mutation | 10,640 / 10,640 | 820,724,927 / 820,724,927 |

5. Selected CN workflows: {"ABSOLUTE LiftOver": 10677, "AscatNGS": 148, "ASCAT3": 418, "ASCAT2": 96}. Priority per sample UUID: ABSOLUTE LiftOver > ASCAT3 > AscatNGS > ASCAT2. A deterministic same-workflow aliquot/file tie break is audited.
8. Clinical: **11,428 cases**; 175,065 biospecimen rows; 44,673 sample/file map rows.
9. Current processed directory at snapshot: **13,751,276,246 bytes**. Segments: 1,070,293 rows; masked MAF: 2,570,542 rows. Full RNA/CN matrices have been generated and validated.
10. Outstanding raw files: **0**. Status counts: {"verified": 44673}. Complete UUID list is in [gdc_data_manifest.csv](../data/manifests/gdc_data_manifest.csv). Earlier DTT TLS errors, a truncated large bundle and a corrected live-tracking indentation error are retained in local logs; none is represented as a successful transfer. Verified raw was retained.
11. Free D-drive space at snapshot: **130,221,985,792 bytes**. Initial conservative required estimate 170,324,542,336 bytes was within 80% of initial free 237,361,651,712 bytes; ongoing transfers also check remaining space.

Segments for 150 selected CN samples are not published in the queried source. 10,527 segment selections use a documented different-workflow fallback, mostly because matching ABSOLUTE segments are not published. Do not interpret these as exact pipeline-matched segments. No purity/ploidy is invented.

Two file-level MAF biospecimen inconsistencies were found. Each source Tumor_Sample_UUID/barcode was resolved against the full official case hierarchy, with case/barcode validation; source tumor identity is canonical and the file API association is retained for audit. See [identity QC](../data/manifests/gdc_maf_identity_qc.json).

Local continuation: `scripts/utils/complete_gdc_pipeline.py` resumes verified files, then performs complete ETL and numerical/identity validation. Progress: `data/raw/tcga/gdc_current_DR46/manifests/live/pipeline_status.json`; detail: `logs/gdc_pipeline.log`. A failed phase is recorded and can be resumed with the same script. Only successful validation marks current complete. Run `scripts/utils/publish_gdc_snapshot.py` explicitly to refresh repository snapshots, and rerun desired R TCGA modes to update example results. Do not run competing coordinators.

Validation: **passed**; 85 recorded checks. Actual raw MD5/SHA256 audit, all selected RNA/CN identities, full biospecimen mapping, original CN workflow priority, source numerical values and exact log2(TPM+1) are checked; no reference substitution. See [validation report](../data/manifests/gdc_validation_report.json) and [raw integrity](../data/manifests/gdc_raw_integrity.json).

| Processed dataset | Rows | Columns including identity fields | Bytes |
|---|---:|---:|---:|
| TCGA_STAR_UnstrandedCounts | 11,505 | 60,669 | 930,848,580 |
| TCGA_STAR_TPM | 11,505 | 60,669 | 2,682,181,013 |
| TCGA_STAR_FPKM | 11,505 | 60,669 | 2,527,920,671 |
| TCGA_STAR_FPKMUQ | 11,505 | 60,669 | 2,539,531,741 |
| TCGA_STAR_log2TPMplus1 | 11,505 | 60,669 | 3,493,702,851 |
| TCGA_GeneLevel_CN | 11,339 | 60,633 | 232,306,410 |
| TCGA_CN_Segments | 1,070,293 | 26 | 94,409,553 |
| TCGA_Masked_Somatic_Mutation | 2,570,542 | 158 | 1,216,475,251 |
| TCGA_Clinical | 11,428 | 20 | 7,335,564 |
| TCGA_Biospecimen | 175,065 | 17 | 17,334,135 |
| TCGA_Sample_Map | 44,673 | 15 | 3,905,657 |

Parquet is samples/aliquots × original gene IDs, ZSTD. STAR has five matrices: counts, TPM, FPKM, FPKM-UQ and log2(TPM+1), with gene annotation and N_* library summaries stored separately. CN is one selected workflow per sample. DuckDB contains views over these Parquet files, without duplicating matrices. RNA/CN do not silently collapse samples to patients.

VPS4B/VPS4A current modules completed under R 4.5.0; independently checked at 2026-10-02T21:45:49+08:00. CN landscape: 11,339 unique selected samples, 11,330 finite VPS4B CN values. CN-expression: 10,542 unique sample UUID pairs, Pearson r=0.223246300476, Spearman rho=0.236553254404. Source workflow/file selection, RNA representatives and current/reference separation passed [case validation](../data/manifests/tcga_case_validation.json). PanCanAtlas GISTIC prevalence remains a separate reference module.

Pan-cancer overall correlation may be influenced by between-cancer differences. Report overall correlation alongside cancer-adjusted CN beta; these expression associations do not establish dependency or causality.

Cancer-adjusted current DR46 Expression ~ CN + CancerType: CN beta **0.201073386062131**, P **1.16268549540614e-170**, N **10542**. Per-cancer correlations, Fisher-z CI and separate BH families are in the case's TCGA_CN_Expression_ByCancer.csv; independent extension checks are in [extension validation](../data/manifests/extension_results_validation.json).

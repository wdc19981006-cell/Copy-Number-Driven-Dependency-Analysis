# GDC current preparation snapshot

Snapshot: **2026-10-02T19:03:52+08:00**. Current status: **incomplete**. Download/ETL/validation continues locally; this report is not a live counter. The VPS4B/VPS4A example records the readiness at its run time.

1. Official release: **46.0**, API version 1, tag 8.5.0; query snapshot 2026-10-02T17:01:45+08:00.
2. TCGA projects: **33**.
3–4, 6–7. File coverage and bytes (MD5 verified / selected):

| Kind | Files verified / selected | Bytes verified / selected |
|---|---:|---:|
| rna | 11,503 / 11,505 | 48,734,853,759 / 48,743,334,626 |
| cn | 819 / 11,339 | 2,813,136,821 / 38,946,368,973 |
| segments | 11,189 / 11,189 | 72,026,450 / 72,026,450 |
| mutation | 10,640 / 10,640 | 820,724,927 / 820,724,927 |

5. Selected CN workflows: {"ABSOLUTE LiftOver": 10677, "AscatNGS": 148, "ASCAT3": 418, "ASCAT2": 96}. Priority per sample UUID: ABSOLUTE LiftOver > ASCAT3 > AscatNGS > ASCAT2. A deterministic same-workflow aliquot/file tie break is audited.
8. Clinical: **11,428 cases**; 175,065 biospecimen rows; 44,673 sample/file map rows.
9. Current processed directory at snapshot: **1,339,461,003 bytes**. Segments: 1,070,293 rows; masked MAF: 2,570,542 rows. Full RNA/CN matrices remain pending while raw downloads are incomplete.
10. Outstanding raw files: **10,522**. Status counts: {"verified": 34151, "not_downloaded": 10522}. Complete UUID list is in [gdc_data_manifest.csv](../data/manifests/gdc_data_manifest.csv). Earlier DTT TLS errors, a truncated large bundle and a corrected live-tracking indentation error are retained in local logs; none is represented as a successful transfer. Verified raw was retained.
11. Free D-drive space at snapshot: **178,764,800,000 bytes**. Initial conservative required estimate 170,324,542,336 bytes was within 80% of initial free 237,361,651,712 bytes; ongoing transfers also check remaining space.

Segments for 150 selected CN samples are not published in the queried source. 10,527 segment selections use a documented different-workflow fallback, mostly because matching ABSOLUTE segments are not published. Do not interpret these as exact pipeline-matched segments. No purity/ploidy is invented.

Two file-level MAF biospecimen inconsistencies were found. Each source Tumor_Sample_UUID/barcode was resolved against the full official case hierarchy, with case/barcode validation; source tumor identity is canonical and the file API association is retained for audit. See [identity QC](../data/manifests/gdc_maf_identity_qc.json).

Local continuation: `scripts/utils/complete_gdc_pipeline.py` resumes verified files, then performs complete ETL and numerical/identity validation. Progress: `data/raw/tcga/gdc_current_DR46/manifests/live/pipeline_status.json`; detail: `logs/gdc_pipeline.log`. A failed phase is recorded and can be resumed with the same script. Only successful validation marks current complete. Run `scripts/utils/publish_gdc_snapshot.py` explicitly to refresh repository snapshots, and rerun desired R TCGA modes to update example results. Do not run competing coordinators.

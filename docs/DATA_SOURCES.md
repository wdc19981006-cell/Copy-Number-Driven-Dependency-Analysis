# Current data sources

As of 2026-10-02T19:03:52+08:00, all **13 user-supplied DepMap 26Q1 exports are active** and organized under `data/raw/depmap/26Q1/`. Before/after SHA256 checks verified unchanged content. No duplicate exports were present. Original names, dimensions, imported timestamps, raw paths and SHA256 are in [DepMap manifest](../data/manifests/depmap_26Q1_manifest.csv). Original download dates and publisher MD5 for these exports were not supplied; local SHA256 is a content fingerprint, not publisher verification.

| Active canonical file | Rows × columns including IDs | Unique models | Gene columns |
|---|---:|---:|---:|
| CopyNumber_WGS_26Q1.csv | 1118 × 18614 | 1118 | 18613 |
| CRISPRScreenMap.csv | 1451 × 2 | 1208 | 0 |
| CRISPR_Chronos_26Q1.csv | 1208 × 18532 | 1208 | 18531 |
| CRISPR_GeneDependency_26Q1.csv | 1208 × 18532 | 1208 | 18531 |
| Mutation_Damaging_26Q1.csv | 1968 × 19506 | 1968 | 19505 |
| Expression_26Q1.csv | 1719 × 19216 | 1719 | 19215 |
| Mutation_Hotspot_26Q1.csv | 1968 × 554 | 1968 | 553 |
| MolecularSubtypes_26Q1.csv | 2044 × 40 | 2044 | 0 |
| Model.csv | 2154 × 49 | 2154 | 0 |
| ModelCondition.csv | 2826 × 16 | 2044 | 0 |
| OmicsProfiles.csv | 4819 × 18 | 2044 | 0 |
| OmicsSignatures_26Q1.csv | 1968 × 7 | 1968 | 0 |
| ScreenSequenceMap.csv | 3346 × 14 | 1257 | 0 |

Exact ModelID overlaps: CN/Chronos **858**, CN/expression **1,105**, expression/Chronos **1,140**, all three **852**. Model metadata contains 2,154 unique models. Metadata/profile tables can have multiple rows per model; matrix ModelIDs are unique. CN is WGS relative CN; Chronos is gene effect, distinct from dependency probability. Expression keeps the original export scale. Gene symbol columns and blank first ID headers are adapted in memory, without changing raw.

The [official DepMap catalog](https://depmap.org/portal/api/no-captcha/download/files) selected Public 26Q1 in the acquisition snapshot. Initial official-file acquisition was skipped at the user's request when links required browser verification; these later local exports now satisfy the R inputs. Several supplied names contain `subsetted`; their completeness against official full files is unverified. Full somatic variant, fusion supplement and all-gene stranded expression were not supplied. Optional combined/array CN was unavailable in that catalog. No older release is substituted.

TCGA current and reference: see [TCGA sources](TCGA_DATA_SOURCES.md) and [dated current preparation report](TCGA_CURRENT_REPORT.md). Current RNA/CN remain incomplete at this snapshot; reference GISTIC is available and explicitly labeled.

Manifest scopes: `depmap_26Q1_manifest.csv` and `depmap_26Q1_local_qc.json` describe the current local exports; `gdc_*` are dated current snapshots; `data_manifest.csv`, `not_downloaded.csv`, `preparation_summary.json`, `depmap_26Q1_qc.json` and the initial preparation report describe the earlier official-download stage, not the current local-export state. [current_missing_data.csv](../data/manifests/current_missing_data.csv) provides the current missing-data summary. Ongoing GDC state lives in the ignored `raw/manifests/live` directory.

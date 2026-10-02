# VPS4B → VPS4A analysis

R 4.5.0 ; DepMap 26Q1 user-supplied exports.

CN groups use analysis-defined thresholds on log2(relative CN + 1); these are not official GISTIC states.
Continuous dependency and adjusted CN analyses use CN_log exactly as in the supplied R code; CN-expression and covariation use CN_relative.
BH adjustments are separate for each complete screen/family; rank orders Wilcoxon FDR then delta median, exactly as supplied.
Export completeness relative to official full-release files is unverified. All supplied Chronos columns are tested.

Deep CN Loss N = 1
Shallow CN Loss N = 87
CN Non-Low N = 1030
Raw CN-Low N = 88 ; CN-NonLow N = 1030
Matched CN-NonLow N = 782
Matched CN-Low N = 76

VPS4A rank = 276 ; delta median = -0.56074 ; Wilcoxon P = 1.1534e-22 ; FDR = 2.1056e-18
Eligible_Rank = 1 (non-NA Wilcoxon_FDR genes in the original sorted order; original Rank/P/FDR/effects unchanged).
The supplied setorder() places NA FDR first: 275 untested rows precede this candidate. Original Rank is retained unchanged. The candidate has the minimum eligible-test FDR: TRUE

CN_vs_expression Pearson 0.45427 P 2.3042e-57 ; Spearman 0.55246 P 2.6288e-89 ; N 1105

expression_vs_dependency Pearson 0.31828 P 2.9762e-28 ; Spearman 0.31802 P 3.3079e-28 ; N 1140

Targeted Pearson 0.48306 P 2.2755e-51 ; Spearman 0.45705 P 1.6495e-45 ; Wilcoxon P 1.1534e-22 ; N 858

Adjusted CN_log beta 0.87585 SE 0.057013 t 15.362 P 4.7671e-47

Adjusted I(CN_binary == "CN-Low")TRUE beta -0.68339 SE 0.041158 t -16.604 P 1.1696e-53

Eligible lineages: 10 ; FDR < 0.05: 8
Forest displays delta median with percentile 95% CI from 1,000 bootstrap draws (seed 1234).
Head and Neck : delta -0.9993 CI [ -1.928 , -0.06649 ] FDR 0.01092
Soft Tissue : delta -0.9855 CI [ -1.7 , -0.2278 ] FDR 0.01226
Breast : delta -0.8453 CI [ -1.792 , -0.08347 ] FDR 0.01263
Lung : delta -0.6753 CI [ -1.051 , -0.3487 ] FDR 3.43e-05
Pancreas : delta -0.6102 CI [ -1.182 , -0.07251 ] FDR 0.0318
Esophagus/Stomach : delta -0.5079 CI [ -0.6896 , -0.2699 ] FDR 0.001214
Uterus : delta -0.3916 CI [ -1.094 , -0.1878 ] FDR 0.01537
Lymphoid : delta -0.3383 CI [ -0.5655 , -0.1975 ] FDR 0.01226
Bowel : delta -0.1678 CI [ -0.801 , 0.2044 ] FDR 0.2906
Ovary/Fallopian Tube : delta -0.1255 CI [ -1.178 , -0.02588 ] FDR 0.1024

forward VPS4B → VPS4A delta median -0.56074 Pearson 0.48306 Wilcoxon P 1.1534e-22

reverse VPS4A → VPS4B delta median 0.0040105 Pearson 0.24466 Wilcoxon P 0.94599

CN covariation is an association/co-deletion signal; it does not establish physical proximity or causality.

GDC current DR46 CN landscape: 11339 unique selected samples; finite gene CN 11330 ; cancer types 33 . One selected CN workflow per sample; no reference fallback.
GDC current DR46 CN-expression: matched N 10542 ; Pearson 0.223246 P 3.4073e-119 ; Spearman 0.236553 P 5.266e-134 . Exact sample UUID join, one representative RNA aliquot per sample, audit retained; RNA is log2(TPM+1).

Module status (latest per mode):
- genomewide_dependency: completed.
- lineage_dependency: completed.
- reverse_dependency: completed.
- cn_covariation: completed.
- qc: completed.
- depmap_cn_expression: completed.
- mutation_dependency: completed.
- adjusted_dependency: completed.
- targeted_dependency: completed.
- tcga_cn_landscape: completed.
- tcga_cn_expression: completed.
- tcga_cna_prevalence: completed.

These are observational cell-line associations. Lineage adjustment reduces measured lineage confounding; it does not establish a causal synthetic-lethal mechanism or clinical benefit.

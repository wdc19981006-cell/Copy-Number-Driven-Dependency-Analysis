# ENO1 → ENO2 analysis

R 4.5.0 ; DepMap 26Q1 user-supplied exports.

CN groups use analysis-defined thresholds on log2(relative CN + 1); these are not official GISTIC states.
Continuous dependency and adjusted CN analyses use CN_log exactly as in the supplied R code; CN-expression and covariation use CN_relative.
BH adjustments are separate for each complete screen/family; rank orders Wilcoxon FDR then delta median, exactly as supplied.
Export completeness relative to official full-release files is unverified. All supplied Chronos columns are tested.


Module status (latest per mode):
- expression_dependency: completed.
- cn_threshold_sensitivity: completed.

These are observational cell-line associations. Lineage adjustment reduces measured lineage confounding; it does not establish a causal synthetic-lethal mechanism or clinical benefit.

CN threshold sensitivity (analysis-defined; not official DepMap GISTIC):
- Threshold 0.585 : N_low/nonlow 13 / 845 ; delta median 2.39242e-02 ; P 2.55183e-01 ; FDR 2.55183e-01 ; eligible TRUE
- Threshold 0.5 : N_low/nonlow 2 / 856 ; delta median 2.16386e-02 ; P NA ; FDR NA ; eligible FALSE
- Threshold 0.4 : N_low/nonlow 0 / 858 ; delta median NA ; P NA ; FDR NA ; eligible FALSE
- Threshold 0.35 : N_low/nonlow 0 / 858 ; delta median NA ; P NA ; FDR NA ; eligible FALSE

Expression-defined dependency (independent of CN-loss and mutation):
Continuous expression -> Chronos: N 1140 ; Pearson -2.85117e-02 P 3.36147e-01 ; Spearman -2.88874e-02 P 3.29816e-01
- Bottom 10 %: cutoff 9.31387e+00 ; N_low/nonlow 114 / 1026 ; delta median 1.21205e-02 ; P 3.43428e-01 ; FDR 6.86855e-01 ; eligible TRUE
- Bottom 20 %: cutoff 9.71266e+00 ; N_low/nonlow 228 / 912 ; delta median 3.48954e-04 ; P 9.52566e-01 ; FDR 9.52566e-01 ; eligible TRUE

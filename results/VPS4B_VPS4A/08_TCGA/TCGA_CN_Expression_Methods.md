GDC current DR46. Expression = log2(TPM + 1); CN = source total gene-level CN.
Pan-cancer overall correlation may be influenced by between-cancer differences.
Report overall correlations alongside cancer-adjusted CN beta; association does not establish causality.
Per-cancer correlations require N >= 20 and nonconstant CN/expression. BH families are separate for Pearson and Spearman.
Fisher-z CI = tanh(atanh(r) +/- qnorm(0.975)/sqrt(N-3)).
Adjusted regressions require known cancer type; standardization uses global sample SD on that same complete regression cohort.

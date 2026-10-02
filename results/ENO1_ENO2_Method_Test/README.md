# ENO1 → ENO2 method test

R 4.5.0; existing user-supplied DepMap 26Q1 processed Parquet. This is a method test without a required positive result. CN threshold sensitivity and expression-defined dependency remain independent; no deficient group combines CN, expression and mutation.

CN_log=log2(relative CN+1); all thresholds are analysis-defined, not official DepMap GISTIC classification. At 0.585, low/non-low N=13/845, delta median +0.0239242208, Wilcoxon P/FDR=0.255182944323942. At 0.50, low N=2; at 0.40/0.35, low N=0. Those three tests are ineligible and have NA P/FDR. This sensitivity FDR family is separate from the historical genome-wide Wilcoxon family.

Expression quantiles use finite matched pairs, R type 7, with expression <= cutoff defined as low. Bottom 10%: N=114/1026, delta +0.01212050595, P=0.343427575025889, FDR=0.686855150051779. Bottom 20%: N=228/912, delta +0.00034895405, P/FDR=0.952565685100337. Continuous N=1140, Pearson=-0.0285116980120274, Spearman=-0.0288874249465348. No significant dependency association appears in these defaults; cutoff tuning was not performed.

Actual tables, PDFs and summary are under 10_CN_Threshold_Sensitivity/, 11_Expression_Dependency/ and Summary/. Independent processed-source verification is recorded in data/manifests/extension_results_validation.json. Full TCGA and other DepMap modules were not run in this method-test folder.

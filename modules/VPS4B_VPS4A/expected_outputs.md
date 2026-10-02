# Output map

Results are under results/VPS4B_VPS4A/. Each core module has its named subdirectory:

- 00_QC/: CN_Distribution.pdf, CN_Group_Counts.csv, QC_Statistics.csv, model alignment.
- 01_CN_Expression/: VPS4B_CN_vs_Expression.pdf, CN_Expression_Statistics.csv.
- 02_GenomeWide_Dependency/: complete GenomeWide_Dependency.csv, volcano PDF, Top100, VPS4A_Candidate_Rank.csv.
- 03_Targeted_Dependency/: VPS4B_VPS4A_Statistics.csv, Group_Summary.csv, CellLines.csv and combined PDF; three-group plot only if every group >=3.
- 04_Lineage_Dependency/: VPS4B_VPS4A_Lineage_Dependency.csv and Lineage_Forest.pdf.
- 05_Adjusted_Dependency/: Continuous_CN_Adjusted.csv, CNLow_Adjusted.csv.
- 06_Reverse_Dependency/: Reverse_Dependency_Summary.csv and nested reverse-targeted outputs.
- 07_CN_Covariation/: VPS4B_CN_Covariation.csv and Top_CN_Covariation.pdf.
- 08_TCGA/: current CN landscape/expression if prepared; separately labeled reference GISTIC prevalence.
- 09_Mutation_Dependency/: eligible damaging/hotspot comparisons, or explicit skip reasons.
- Summary/: VPS4B_VPS4A_Summary.md, Key_Statistics.csv, Analysis_Metadata.json, Module_Runs.csv, Resource_Monitor.json and SessionInfo.txt.

Missing data or inadequate group sizes produce documented skips; they never produce fabricated values. Resource_Monitor records observed native process exits and approximate RSS, sampled every 0.5 s. No file above 50 MB is uploaded.

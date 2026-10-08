# Exercise both public output contracts with in-memory statistics and tiny
# artifacts. Structural migration validation must never read a data matrix.
source("scripts/R/validate_final_workflow.R",encoding="UTF-8")
for(workflow in c("geneA_screen","geneA_geneB")) {
 WORKFLOW<-workflow;GENE_B_PROVIDED<-workflow=="geneA_geneB"
 PATHS<-workflow_paths(WORKFLOW,GENE_A,GENE_B)
 RESULT_ROOT<-gsub("\\\\","/",tempfile("wf_layout_",tmpdir=".runtime"))
 wf_create_directories()
 stopifnot(!dir.exists(file.path(RESULT_ROOT,"Tables")),
  dirname(wf_table("TCGA_CNA_Percentage.csv"))==dirname(file.path(RESULT_ROOT,PATHS["cna"])),
  dirname(wf_table("TCGA_Cancer_Order.csv"))==dirname(file.path(RESULT_ROOT,PATHS["landscape"])),
  dirname(wf_table("TCGA_CN_mRNA_CNA_Counts.csv"))==file.path(RESULT_ROOT,PATHS["rna"]),
  dirname(wf_table("Top_Dependency_Candidates.csv"))==dirname(file.path(RESULT_ROOT,PATHS["volcano"])),
  dirname(wf_table("Top_CN_Covariation.csv"))==dirname(file.path(RESULT_ROOT,PATHS["covariation"])),
  grepl("Provenance/Data_Audit/",wf_table("TCGA_Current_Samples.csv"),fixed=TRUE),
  inherits(try(wf_table("Unassigned.csv"),silent=TRUE),"try-error"))
 expected<-unname(unlist(lapply(workflow_modules(WORKFLOW),wf_expected)))
 stopifnot(!anyDuplicated(expected),!any(grepl("Tables/",expected,fixed=TRUE)))
 for(path in expected) {
  full<-file.path(RESULT_ROOT,path)
  switch(tools::file_ext(path),pdf=writeLines("%PDF-1.4 synthetic path fixture",full),
   csv=data.table::fwrite(data.table(Column=1),full),json=jsonlite::write_json(list(fixture=TRUE),full),
   txt=writeLines("fixture",full))
 }
 wf_write(current_prevalence(rna,all_order),"TCGA_CNA_Percentage.csv")
 wf_write(all_stats,paste0("TCGA_",GENE_A,"_CN_mRNA_AllCancer_Statistics.csv"))
 wf_write(empty_screen_statistics(),"GenomeWide_Dependency.csv")
 wf_write(dependency_top(empty_screen_statistics()),"Top_Dependency_Candidates.csv")
 status<-list(Status="SKIPPED",N_low=1L,N_nonlow=11L,min_group_n=3L,CN_log_threshold=.585)
 jsonlite::write_json(status,wf_provenance("GenomeWide_Dependency_Status.json"),auto_unbox=TRUE)
 wf_write_screen_status(status)
 if(GENE_B_PROVIDED) {
  wf_write(s,paste0(GENE_A,"_",GENE_B,"_Targeted_Statistics.csv"))
  wf_write(wf_group_statistics(s),paste0(GENE_A,"_",GENE_B,"_CNlow_vs_Normal_Statistics.csv"))
  wf_write(data.table(Lineage=character(),FDR=numeric(),Delta_median=numeric()),"Lineage_Dependency.csv")
  stopifnot(grepl("04_DepMap_CN_vs_Dependency/",wf_table(paste0(GENE_A,"_",GENE_B,"_CellLines.csv")),fixed=TRUE),
   grepl("05_DepMap_Dependency_Waterfall/",wf_table(paste0(GENE_A,"_",GENE_B,"_Waterfall_Order.csv")),fixed=TRUE),
   grepl("06_DepMap_Dependency_CNlow_vs_Normal/",wf_table(paste0(GENE_A,"_",GENE_B,"_CNlow_vs_Normal_Statistics.csv")),fixed=TRUE),
   grepl("Supplementary/01_Lineage/",wf_table("Lineage_Eligibility.csv"),fixed=TRUE),
   grepl("Supplementary/02_Adjusted/",wf_table("Continuous_CN_Adjusted.csv"),fixed=TRUE),
   grepl("Supplementary/03_CN_Threshold_Sensitivity/",wf_table("CN_Threshold_Sensitivity.csv"),fixed=TRUE))
 }
 wf_summary()
 read_gene<-function(...)stop("Structural validation must not read data")
 stopifnot(validate_result_layout())
 # Reject missing cancers, extra main files, stale indexes and old Tables.
 missing<-file.path(RESULT_ROOT,PATHS["rna"],paste0("COAD_",GENE_A,"_CN_mRNA.pdf"))
 bytes<-readBin(missing,"raw",n=file.info(missing)$size);unlink(missing)
 stopifnot(inherits(try(validate_result_layout(),silent=TRUE),"try-error"))
 writeBin(bytes,missing)
 extra<-file.path(RESULT_ROOT,PATHS["rna"],"Unplanned.pdf");writeLines("%PDF-1.4 extra",extra)
 stopifnot(inherits(try(validate_result_layout(),silent=TRUE),"try-error"));unlink(extra)
 index<-data.table::fread(wf_provenance("File_Index.csv"))
 data.table::fwrite(index[-1],wf_provenance("File_Index.csv"))
 stopifnot(inherits(try(validate_result_layout(),silent=TRUE),"try-error"))
 data.table::fwrite(index,wf_provenance("File_Index.csv"))
 dir.create(file.path(RESULT_ROOT,"Tables"))
 stopifnot(inherits(try(validate_result_layout(),silent=TRUE),"try-error"))
 stopifnot(dirname(normalizePath(file.path(RESULT_ROOT,"Tables"),winslash="/"))==normalizePath(RESULT_ROOT,winslash="/"),
  startsWith(normalizePath(RESULT_ROOT,winslash="/"),paste0(normalizePath(".runtime",winslash="/"),"/")))
 unlink(file.path(RESULT_ROOT,"Tables"),recursive=TRUE)
 stopifnot(validate_result_layout())
 unlink(RESULT_ROOT,recursive=TRUE)
}
cat("PASS: both module-folder contracts, table ownership, extracted group statistics, exact 33 cancers, summaries and negative layout checks.\n")

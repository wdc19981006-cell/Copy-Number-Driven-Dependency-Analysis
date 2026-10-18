invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
source(".Rprofile")
library(data.table);library(dplyr);library(ggplot2)
PROJECT_ROOT<-normalizePath(".",winslash="/")
for(s in c("data_access.R","plotting.R","dependency_screen.R","lineage_analysis.R","adapters.R",
 "additional_modes.R","tcga_analysis.R","extensions_statistics.R","extension_modules.R","extension_summary.R",
 "local_data.R","workflow_statistics.R","workflow_plots.R","final_workflow.R"))source(file.path("scripts/R",s),encoding="UTF-8")
GENE_A<-"FIXTURE_A";GENE_B<-"";GENE_B_PROVIDED<-FALSE;WORKFLOW<-"geneA_screen"
MIN_N<-3L;BOOT_R<-25L;opt<-list(force=FALSE);DATA_MODE<-"local"
TCGA_CANCERS<-sort(unique(unlist(jsonlite::fromJSON("config/tcga_cancer_types.json")$mapping)))
RESULT_ROOT<-tempfile("geneA_project_",tmpdir=".runtime")
dir.create(RESULT_ROOT,recursive=TRUE)
input<-tempfile("fixture_input_",tmpdir=".runtime",fileext=".parquet");writeLines("synthetic release 1",input)
wf_input_paths<-function(module)input
calls<-setNames(integer(7),workflow_modules("geneA_geneB"))
pair<-as.data.table(prepare_cn(data.frame(ModelID=paste0("m",1:30),
 CN_relative=c(rep(0,10),seq(.00001,.009,length.out=10),seq(.6,2,length.out=10)))))
pair[,Chronos:=seq(-2,.2,length.out=.N)+rep(c(0,.1),15)]
original_pair<-copy(pair)
fixture_module<-function(module) {
 calls[module]<<-calls[module]+1L
 for(p in unname(wf_expected(module))) {
  path<-file.path(RESULT_ROOT,p);dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  switch(tools::file_ext(p),pdf=writeLines("%PDF-1.4 synthetic fixture",path),
   csv=data.table::fwrite(data.table(Fixture=1),path),json=jsonlite::write_json(list(fixture=TRUE),path),txt=writeLines("fixture",path))
 }
 if(module=="tcga") {
  wf_write(data.table(CancerType="ACC",CNA=CNA_STATES[1],Percentage=100),"TCGA_CNA_Percentage.csv")
  wf_write(data.table(CancerType=TCGA_CANCERS,Pearson_FDR=1),paste0("TCGA_",GENE_A,"_CN_mRNA_AllCancer_Statistics.csv"))
 }
 if(module=="genomewide_dependency") {
  status<-list(Status="COMPLETED",N_low=20L,N_nonlow=10L,min_group_n=MIN_N,CN_log_threshold=.585)
  jsonlite::write_json(status,wf_provenance("GenomeWide_Dependency_Status.json"),auto_unbox=TRUE);wf_write_screen_status(status)
 }
 if(module=="targeted_dependency") {
  s<-targeted_statistics(pair)
  wf_write(pair,paste0(GENE_A,"_",GENE_B,"_CellLines.csv"))
  wf_write(s,paste0(GENE_A,"_",GENE_B,"_Targeted_Statistics.csv"))
  wf_write(wf_group_statistics(s),paste0(GENE_A,"_",GENE_B,"_CNlow_vs_Normal_Statistics.csv"))
 }
 invisible(NULL)
}
wf_run_tcga<-function()fixture_module("tcga")
wf_run_screen<-function()fixture_module("genomewide_dependency")
wf_run_covariation<-function()fixture_module("cn_covariation")
saved_targeted_run<-wf_run_targeted
wf_run_targeted<-function()fixture_module("targeted_dependency")
wf_run_supplement<-function(module)fixture_module(module)
# Execute the real orchestration/cache/summary/index paths. Statistical source
# validation is tested separately; this fixture has no TCGA or DepMap matrices.
source<-function(file,...) {
 if(grepl("validate_final_workflow.R$",file)) {
  base::source(file,...)
  validate_final_workflow<<-function(){validate_result_layout();invisible(TRUE)}
 }else base::source(file,...)
}
run_fixture<-function(b="") {
 GENE_B<<-b;GENE_B_PROVIDED<<-nzchar(b);WORKFLOW<<-if(nzchar(b))"geneA_geneB" else "geneA_screen"
 run_final_workflow()
}
run_fixture()
base_outputs<-unname(unlist(lapply(wf_base_modules(),wf_expected)))
base_hashes<-wf_hash_outputs(base_outputs)
base_summary<-readLines(wf_provenance("GeneA_Summary.txt"),encoding="UTF-8")
base_keys<-lapply(wf_base_modules(),function(m)jsonlite::fromJSON(wf_provenance(paste0("Cache_",m,".json")))$key)
run_fixture("B1")
first_paths<-unname(unlist(lapply(wf_geneB_modules(),wf_expected)))
first_hashes<-wf_hash_outputs(first_paths)
run_fixture("B2")
stopifnot(all(calls[wf_base_modules()]==1L),all(calls[wf_geneB_modules()]==2L),
 identical(base_hashes,wf_hash_outputs(base_outputs)),identical(first_hashes,wf_hash_outputs(first_paths)),
 identical(wf_geneBs(),c("B1","B2")),validate_result_layout(),
 identical(base_summary,readLines(wf_provenance("GeneA_Summary.txt"),encoding="UTF-8")))
summary<-readLines(file.path(RESULT_ROOT,"00_Analysis_Summary.txt"),encoding="UTF-8")
stopifnot(any(grepl("06. FIXTURE_A → B1 Dependency",summary,fixed=TRUE)),
 any(grepl("06. FIXTURE_A → B2 Dependency",summary,fixed=TRUE)))
run_fixture();run_fixture("B1")
stopifnot(all(calls[wf_base_modules()]==1L),all(calls[wf_geneB_modules()]==2L),validate_result_layout(),
 identical(base_keys,lapply(wf_base_modules(),function(m)jsonlite::fromJSON(wf_provenance(paste0("Cache_",m,".json")))$key)))
# A parameter/code/file change must invalidate only the affected modules.
cache_path<-wf_module_provenance("targeted_dependency","Cache_targeted_dependency.json")
old<-jsonlite::fromJSON(cache_path,simplifyVector=FALSE)
key<-old$key;paths<-unname(wf_expected("targeted_dependency"))
writeLines("tampered",file.path(RESULT_ROOT,paths[1]))
stopifnot(!wf_cache_valid(old,key,paths))
run_fixture("B1")
stopifnot(calls["targeted_dependency"]==3L,all(calls[wf_base_modules()]==1L))
BOOT_R<-BOOT_R+1L;run_fixture("B1")
stopifnot(calls["lineage_dependency"]==3L,all(calls[wf_base_modules()]==1L),calls["adjusted_dependency"]==2L)
# Exercise the production plot-only branch with source-data access forbidden.
load_target_pair<-function(...)stop("Plot-only refresh must not access data")
old<-jsonlite::fromJSON(cache_path,simplifyVector=FALSE)
old$statistics_key<-wf_targeted_statistics_key()
jsonlite::write_json(old,cache_path,pretty=TRUE,auto_unbox=TRUE)
preserved<-setdiff(unname(wf_expected("targeted_dependency")),unname(PATHS["scatter"]))
preserved_hashes<-wf_hash_outputs(preserved)
base_code<-lapply(wf_base_modules(),function(m)wf_code_signature(m)$hash)
target_code<-wf_code_signature("targeted_dependency")$hash
statistics_key<-wf_targeted_statistics_key()
original_scatter<-targeted_scatter_plot
targeted_scatter_plot<-function(pair,stat)original_scatter(pair,stat)
stopifnot(!identical(target_code,wf_code_signature("targeted_dependency")$hash),
 identical(statistics_key,wf_targeted_statistics_key()),
 identical(base_code,lapply(wf_base_modules(),function(m)wf_code_signature(m)$hash)),
 identical(saved_targeted_run(),"plot_regenerated"),identical(preserved_hashes,wf_hash_outputs(preserved)))
targeted_scatter_plot<-original_scatter
# The true scatter, including the near-zero inset, never changes coordinates,
# CN_log, classification, correlations, group statistics or the RNG state.
s<-targeted_statistics(pair);set.seed(814);seed<-.Random.seed
plot<-targeted_scatter_plot(pair,s);built<-ggplot_build(if(inherits(plot,"patchwork"))plot[[1]] else plot)
stopifnot(identical(pair,original_pair),identical(.Random.seed,seed),
 identical(built$data[[1]]$x,pair$CN_relative),identical(built$data[[1]]$y,pair$Chronos),
 identical(targeted_statistics(pair),s),identical(built$data[[2]]$yintercept,0),
 identical(built$data[[3]]$xintercept,2^.585-1),identical(built$data[[4]]$xintercept,2^.35-1),
 all(built$data[[1]]$x>=0),isTRUE(all.equal(s$Relative_CN_Pearson_r,unname(cor(pair$CN_relative,pair$Chronos)))),
 isTRUE(all.equal(s$Relative_CN_Spearman_rho,unname(cor(pair$CN_relative,pair$Chronos,method="spearman")))))
pdf<-file.path(RESULT_ROOT,"scatter_fixture.pdf");workflow_save(plot,pdf)
stopifnot(wf_validate_artifact(pdf))
# Verify entry point uses only Gene A and rejects workflow folder overrides.
entry<-readLines("scripts/R/run_analysis.R")
stopifnot(any(grepl('paste0(GENE_A,"_Analysis")',entry,fixed=TRUE)),
 !any(grepl('paste0(GENE_A,if(GENE_B_PROVIDED)',entry,fixed=TRUE)))
stopifnot(startsWith(normalizePath(RESULT_ROOT,winslash="/"),paste0(normalizePath(".runtime",winslash="/"),"/")))
unlink(RESULT_ROOT,recursive=TRUE)
# Starting with A+B must create complete A provenance and compute A only once.
RESULT_ROOT<-tempfile("geneA_first_pair_",tmpdir=".runtime");calls[]<-0L
run_fixture("B1")
stopifnot(all(calls==1L),file.exists(wf_provenance("Run_Metadata.json")),
 file.exists(wf_provenance("Input_Versions.json")),file.exists(wf_provenance("Module_Runs.csv")),
 is.null(jsonlite::fromJSON(wf_provenance("Run_Metadata.json"))$geneB),
 identical(jsonlite::fromJSON(wf_provenance("Run_Metadata.json"))$workflow,"geneA_screen"))
run_fixture("B2");stopifnot(all(calls[wf_base_modules()]==1L),validate_result_layout())
stopifnot(startsWith(normalizePath(RESULT_ROOT,winslash="/"),paste0(normalizePath(".runtime",winslash="/"),"/")))
unlink(RESULT_ROOT,recursive=TRUE)
unlink(input)
cat("PASS: real orchestration reuses A across screen/B1/B2; appends preserve prior results/summaries; corruption invalidates one module; scatter coordinates/RNG/statistics/thresholds unchanged.\n")

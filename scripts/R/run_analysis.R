#!/usr/bin/env Rscript
Sys.setlocale("LC_CTYPE", "English_United States.utf8")
if(.Platform$OS.type=="windows"&&!nzchar(Sys.getenv("PROCESSOR_ARCHITECTURE")))
 Sys.setenv(PROCESSOR_ARCHITECTURE=switch(R.version$arch,x86_64="AMD64",aarch64="ARM64",stop("Unknown runtime architecture")))
if (as.character(getRversion()) != "4.5.0") stop("Use R 4.5.0 exactly; current: ",getRversion())
root_arg <- commandArgs(trailingOnly=TRUE)
project_i <- match("--project",root_arg)
PROJECT_ROOT <- normalizePath(if(!is.na(project_i))root_arg[project_i+1L] else ".",winslash="/",mustWork=TRUE)
setwd(PROJECT_ROOT)
Sys.setenv(RENV_PATHS_CACHE=file.path(PROJECT_ROOT,".runtime/renv-cache"),RENV_PATHS_ROOT=file.path(PROJECT_ROOT,".runtime/renv-root"),
 RENV_PATHS_LIBRARY=file.path(PROJECT_ROOT,"renv/library"),RENV_PATHS_SANDBOX=file.path(PROJECT_ROOT,".runtime/renv-sandbox"))
source(file.path(PROJECT_ROOT,"renv/activate.R"))
suppressPackageStartupMessages({
 library(data.table);library(dplyr);library(ggplot2);library(ggpubr);library(ggrepel)
 library(patchwork);library(broom);library(scales);library(stringr);library(optparse)
})
options(stringsAsFactors=FALSE,scipen=999)
data.table::setDTthreads(4)
option_list <- list(make_option("--mode",default=NULL),make_option("--workflow",default=NULL),make_option("--geneA",default=NULL),
 make_option("--geneB",default=NULL),make_option("--project",default="."),make_option("--data_mode",default="local"),
 make_option("--force",action="store_true",default=FALSE),
 make_option("--bootstrap",type="integer",default=1000),make_option("--min_group_n",type="integer",default=3),make_option("--output_case",default=NULL))
opt <- parse_args(OptionParser(option_list=option_list))
WORKFLOW <- opt$workflow
if(!is.null(WORKFLOW)&&!is.null(opt$mode))stop("Choose --workflow or --mode")
if(!identical(opt$data_mode,"local"))stop("Ordinary analysis requires DATA_MODE=local. Use separate maintenance commands.")
if(!is.null(WORKFLOW)&&!WORKFLOW %in% c("geneA_screen","geneA_geneB"))stop("Unsupported workflow")
if(!is.null(WORKFLOW)&&is.null(opt$geneA))stop("--geneA is required for a workflow")
GENE_B_PROVIDED <- !is.null(opt$geneB)
if(identical(WORKFLOW,"geneA_geneB")&&!GENE_B_PROVIDED)stop("geneA_geneB requires --geneB")
if(identical(WORKFLOW,"geneA_screen")&&GENE_B_PROVIDED)stop("Use geneA_geneB when providing --geneB")
MODE <- if(is.null(opt$mode))"targeted_dependency" else opt$mode
GENE_A <- toupper(if(is.null(opt$geneA))"ENO1" else opt$geneA)
GENE_B <- toupper(if(GENE_B_PROVIDED)opt$geneB else if(is.null(WORKFLOW))"ENO2" else "")
if(MODE=="genomewide_adjusted_dependency"&&!GENE_B_PROVIDED)GENE_B<-"ALL"
BOOT_R <- opt$bootstrap;MIN_N <- opt$min_group_n
stopifnot(BOOT_R>0,MIN_N>=3,grepl("^[A-Z0-9_.-]+$",GENE_A),!nzchar(GENE_B)||grepl("^[A-Z0-9_.-]+$",GENE_B))
RESULT_ROOT <- file.path(PROJECT_ROOT,"results",paste0(GENE_A,"_",GENE_B))
if(!is.null(WORKFLOW))RESULT_ROOT<-file.path(PROJECT_ROOT,"results",paste0(GENE_A,if(GENE_B_PROVIDED)paste0("_",GENE_B),"_Analysis"))
if(MODE=="genomewide_adjusted_dependency"&&!GENE_B_PROVIDED)RESULT_ROOT<-file.path(PROJECT_ROOT,"results",paste0(GENE_A,"_GenomeWide"))
if(!is.null(opt$output_case)){
 if(!grepl("^[A-Za-z0-9_.-]+$",opt$output_case)||opt$output_case %in% c(".",".."))stop("Invalid output case folder")
 RESULT_ROOT<-file.path(PROJECT_ROOT,"results",opt$output_case)
}
for(script in c("data_access.R","plotting.R","dependency_screen.R","lineage_analysis.R","adapters.R","additional_modes.R","tcga_analysis.R","extensions_statistics.R","extension_modules.R","extension_summary.R"))
 source(file.path(PROJECT_ROOT,"scripts/R",script),encoding="UTF-8")
source(file.path(PROJECT_ROOT,"scripts/R/local_data.R"),encoding="UTF-8")
FILES <- depmap_files(PROJECT_ROOT)
check_file <- local_required
if(!is.null(WORKFLOW)) {
 for(script in c("workflow_statistics.R","workflow_plots.R","final_workflow.R"))
  source(file.path(PROJECT_ROOT,"scripts/R",script),encoding="UTF-8")
 run_final_workflow()
 quit(save="no",status=0)
}
dir.create(RESULT_ROOT,recursive=TRUE,showWarnings=FALSE)
for(folder in unname(OUTPUT_DIRS))dir.create(file.path(RESULT_ROOT,folder),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(RESULT_ROOT,"Summary"),showWarnings=FALSE)
run_metadata <- list(geneA=GENE_A,geneB=GENE_B,DepMap_release="26Q1",R_version=as.character(getRversion()),
 thresholds="Analysis-defined; not official DepMap GISTIC: CN_log=log2(relative CN+1), low<0.585, deep<0.35, shallow=[0.35,0.585).",
 core_continuous_variable="CN_log (supplied core); CN-expression uses CN_relative.",
 bootstrap=BOOT_R,min_group_n=MIN_N,seed=1234,export_completeness="All supplied gene columns; full-release completeness unverified.",
 source_sha256=jsonlite::fromJSON(file.path(PROJECT_ROOT,"docs/R_CORE_PROVENANCE.json"))$input_sha256,
 input_manifest_sha256=digest::digest(file=file.path(PROJECT_ROOT,"data/manifests/depmap_26Q1_manifest.csv"),algo="sha256"))
run_metadata$tcga_layers<-list(current="GDC DR46",reference="PanCanAtlas/Xena GISTIC; separate from current")
jsonlite::write_json(run_metadata,file.path(RESULT_ROOT,"Summary/Analysis_Metadata.json"),pretty=TRUE,auto_unbox=TRUE)
fun <- list(qc=run_qc,depmap_cn_expression=run_cn_expression,genomewide_dependency=run_genomewide,
 targeted_dependency=run_targeted,lineage_dependency=run_lineage,adjusted_dependency=run_adjusted,
 reverse_dependency=run_reverse_safe,mutation_dependency=run_mutation,cn_covariation=run_cn_covariation,
 tcga_cn_landscape=run_tcga_landscape,tcga_cna_prevalence=run_tcga_prevalence,tcga_cn_expression=run_tcga_expression,
 cn_threshold_sensitivity=run_threshold_sensitivity,expression_dependency=run_expression_dependency,
 genomewide_expression_dependency=run_genomewide_expression,genomewide_adjusted_dependency=run_genomewide_adjusted)
modes <- if(MODE=="full")setdiff(names(fun),c("genomewide_expression_dependency","genomewide_adjusted_dependency")) else MODE
if(any(!modes %in% names(fun)))stop("Unsupported mode: ",MODE)
msg(R.version.string," | ",GENE_A," -> ",GENE_B," | ",MODE)
for(mode in modes){
 required_inputs<-switch(mode,qc=c("cn","expression","chronos","model"),
  depmap_cn_expression=c("cn","expression","chronos"),genomewide_dependency=c("cn","chronos"),
  targeted_dependency=c("cn","chronos"),lineage_dependency=c("cn","chronos","model"),
  adjusted_dependency=c("cn","chronos","model"),reverse_dependency=c("cn","chronos"),
  mutation_dependency=c("damaging","hotspot","chronos"),cn_covariation="cn",
  cn_threshold_sensitivity=c("cn","chronos"),expression_dependency=c("expression","chronos"),
  genomewide_expression_dependency=c("expression","chronos"),
  genomewide_adjusted_dependency=c("cn","chronos","model"),character())
 for(key in required_inputs)local_required(FILES[[key]])
 start <- Sys.time();state <- "completed";detail <- ""
 tryCatch({fun[[mode]]();finalize_outputs(mode)},error=function(e){state<<-"failed";detail<<-conditionMessage(e)})
 if(!is.null(MODULE_SKIP)){state<-"skipped";detail<-MODULE_SKIP;MODULE_SKIP<-NULL}
 record <- data.table(mode=mode,status=state,seconds=as.numeric(difftime(Sys.time(),start,units="secs")),
  completed_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),detail=detail)
 timing_file <- file.path(RESULT_ROOT,"Summary/Module_Runs.csv")
 fwrite(record,timing_file,append=file.exists(timing_file))
 msg(mode," ",state," (",round(record$seconds,1)," s) ",detail)
 gc()
 if(state=="failed")stop(mode,": ",detail)
}
write_case_summary()
extend_case_summary()
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(RESULT_ROOT,"Summary/SessionInfo.txt"))
msg("Analysis completed.")

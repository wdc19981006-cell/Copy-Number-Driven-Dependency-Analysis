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
option_list <- list(make_option("--mode",default="targeted_dependency"),make_option("--geneA",default="ENO1"),
 make_option("--geneB",default="ENO2"),make_option("--project",default="."),
 make_option("--bootstrap",type="integer",default=1000),make_option("--min_group_n",type="integer",default=3))
opt <- parse_args(OptionParser(option_list=option_list))
MODE <- opt$mode;GENE_A <- toupper(opt$geneA);GENE_B <- toupper(opt$geneB)
BOOT_R <- opt$bootstrap;MIN_N <- opt$min_group_n
stopifnot(BOOT_R>0,MIN_N>=3,grepl("^[A-Z0-9_.-]+$",GENE_A),grepl("^[A-Z0-9_.-]+$",GENE_B))
DEPMAP_DIR <- file.path(PROJECT_ROOT,"data/raw/depmap/26Q1")
RESULT_ROOT <- file.path(PROJECT_ROOT,"results",paste0(GENE_A,"_",GENE_B))
canonical <- c(model="Model.csv",condition="ModelCondition.csv",profiles="OmicsProfiles.csv",cn="CopyNumber_WGS_26Q1.csv",
 expression="Expression_26Q1.csv",chronos="CRISPR_Chronos_26Q1.csv",dependency="CRISPR_GeneDependency_26Q1.csv",
 damaging="Mutation_Damaging_26Q1.csv",hotspot="Mutation_Hotspot_26Q1.csv",signatures="OmicsSignatures_26Q1.csv",subtype="MolecularSubtypes_26Q1.csv")
FILES <- as.list(setNames(file.path(DEPMAP_DIR,canonical),names(canonical)))
for(script in c("data_access.R","plotting.R","dependency_screen.R","lineage_analysis.R","adapters.R","additional_modes.R","tcga_analysis.R"))
 source(file.path(PROJECT_ROOT,"scripts/R",script),encoding="UTF-8")
dir.create(RESULT_ROOT,recursive=TRUE,showWarnings=FALSE)
for(folder in unname(OUTPUT_DIRS))dir.create(file.path(RESULT_ROOT,folder),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(RESULT_ROOT,"Summary"),showWarnings=FALSE)
run_metadata <- list(geneA=GENE_A,geneB=GENE_B,DepMap_release="26Q1",R_version=as.character(getRversion()),
 thresholds="Analysis-defined; not official DepMap GISTIC: CN_log=log2(relative CN+1), low<0.585, deep<0.35, shallow=[0.35,0.585).",
 core_continuous_variable="CN_log (supplied core); CN-expression uses CN_relative.",
 bootstrap=BOOT_R,min_group_n=MIN_N,seed=1234,export_completeness="All supplied gene columns; full-release completeness unverified.",
 source_sha256=jsonlite::fromJSON(file.path(PROJECT_ROOT,"docs/R_CORE_PROVENANCE.json"))$input_sha256,
 input_manifest_sha256=digest::digest(file=file.path(PROJECT_ROOT,"data/manifests/depmap_26Q1_manifest.csv"),algo="sha256"))
jsonlite::write_json(run_metadata,file.path(RESULT_ROOT,"Summary/Analysis_Metadata.json"),pretty=TRUE,auto_unbox=TRUE)
fun <- list(qc=run_qc,depmap_cn_expression=run_cn_expression,genomewide_dependency=run_genomewide,
 targeted_dependency=run_targeted,lineage_dependency=run_lineage,adjusted_dependency=run_adjusted,
 reverse_dependency=run_reverse_safe,mutation_dependency=run_mutation,cn_covariation=run_cn_covariation,
 tcga_cn_landscape=run_tcga_landscape,tcga_cna_prevalence=run_tcga_prevalence,tcga_cn_expression=run_tcga_expression)
modes <- if(MODE=="full")names(fun) else MODE
if(any(!modes %in% names(fun)))stop("Unsupported mode: ",MODE)
msg(R.version.string," | ",GENE_A," -> ",GENE_B," | ",MODE)
for(mode in modes){
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
writeLines(capture.output(sessionInfo()),file.path(RESULT_ROOT,"Summary/SessionInfo.txt"))
msg("Analysis completed.")

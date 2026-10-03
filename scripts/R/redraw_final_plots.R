#!/usr/bin/env Rscript
# Developer repair utility: reuse verified saved statistics, draw affected PDFs
# and refresh just their module caches. No screen or sensitivity recomputation.
expressions<-parse("scripts/R/run_analysis.R")
stopat<-which(vapply(as.list(expressions),function(e)is.call(e)&&identical(e[[1]],as.name("if"))&&
 grepl("run_final_workflow()",paste(deparse(e),collapse=" "),fixed=TRUE),logical(1)))
stopifnot(length(stopat)==1L)
for(e in as.list(expressions)[seq_len(stopat-1L)])eval(e,envir=.GlobalEnv)
if(is.null(WORKFLOW))stop("Supply the existing --workflow and gene arguments")
for(s in c("workflow_statistics.R","workflow_plots.R","final_workflow.R"))source(file.path(PROJECT_ROOT,"scripts/R",s),encoding="UTF-8")
wf_route_supplements()
PATHS<-workflow_paths(WORKFLOW,GENE_A,GENE_B)
TCGA_CANCERS<-sort(unique(unname(unlist(jsonlite::fromJSON(file.path(PROJECT_ROOT,"config/tcga_cancer_types.json"))$mapping))))
dat<-fread(wf_table("TCGA_Reference_Samples.csv"));order<-fread(wf_table("TCGA_Cancer_Order.csv"))
stats<-fread(wf_table(paste0("TCGA_",GENE_A,"_CN_mRNA_AllCancer_Statistics.csv")))
for(cancer in order$CancerType)workflow_save(tcga_expression_plot(dat,cancer,stats[CancerType==cancer]),
 file.path(RESULT_ROOT,PATHS["rna"],paste0(cancer,"_",GENE_A,"_CN_mRNA.pdf")),width=9,height=6.5)
top<-fread(wf_table("Top_CN_Covariation.csv"))
workflow_save(covariation_plot(top),file.path(RESULT_ROOT,PATHS["covariation"]),width=12,height=8)
for(module in c("tcga","cn_covariation")) {
 path<-wf_provenance(paste0("Cache_",module,".json"));cache<-jsonlite::fromJSON(path,simplifyVector=FALSE)
 versions<-wf_input_versions(wf_input_paths(module));code<-wf_code_signature(module)
 cache$code_hash<-code$hash
 params<-wf_parameters(module)
 cache$key<-digest::digest(list(inputs=versions,params=params,code=code$hash),algo="sha256")
 cache$parameters<-params
 cache$outputs<-wf_hash_outputs(unname(wf_expected(module)))
 cache$redraw_only<-TRUE
 jsonlite::write_json(cache,path,pretty=TRUE,auto_unbox=TRUE)
}
metadata<-jsonlite::fromJSON(wf_provenance("Run_Metadata.json"),simplifyVector=FALSE)
metadata$code_files<-setNames(lapply(list.files("scripts/R",pattern="\\.R$",full.names=TRUE),function(f)digest::digest(file=f,algo="sha256")),list.files("scripts/R",pattern="\\.R$",full.names=TRUE))
metadata$plot_only_repair<-"CN-mRNA annotation/inset readability and covariation facet scales; saved statistics unchanged"
jsonlite::write_json(metadata,wf_provenance("Run_Metadata.json"),pretty=TRUE,auto_unbox=TRUE)
wf_summary()
source("scripts/R/validate_final_workflow.R",encoding="UTF-8");validate_final_workflow();wf_summary()
msg("Affected PDFs redrawn from saved statistics; no statistical modules rerun.")

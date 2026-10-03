#!/usr/bin/env Rscript
# Independent acceptance of saved results and cache keys; does not run modules.
expressions<-parse("scripts/R/run_analysis.R")
stopat<-which(vapply(as.list(expressions),function(e)is.call(e)&&identical(e[[1]],as.name("if"))&&
 grepl("run_final_workflow()",paste(deparse(e),collapse=" "),fixed=TRUE),logical(1)))
stopifnot(length(stopat)==1L)
for(e in as.list(expressions)[seq_len(stopat-1L)])eval(e,envir=.GlobalEnv)
for(s in c("workflow_statistics.R","workflow_plots.R","final_workflow.R"))source(file.path(PROJECT_ROOT,"scripts/R",s),encoding="UTF-8")
wf_route_supplements()
PATHS<-workflow_paths(WORKFLOW,GENE_A,GENE_B)
TCGA_CANCERS<-sort(unique(unname(unlist(jsonlite::fromJSON(file.path(PROJECT_ROOT,"config/tcga_cancer_types.json"))$mapping))))
for(module in workflow_modules(WORKFLOW)) {
 old<-jsonlite::fromJSON(wf_provenance(paste0("Cache_",module,".json")),simplifyVector=FALSE)
 key<-digest::digest(list(inputs=wf_input_versions(wf_input_paths(module)),params=wf_parameters(module),
                         code=wf_code_signature(module)$hash),algo="sha256")
 stopifnot(wf_cache_valid(old,key,unname(wf_expected(module))))
 msg("PASS: saved ",module," cache key and output hashes")
}
source("scripts/R/validate_final_workflow.R",encoding="UTF-8");validate_final_workflow()
metadata<-jsonlite::fromJSON(wf_provenance("Run_Metadata.json"),simplifyVector=FALSE)
metadata$code_files<-setNames(lapply(list.files("scripts/R",pattern="\\.R$",full.names=TRUE),function(f)digest::digest(file=f,algo="sha256")),list.files("scripts/R",pattern="\\.R$",full.names=TRUE))
metadata$acceptance<-"All saved module caches verified against final code/parameters/inputs; no analysis modules executed during acceptance"
jsonlite::write_json(metadata,wf_provenance("Run_Metadata.json"),pretty=TRUE,auto_unbox=TRUE)
wf_summary()

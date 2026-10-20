# Maintenance only: after offline numeric equivalence, preserve valid statistics
# and PDFs while updating input manifests and I/O signatures. Never run modules.
expressions<-parse("scripts/R/run_analysis.R")
stopat<-which(vapply(as.list(expressions),function(e)is.call(e)&&identical(e[[1]],as.name("if"))&&
 grepl("run_final_workflow()",paste(deparse(e),collapse=" "),fixed=TRUE),logical(1)))
stopifnot(length(stopat)==1L)
for(e in as.list(expressions)[seq_len(stopat-1L)])eval(e,envir=.GlobalEnv)
for(s in c("workflow_statistics.R","workflow_plots.R","final_workflow.R"))source(file.path(PROJECT_ROOT,"scripts/R",s),encoding="UTF-8")
evidence<-jsonlite::fromJSON("data/manifests/analysis_ready_regression.json",simplifyVector=FALSE)
stopifnot(identical(evidence$status,"PASS"),isTRUE(evidence$all_saved_CSVs_equal),isTRUE(evidence$base_modules_reused))
validated_case<-paste0(evidence$geneA,"_Analysis")
prepared<-jsonlite::fromJSON("data/manifests/analysis_ready_validation.json",simplifyVector=FALSE)
for(p in names(prepared$core))stopifnot(identical(digest::digest(file=p,algo="sha256"),prepared$core[[p]]$SHA256))
TCGA_CANCERS<-sort(unique(unname(unlist(jsonlite::fromJSON("config/tcga_cancer_types.json")$mapping))))
wf_route_supplements()
records<-list();skipped<-list()
normalize_empty_json<-function(x) {
 if(is.list(x)) {
  if(!length(x))return(NULL)
  return(lapply(x,normalize_empty_json))
 }
 x
}
for(case in list.dirs(file.path(PROJECT_ROOT,"results"),recursive=FALSE)) {
 meta_path<-file.path(case,"Provenance/Run_Metadata.json")
 if(!file.exists(meta_path))next
 meta<-jsonlite::fromJSON(meta_path)
 if(is.null(meta$geneA)||!identical(basename(case),paste0(meta$geneA,"_Analysis")))next
 RESULT_ROOT<-case;GENE_A<-meta$geneA
 for(b in c("",wf_geneBs())) {
  GENE_B<-b;GENE_B_PROVIDED<-nzchar(b);WORKFLOW<-if(nzchar(b))"geneA_geneB" else "geneA_screen"
  PATHS<-workflow_paths(WORKFLOW,GENE_A,GENE_B)
  modules<-if(nzchar(b))wf_geneB_modules() else wf_base_modules()
  for(module in modules) {
   cache_path<-wf_module_provenance(module,paste0("Cache_",module,".json"))
   if(!file.exists(cache_path))next
   old<-jsonlite::fromJSON(cache_path,simplifyVector=FALSE)
   outputs<-unname(wf_expected(module))
   # Missing/corrupt results stay invalid; never bless them with a new key.
   if(!wf_cache_valid(old,old$key,outputs)) {
    if(basename(case)==validated_case)stop("Validated case cache outputs invalid: ",cache_path)
    skipped[[length(skipped)+1L]]<-list(case=basename(case),geneB=b,module=module,
     reason="Pre-existing cache output names/hashes do not match current layout; retained unchanged and not declared equivalent")
    next
   }
   versions<-wf_input_versions(wf_input_paths(module))
   identical_inputs<-TRUE
   for(prior in old$data_versions) {
    hit<-Filter(function(x)identical(x$path,prior$path),versions)
    if(length(hit)!=1L||!isTRUE(all.equal(normalize_empty_json(hit[[1]]),
      normalize_empty_json(prior),tolerance=0,check.attributes=FALSE)))identical_inputs<-FALSE
   }
   if(!identical_inputs) {
    if(basename(case)==validated_case)stop("Validated case inputs differ: ",cache_path)
    skipped[[length(skipped)+1L]]<-list(case=basename(case),geneB=b,module=module,
     reason="Pre-existing input metadata mismatch; retained unchanged and not declared equivalent")
    next
   }
   params<-wf_parameters(module)
   if(!identical(old$parameter_hash,digest::digest(params,algo="sha256"))) {
    skipped[[length(skipped)+1L]]<-list(case=basename(case),geneB=b,module=module,
     reason="Existing cache parameter hash predates current framework; retained unchanged and not declared equivalent")
    next
   }
   code<-wf_code_signature(module)
   key<-digest::digest(list(inputs=versions,params=params,code=code$hash),algo="sha256")
   records[[length(records)+1L]]<-list(case=basename(case),geneB=b,module=module,
      statistics_recomputed=FALSE,outputs_unchanged=TRUE,previous_key=old$key,key=key,
      previous_code_hash=old$code_hash,code_hash=code$hash,data_versions=versions)
  }
 }
}
jsonlite::write_json(list(validation_SHA256=digest::digest(file="data/manifests/analysis_ready_regression.json",algo="sha256"),
 records=records,retained_without_migration=skipped),"data/manifests/analysis_storage_cache_migration.json",pretty=TRUE,auto_unbox=TRUE)
cat("PASS: migrated",length(records),"verified module caches; retained",length(skipped),"older parameter caches unchanged; no statistics or plots executed.\n")

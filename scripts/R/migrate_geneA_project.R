# Explicit maintenance only: consolidate a legacy A+B case without running analysis.
migrate_geneA_project <- function(previous_source_dir,legacy_root) {
 stopifnot(as.character(getRversion())=="4.5.0",GENE_B_PROVIDED,
  normalizePath(legacy_root,winslash="/")!=normalizePath(RESULT_ROOT,winslash="/"))
 original_root<-RESULT_ROOT
 legacy<-new.env(parent=.GlobalEnv)
 sys.source(file.path(previous_source_dir,"final_workflow.R"),legacy)
 sys.source(file.path(previous_source_dir,"workflow_plots.R"),legacy)
 # Reproduce the exact, formerly signed supplementary routing and CSV adapter.
 legacy$wf_route_supplements()
 fwrite<<-function(x,file,...) {
  if(grepl("/Supplementary/",gsub("\\\\","/",file),fixed=TRUE)&&grepl("\\.csv$",file))file<-wf_table(basename(file))
  data.table::fwrite(x,file,...,scipen=0)
 }
 baseline<-list();prepared<-list()
 for(which in c("base","pair")) {
  legacy$RESULT_ROOT<-if(which=="base")RESULT_ROOT else legacy_root
  legacy$WORKFLOW<-if(which=="base")"geneA_screen" else "geneA_geneB"
  legacy$GENE_B_PROVIDED<-which=="pair";legacy$GENE_B<-if(which=="pair")GENE_B else ""
  legacy$PATHS<-legacy$workflow_paths(legacy$WORKFLOW,GENE_A,legacy$GENE_B)
  for(module in workflow_modules(legacy$WORKFLOW)) {
   path<-file.path(legacy$RESULT_ROOT,"Provenance",paste0("Cache_",module,".json"))
   old<-jsonlite::fromJSON(path,simplifyVector=FALSE)
   versions<-legacy$wf_input_versions(legacy$wf_input_paths(module))
   params<-legacy$wf_parameters(module);code<-legacy$wf_code_signature(module)
   key<-digest::digest(list(inputs=versions,params=params,code=code$hash),algo="sha256")
   if(!legacy$wf_cache_valid(old,key,unname(legacy$wf_expected(module))))
    stop("Legacy cache verification failed: ",which," / ",module,". No analysis was run.")
   if(which=="base")baseline[[module]]<-old
   else prepared[[module]]<-old
  }
 }
 # Compare every duplicated statistical/sample CSV by SHA256 AND parsed content.
 comparisons<-list()
 for(module in wf_base_modules()) {
  a<-names(baseline[[module]]$outputs);b<-names(prepared[[module]]$outputs)
  a<-a[grepl("\\.csv$",a)];b<-b[grepl("\\.csv$",b)]
  for(p in a) {
   q<-b[basename(b)==basename(p)];stopifnot(length(q)==1L)
   fa<-file.path(RESULT_ROOT,p);fb<-file.path(legacy_root,q)
   ha<-digest::digest(file=fa,algo="sha256");hb<-digest::digest(file=fb,algo="sha256")
   stopifnot(identical(ha,hb),identical(data.table::fread(fa),data.table::fread(fb)))
   comparisons[[p]]<-list(legacy_path=q,sha256=ha,content_identical=TRUE)
  }
 }
 old_files<-sort(list.files(legacy_root,recursive=TRUE))
 hashes<-setNames(lapply(old_files,function(p)digest::digest(file=file.path(legacy_root,p),algo="sha256")),old_files)
 base_files<-sort(list.files(file.path(RESULT_ROOT,"Main_Results"),recursive=TRUE,full.names=TRUE))
 base_hashes<-setNames(lapply(base_files,function(p)digest::digest(file=p,algo="sha256")),base_files)
 # Store the prior base narrative before any Gene B entry is appended.
 GENE_B_PROVIDED<<-FALSE;wf_summary();GENE_B_PROVIDED<<-TRUE
 wf_create_directories()
 mapping<-setNames(rep(NA_character_,length(old_files)),old_files)
 legacy$PATHS<-legacy$workflow_paths("geneA_geneB",GENE_A,GENE_B)
 for(module in wf_geneB_modules()) {
  from<-unname(legacy$wf_expected(module));to<-unname(wf_expected(module))
  stopifnot(identical(basename(from),basename(to)),length(from)==length(to))
  mapping[from]<-to
 }
 # Preserve legacy provenance/summary, including original caches and validation.
 archive<-file.path("Provenance/GeneB",GENE_B,"Legacy")
 for(p in old_files[startsWith(old_files,"Provenance/")&!startsWith(old_files,"Provenance/Data_Audit/")])
  mapping[p]<-file.path(archive,sub("^Provenance/","",p))
 mapping["00_Analysis_Summary.txt"]<-file.path(archive,"00_Analysis_Summary.txt")
 for(p in names(mapping)[!is.na(mapping)]) {
  target<-file.path(RESULT_ROOT,mapping[[p]]);dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  if(file.exists(target))stop("Migration destination exists: ",target)
  stopifnot(file.copy(file.path(legacy_root,p),target),identical(hashes[[p]],digest::digest(file=target,algo="sha256")))
 }
 # Duplicated PDFs may differ in timestamps or Gene B volcano highlighting.
 # Record both hashes; retain the original Gene A discovery PDFs unchanged.
 pdf_comparisons<-list()
 for(module in wf_base_modules()) {
  a<-names(baseline[[module]]$outputs);b<-names(prepared[[module]]$outputs)
  a<-a[grepl("\\.pdf$",a)];b<-b[grepl("\\.pdf$",b)]
  for(i in seq_along(a))pdf_comparisons[[a[i]]]<-list(sha256_base=baseline[[module]]$outputs[[a[i]]],
   sha256_legacy=prepared[[module]]$outputs[[b[i]]],statistics_identical=TRUE,
   retained="Existing Gene A PDF; no regeneration. Legacy PDF may include timestamps / candidate highlighting.")
 }
 # Now activate new routing. Restore original adapter definitions first.
 source("scripts/R/lineage_analysis.R",encoding="UTF-8")
 source("scripts/R/adapters.R",encoding="UTF-8")
 wf_route_supplements()
 fwrite<<-function(x,file,...) {
  if(grepl("/Supplementary/",gsub("\\\\","/",file),fixed=TRUE)&&grepl("\\.csv$",file))file<-wf_table(basename(file))
  data.table::fwrite(x,file,...,scipen=0)
 }
 pair<-data.table::fread(wf_table(paste0(GENE_A,"_",GENE_B,"_CellLines.csv")))
 s<-data.table::fread(wf_table(paste0(GENE_A,"_",GENE_B,"_Targeted_Statistics.csv")))
 # The only redraw in this maintenance operation; all CSV bytes stay intact.
 workflow_save(targeted_scatter_plot(pair,s),file.path(RESULT_ROOT,PATHS["scatter"]))
 for(module in workflow_modules(WORKFLOW)) {
  old<-if(module %in% wf_base_modules())baseline[[module]] else prepared[[module]]
  code<-wf_code_signature(module);params<-wf_parameters(module);versions<-wf_input_versions(wf_input_paths(module))
  old$previous_key<-old$key;old$previous_code_hash<-old$code_hash
  old$parameters<-params;old$parameter_hash<-digest::digest(params,algo="sha256")
  old$key<-digest::digest(list(inputs=versions,params=params,code=code$hash),algo="sha256")
  old$code_hash<-code$hash;old$code_functions<-code$functions
  old$outputs<-wf_hash_outputs(unname(wf_expected(module)))
  if(module=="targeted_dependency")old$statistics_key<-wf_targeted_statistics_key()
  stopifnot(wf_cache_valid(old,old$key,unname(wf_expected(module))))
  jsonlite::write_json(old,wf_module_provenance(module,paste0("Cache_",module,".json")),pretty=TRUE,auto_unbox=TRUE)
 }
 preserved<-names(mapping)[!is.na(mapping)&!grepl("CN_vs_.*_Dependency\\.pdf$",mapping)]
 stopifnot(all(vapply(preserved,function(p)identical(hashes[[p]],digest::digest(file=file.path(RESULT_ROOT,mapping[[p]]),algo="sha256")),logical(1))),
  identical(base_hashes,setNames(lapply(base_files,function(p)digest::digest(file=p,algo="sha256")),base_files)))
 record<-list(status="PASS_AWAITING_VISUAL_REVIEW",legacy_root=normalizePath(legacy_root,winslash="/"),
  statistics_recomputed=FALSE,redrawn_PDF=unname(PATHS["scatter"]),base_outputs_unchanged=TRUE,
  duplicated_statistics=comparisons,duplicated_PDFs=pdf_comparisons,
  files=lapply(old_files,function(p)list(before=p,after=mapping[[p]],sha256_before=hashes[[p]],
   disposition=if(is.na(mapping[[p]]))"Duplicate Gene A artifact; canonical statistics compared" else "Copied and SHA256 verified")))
 jsonlite::write_json(record,wf_provenance("GeneA_Project_Migration.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
 wf_summary()
 source("scripts/R/validate_final_workflow.R",encoding="UTF-8");validate_result_layout()
 message("PASS: consolidated without statistical execution; only scatter redrawn. Legacy retained for visual review.")
 invisible(record)
}
finalize_geneA_migration <- function() {
 wf_route_supplements()
 report_path<-wf_provenance("GeneA_Project_Migration.json")
 record<-jsonlite::fromJSON(report_path,simplifyVector=FALSE)
 old<-normalizePath(record$legacy_root,winslash="/",mustWork=FALSE)
 results<-normalizePath(file.path(PROJECT_ROOT,"results"),winslash="/",mustWork=TRUE)
 stopifnot(dirname(old)==results,basename(old)==paste0(GENE_A,"_",GENE_B,"_Analysis"),
           old!=normalizePath(RESULT_ROOT,winslash="/",mustWork=TRUE))
 for(item in record$files) {
  original<-file.path(old,item$before)
  if(file.exists(original))stopifnot(identical(item$sha256_before,digest::digest(file=original,algo="sha256"))) else
   stopifnot(identical(record$status,"DELETION_PENDING")||isTRUE(record$legacy_removed))
  if(is.character(item$after)&&length(item$after)==1L&&!identical(item$after,unname(PATHS["scatter"])))
   stopifnot(identical(item$sha256_before,digest::digest(file=file.path(RESULT_ROOT,item$after),algo="sha256")))
 }
 source("scripts/R/validate_final_workflow.R",encoding="UTF-8");validate_result_layout()
 for(module in workflow_modules(WORKFLOW)) {
  cache<-jsonlite::fromJSON(wf_module_provenance(module,paste0("Cache_",module,".json")),simplifyVector=FALSE)
  key<-digest::digest(list(inputs=wf_input_versions(wf_input_paths(module)),params=wf_parameters(module),
   code=wf_code_signature(module)$hash),algo="sha256")
  stopifnot(wf_cache_valid(cache,key,unname(wf_expected(module))))
 }
 record$status<-"DELETION_PENDING";record$visual_review<-"PASS";record$legacy_removed<-FALSE
 record$verified_caches<-workflow_modules(WORKFLOW)
 validation<-jsonlite::fromJSON(wf_provenance("Validation.json"))
 record$independent_statistics_validation<-list(status=validation$status,checks=validation$number_of_checks,
  report="Provenance/Validation.json",analysis_modules_run=FALSE)
 jsonlite::write_json(record,report_path,pretty=TRUE,auto_unbox=TRUE)
 # The resolved, explicitly named sibling directory has been verified above.
 unlink(old,recursive=TRUE,force=TRUE)
 if(dir.exists(old))stop("Legacy directory removal failed; verified migrated results are intact: ",old)
 record$status<-"PASS";record$legacy_removed<-TRUE
 jsonlite::write_json(record,report_path,pretty=TRUE,auto_unbox=TRUE)
 wf_summary();validate_result_layout()
 message("PASS: legacy directory removed after artifact and visual verification.")
}
refresh_saved_dependency_scatter <- function() {
 wf_route_supplements()
 module<-"targeted_dependency";path<-wf_module_provenance(module,"Cache_targeted_dependency.json")
 old<-jsonlite::fromJSON(path,simplifyVector=FALSE)
 stopifnot(identical(old$statistics_key,wf_targeted_statistics_key()),identical(wf_run_targeted(),"plot_regenerated"))
 code<-wf_code_signature(module);params<-wf_parameters(module);versions<-wf_input_versions(wf_input_paths(module))
 old$code_hash<-code$hash;old$code_functions<-code$functions
 old$key<-digest::digest(list(inputs=versions,params=params,code=code$hash),algo="sha256")
 old$outputs<-wf_hash_outputs(unname(wf_expected(module)))
 jsonlite::write_json(old,path,pretty=TRUE,auto_unbox=TRUE)
 wf_summary()
 message("PASS: only the saved scatter regenerated; CSVs and other PDFs reused.")
}
if(sys.nframe()==0L) {
 expressions<-as.list(parse("scripts/R/run_analysis.R"))
 stopat<-which(vapply(expressions,function(e)is.call(e)&&identical(e[[1]],as.name("if"))&&
  grepl("run_final_workflow()",paste(deparse(e),collapse=" "),fixed=TRUE),logical(1)))
 for(e in expressions[seq_len(stopat-1L)])eval(e,envir=.GlobalEnv)
 for(s in c("workflow_statistics.R","workflow_plots.R","final_workflow.R"))source(file.path("scripts/R",s),encoding="UTF-8")
 PATHS<-workflow_paths(WORKFLOW,GENE_A,GENE_B)
 TCGA_CANCERS<-sort(unique(unname(unlist(jsonlite::fromJSON("config/tcga_cancer_types.json")$mapping))))
 if(identical(Sys.getenv("MIGRATION_ACTION"),"redraw"))refresh_saved_dependency_scatter() else
 if(identical(Sys.getenv("MIGRATION_VISUAL_REVIEW"),"PASS"))finalize_geneA_migration() else
  migrate_geneA_project(Sys.getenv("PREVIOUS_WORKFLOW_DIR"),file.path(PROJECT_ROOT,"results",paste0(GENE_A,"_",GENE_B,"_Analysis")))
}

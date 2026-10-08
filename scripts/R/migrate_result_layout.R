# Path-only migration of an existing final-workflow case. No analysis entry point
# or data-matrix reader is called. Keep the previous workflow source outside results.
wf_layout_previous_code <- function(module,previous_file_sources,envir=.GlobalEnv) {
 # Replay only a reviewed path-only helper edit when the old workflow signature
 # included that file. All other functions/files are hashed normally.
 if(length(previous_file_sources)) {
  stopifnot(identical(names(previous_file_sources),"scripts/utils/export_workflow_tcga.py"))
  before<-paste(readLines(previous_file_sources[[1]],warn=FALSE),collapse="\n")
  after<-paste(readLines(file.path(PROJECT_ROOT,names(previous_file_sources)),warn=FALSE),collapse="\n")
  old_path<-"target.parent.parent / 'Provenance/TCGA_Baseline_Method.json'"
  new_path<-"ROOT / 'results' / target.relative_to((ROOT / 'results').resolve()).parts[0] / 'Provenance/TCGA_Baseline_Method.json'"
  stopifnot(grepl(old_path,before,fixed=TRUE),identical(gsub(old_path,new_path,before,fixed=TRUE),after))
 }
 scope<-new.env(parent=envir)
 scope$file.path<-function(...) {
  args<-list(...)
  if(length(args)==2L&&identical(args[[1]],PROJECT_ROOT)&&args[[2]] %in% names(previous_file_sources))
   return(previous_file_sources[[args[[2]]]])
  base::file.path(...)
 }
 f<-get("wf_code_signature",envir=envir);environment(f)<-scope
 f(module)
}
wf_refresh_layout_cache <- function(previous_file_sources) {
 # For a case already migrated before a discovered helper-path correction.
 # Verify the old key against the recorded cache before refreshing its signature.
 module<-"tcga";path<-wf_provenance("Cache_tcga.json")
 cache<-jsonlite::fromJSON(path,simplifyVector=FALSE)
 old_code<-wf_layout_previous_code(module,previous_file_sources)
 versions<-wf_input_versions(wf_input_paths(module));params<-wf_parameters(module)
 old_key<-digest::digest(list(inputs=versions,params=params,code=old_code$hash),algo="sha256")
 stopifnot(identical(cache$code_hash,old_code$hash),wf_cache_valid(cache,old_key,unname(wf_expected(module))))
 code<-wf_code_signature(module)
 cache$key<-digest::digest(list(inputs=versions,params=params,code=code$hash),algo="sha256")
 cache$code_hash<-code$hash
 proof<-list(previous_code_hash=old_code$hash,previous_key=old_key,
  file=names(previous_file_sources),before_sha256=digest::digest(file=previous_file_sources[[1]],algo="sha256"),
  after_sha256=digest::digest(file=file.path(PROJECT_ROOT,names(previous_file_sources)),algo="sha256"),
  proof="Exact source comparison: only the baseline JSON destination changed; no sample/statistical expression changed.")
 cache$layout_migration$output_path_helper<-proof
 stopifnot(wf_cache_valid(cache,cache$key,unname(wf_expected(module))))
 jsonlite::write_json(cache,path,pretty=TRUE,auto_unbox=TRUE)
 record_path<-wf_provenance("Result_Layout_Migration.json")
 record<-jsonlite::fromJSON(record_path,simplifyVector=FALSE);record$output_path_helper<-proof
 jsonlite::write_json(record,record_path,pretty=TRUE,auto_unbox=TRUE)
 invisible(TRUE)
}
migrate_result_layout <- function(previous_workflow,previous_file_sources=list()) {
 stopifnot(as.character(getRversion())=="4.5.0",dir.exists(file.path(RESULT_ROOT,"Tables")))
 legacy<-new.env(parent=.GlobalEnv);sys.source(previous_workflow,envir=legacy)
 legacy$PATHS<-legacy$workflow_paths(WORKFLOW,GENE_A,GENE_B)
 old_files<-sort(list.files(RESULT_ROOT,recursive=TRUE))
 old_hashes<-wf_hash_outputs(old_files)
 # Remove only the new metadata-writing statements for a source-text proof that
 # every pre-existing statistical/plot function is unchanged.
 strip_layout_writes<-function(x) {
  if(missing(x))return(quote(expr=))
  if(!is.call(x))return(x)
  if(identical(x[[1]],as.name("wf_write_screen_status")))return(NULL)
  if(identical(x[[1]],as.name("wf_write"))&&is.call(x[[2]])&&identical(x[[2]][[1]],as.name("wf_group_statistics")))return(NULL)
  parts<-lapply(as.list(x),strip_layout_writes)
  if(identical(x[[1]],as.name("{")))parts<-Filter(Negate(is.null),parts)
  as.call(parts)
 }
 text<-function(f)paste(deparse(body(f),width.cutoff=500L),collapse="\n")
 modules<-workflow_modules(WORKFLOW);prepared<-list()
 for(module in modules) {
  old<-jsonlite::fromJSON(wf_provenance(paste0("Cache_",module,".json")),simplifyVector=FALSE)
  outputs<-unname(legacy$wf_expected(module))
  versions<-wf_input_versions(wf_input_paths(module));params<-wf_parameters(module)
  old_code<-wf_layout_previous_code(module,previous_file_sources,envir=legacy)
  old_key<-digest::digest(list(inputs=versions,params=params,code=old_code$hash),algo="sha256")
  if(!identical(old_code$hash,old$code_hash)||!wf_cache_valid(old,old_key,outputs))
   stop("Cannot safely migrate cache: ",module,"; previous code, inputs or output hashes differ. No analysis was run.")
  for(n in setdiff(old_code$functions,"workflow_paths")) {
   a<-get(n,envir=legacy);b<-get(n,envir=.GlobalEnv)
   if(n %in% c("wf_run_screen","wf_run_targeted"))body(b)<-strip_layout_writes(body(b))
   stopifnot(identical(formals(a),formals(b)),identical(text(a),text(b)))
  }
  # Routing may change directories only, never the existing artifact basenames.
  new_outputs<-unname(wf_expected(module))
  stopifnot(all(basename(outputs) %in% basename(new_outputs)))
  new_code<-wf_code_signature(module)
  prepared[[module]]<-list(old=old,code=new_code,
   key=digest::digest(list(inputs=versions,params=params,code=new_code$hash),algo="sha256"))
 }
 map<-setNames(old_files,old_files)
 for(p in old_files) {
  if(startsWith(p,"Tables/"))map[[p]]<-wf_table_relative(basename(p))
  if(startsWith(p,"Main_Results/")&&tools::file_ext(p)=="pdf") {
   i<-match(p,legacy$PATHS)
   if(!is.na(i))map[[p]]<-unname(PATHS[i])
  }
 }
 stopifnot(!anyDuplicated(unname(map)))
 targets<-names(map)[map!=names(map)]
 stopifnot(!any(file.exists(file.path(RESULT_ROOT,map[targets]))))
 # Persist the complete pre-move record before the first rename.
 record<-list(status="RECORDED_BEFORE_MOVE",previous_workflow_sha256=digest::digest(file=previous_workflow,algo="sha256"),
  old_file_index_sha256=old_hashes[["Provenance/File_Index.csv"]],analysis_recomputed=FALSE,PDF_redrawn=FALSE,
  index_exception="File_Index.csv is intentionally rewritten to reflect relocated files; result/sample CSVs and Module_Runs.csv are byte-preserved.",
  files=lapply(old_files,function(p)list(before=p,after=map[[p]],sha256_before=old_hashes[[p]])))
 jsonlite::write_json(record,wf_provenance("Result_Layout_Migration.json"),pretty=TRUE,auto_unbox=TRUE)
 wf_create_directories()
 for(p in targets)if(!file.rename(file.path(RESULT_ROOT,p),file.path(RESULT_ROOT,map[[p]])))stop("Cannot move result: ",p)
 # Independent re-read of bytes after every file has reached its destination.
 preserved<-setdiff(old_files,"Provenance/File_Index.csv")
 after_hashes<-wf_hash_outputs(unname(map[preserved]))
 stopifnot(identical(unname(old_hashes[preserved]),unname(after_hashes)))
 stopifnot(length(list.files(file.path(RESULT_ROOT,"Tables"),all.files=TRUE,no..=TRUE))==0L)
 stopifnot(dirname(normalizePath(file.path(RESULT_ROOT,"Tables"),winslash="/"))==normalizePath(RESULT_ROOT,winslash="/"),
  startsWith(normalizePath(RESULT_ROOT,winslash="/"),paste0(normalizePath(file.path(PROJECT_ROOT,"results"),winslash="/"),"/")))
 unlink(file.path(RESULT_ROOT,"Tables"),recursive=TRUE)
 wf_write_screen_status(jsonlite::fromJSON(wf_provenance("GenomeWide_Dependency_Status.json")))
 source_metadata<-readLines(wf_provenance("Source_Metadata.txt"),warn=FALSE)
 writeLines(gsub("audit in Tables.","audit in Provenance/Data_Audit.",source_metadata,fixed=TRUE),wf_provenance("Source_Metadata.txt"))
 if(GENE_B_PROVIDED)wf_write(wf_group_statistics(data.table::fread(wf_table(paste0(GENE_A,"_",GENE_B,"_Targeted_Statistics.csv")))),
  paste0(GENE_A,"_",GENE_B,"_CNlow_vs_Normal_Statistics.csv"))
 for(module in modules) {
  item<-prepared[[module]];cache<-item$old
  cache$key<-item$key;cache$code_hash<-item$code$hash;cache$code_functions<-item$code$functions
  cache$outputs<-wf_hash_outputs(unname(wf_expected(module)))
  cache$layout_migration<-list(previous_key=item$old$key,previous_code_hash=item$old$code_hash,
   proof="Previous code/input key and output hashes verified; signed statistical and plot bodies unchanged; old artifact basenames and bytes preserved.")
  stopifnot(wf_cache_valid(cache,item$key,unname(wf_expected(module))))
  jsonlite::write_json(cache,wf_provenance(paste0("Cache_",module,".json")),pretty=TRUE,auto_unbox=TRUE)
 }
 record$status<-"PASS";record$modules<-modules;record$all_original_result_PDF_CSV_byte_preserved<-TRUE
 record$preserved_PDFs<-sum(tools::file_ext(preserved)=="pdf")
 record$preserved_CSVs<-sum(tools::file_ext(preserved)=="csv")
 for(i in seq_along(record$files)) {
  p<-record$files[[i]]$before
  # Cache metadata and summary are intentionally updated after artifact checks.
  record$files[[i]]$sha256_after_move<-if(p=="Provenance/File_Index.csv")NULL else after_hashes[[map[[p]]]]
 }
 jsonlite::write_json(record,wf_provenance("Result_Layout_Migration.json"),pretty=TRUE,auto_unbox=TRUE)
 wf_summary()
 validate_result_layout()
 # Preserve the original independent statistical validation instead of claiming
 # that it was rerun. Append the separately executed structural migration audit.
 validation<-jsonlite::fromJSON(wf_provenance("Validation.json"),simplifyVector=FALSE)
 validation$result_layout_migration<-list(status="PASS",statistics_recomputed=FALSE,PDFs_redrawn=FALSE,
  artifact_hashes_preserved=TRUE,exact_file_index=TRUE,exact_33_cancer_PDFs=TRUE,
  caches_verified=modules,report="Provenance/Result_Layout_Migration.json")
 jsonlite::write_json(validation,wf_provenance("Validation.json"),pretty=TRUE,auto_unbox=TRUE)
 # Re-read serialized caches, then independently recompute the current key.
 for(module in modules) {
  cache<-jsonlite::fromJSON(wf_provenance(paste0("Cache_",module,".json")),simplifyVector=FALSE)
  code<-wf_code_signature(module)
  key<-digest::digest(list(inputs=wf_input_versions(wf_input_paths(module)),params=wf_parameters(module),code=code$hash),algo="sha256")
  stopifnot(wf_cache_valid(cache,key,unname(wf_expected(module))))
 }
 validate_result_layout()
 message("PASS: path-only migration; all result PDF/CSV bytes and verified caches preserved.")
 invisible(record)
}

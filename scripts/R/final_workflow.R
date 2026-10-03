# Orchestration, module caches and the final user-facing result contract.
workflow_paths <- function(workflow,geneA,geneB=NULL) {
 pair<-identical(workflow,"geneA_geneB")
 c(cna=paste0("Main_Results/01_TCGA_",geneA,"_CNA_Percentage.pdf"),
 landscape=paste0("Main_Results/02_TCGA_",geneA,"_CopyNumber_Landscape.pdf"),
 rna="Main_Results/03_TCGA_CN_mRNA",
 if(pair)c(scatter=paste0("Main_Results/04_DepMap_",geneA,"_CN_vs_",geneB,"_Dependency.pdf"),
 waterfall=paste0("Main_Results/05_DepMap_",geneB,"_Dependency_Waterfall.pdf"),
 box=paste0("Main_Results/06_DepMap_",geneB,"_Dependency_CNlow_vs_Normal.pdf")),
 volcano=paste0("Main_Results/",if(pair)"07" else "04","_DepMap_",geneA,"_GenomeWide_Dependency_Volcano.pdf"),
 covariation=paste0("Main_Results/",if(pair)"08" else "05","_DepMap_",geneA,"_CN_Covariation.pdf"))
}
wf_table <- function(name)file.path(RESULT_ROOT,"Tables",name)
wf_provenance <- function(name)file.path(RESULT_ROOT,"Provenance",name)
wf_write <- function(x,name) data.table::fwrite(x,wf_table(name),scipen=0)
wf_route_supplements <- function() {
 # adapters.R had already substituted legacy folder literals in these supplied
 # functions. Route those literals to final supplementary folders at runtime.
 for(nm in c("run_lineage","run_adjusted")) {
  f<-get(nm);code<-paste(deparse(body(f),width.cutoff=500L),collapse="\n")
  for(pair in list(c("04_Lineage_Dependency","Supplementary/01_Lineage"),
                  c("05_Adjusted_Dependency","Supplementary/02_Adjusted"),
                  c("lineage_dependency","Supplementary/01_Lineage"),
                  c("adjusted_dependency","Supplementary/02_Adjusted")))
   code<-gsub(paste0('"',pair[1],'"'),paste0('"',pair[2],'"'),code,fixed=TRUE)
  body(f)<-parse(text=code)[[1]];assign(nm,f,envir=.GlobalEnv)
 }
}
wf_input_paths <- function(module) {
 ref<-file.path(PROJECT_ROOT,"data/processed/tcga/pancanatlas_reference")
 switch(module,tcga=c(file.path(ref,paste0(c("CN","GISTIC","Expression","Metadata"),".parquet")),
                      file.path(ref,paste0(rep(c("CN","GISTIC","Expression"),each=2),c(".genes.parquet",".samples.parquet")))),
 genomewide_dependency=c(FILES$cn,FILES$chronos),cn_covariation=FILES$cn,
 c(FILES$cn,FILES$chronos,if(module %in% c("lineage_dependency","adjusted_dependency"))FILES$model))
}
wf_input_versions <- function(paths) {
 lapply(paths,function(p) {
  local_required(p);s<-file.info(p);side<-paste0(p,".provenance.json")
  meta<-if(file.exists(side))jsonlite::fromJSON(side,simplifyVector=FALSE) else NULL
  list(path=substring(p,nchar(PROJECT_ROOT)+2),bytes=s$size,mtime=format(s$mtime,"%Y-%m-%dT%H:%M:%OS6%z"),
       prepared_sha256=if(!is.null(meta$SHA256))meta$SHA256 else meta$output_sha256,
       provenance=meta)
 })
}
wf_code_signature <- function(module) {
 functions<-switch(module,
 tcga=c("cancer_order","ordered_cancer","reference_cancer_statistics","reference_prevalence",
        "tcga_prevalence_plot","tcga_landscape_plot","tcga_expression_plot","wf_run_tcga","pair_correlations"),
 targeted_dependency=c("targeted_statistics","targeted_plots","waterfall_order","wf_run_targeted","pair_correlations","group_effect"),
 genomewide_dependency=c("run_genomewide","screen_core_statistics","eligible_rank","dependency_top","volcano_plot","wf_run_screen"),
 cn_covariation=c("covariation_statistics","covariation_top","covariation_plot","wf_run_covariation"),
 lineage_dependency=c("run_lineage","wf_run_supplement"),adjusted_dependency=c("run_adjusted","wf_run_supplement","coefficient_table"),
 cn_threshold_sensitivity=c("run_threshold_sensitivity","cn_threshold_statistics","group_effect","wf_run_supplement"))
 common<-c("prepare_cn","load_target_pair","read_gene","load_model","fread","find_gene_col","clean_gene",
           "processed_depmap","processed_pair","wf_write","workflow_save","workflow_theme","save_pdf",
           "wf_expected","workflow_paths","wf_validate_artifact","p_text","p_star")
 # Hash source text rather than serialized language objects: runtime/JIT or
 # package expression caches can change attributes of an otherwise identical
 # function body after plotting. Text stays stable across cold and warm calls.
 code<-lapply(unique(c(functions,common)),function(n)list(name=n,
  args=paste(deparse(formals(get(n)),width.cutoff=500L),collapse="\n"),
  body=paste(deparse(body(get(n)),width.cutoff=500L),collapse="\n")))
 py<-c("scripts/R/local_data.R","config/workflows.json")
 if(module=="tcga")py<-c(py,"scripts/utils/export_workflow_tcga.py","scripts/data_access.py","config/tcga_cancer_types.json")
 list(hash=digest::digest(list(code=code,files=lapply(py,function(p)digest::digest(file=file.path(PROJECT_ROOT,p),algo="sha256")),
               constants=list(states=CNA_STATES,colors=CNA_COLORS,tcga_min_n=TCGA_MIN_N)),algo="sha256"),functions=unique(c(functions,common)))
}
wf_expected <- function(module) {
 switch(module,
 tcga=c(PATHS[c("cna","landscape")],paste0(PATHS["rna"],"/",TCGA_CANCERS,"_",GENE_A,"_CN_mRNA.pdf"),
        "Tables/TCGA_Cancer_Order.csv","Tables/TCGA_CNA_Percentage.csv","Tables/TCGA_Reference_Samples.csv",
        paste0("Tables/TCGA_",GENE_A,"_CN_mRNA_AllCancer_Statistics.csv"),"Provenance/TCGA_Plot_Order.json","Tables/TCGA_CN_mRNA_CNA_Counts.csv"),
 targeted_dependency=c(PATHS[c("scatter","waterfall","box")],paste0("Tables/",GENE_A,"_",GENE_B,
          c("_Targeted_Statistics.csv","_CellLines.csv","_Waterfall_Order.csv"))),
 genomewide_dependency=c(PATHS["volcano"],"Tables/GenomeWide_Dependency.csv","Tables/Top_Dependency_Candidates.csv"),
 cn_covariation=c(PATHS["covariation"],paste0("Tables/",GENE_A,"_CN_Covariation.csv"),"Tables/Top_CN_Covariation.csv"),
 lineage_dependency=c("Supplementary/01_Lineage/Lineage_Dependency_Forest.pdf","Tables/Lineage_Dependency.csv","Tables/Lineage_Eligibility.csv"),
 adjusted_dependency=c("Supplementary/02_Adjusted/Adjusted_CN_Coefficient.pdf","Tables/Continuous_CN_Adjusted.csv","Tables/CNLow_Adjusted.csv"),
 cn_threshold_sensitivity=c("Supplementary/03_CN_Threshold_Sensitivity/CN_Threshold_Sensitivity.pdf",
                            "Tables/CN_Threshold_Sensitivity.csv","Tables/CN_Threshold_Continuous_Statistics.csv"))
}
wf_validate_artifact <- function(path) {
 if(!file.exists(path)||is.na(file.info(path)$size)||file.info(path)$size==0)return(FALSE)
 ext<-tolower(tools::file_ext(path))
 if(ext=="pdf") {
  con<-file(path,"rb");on.exit(close(con));return(identical(rawToChar(readBin(con,"raw",n=5)),"%PDF-"))
 }
 if(ext=="csv")return(tryCatch(length(names(data.table::fread(path,nrows=0)))>0,error=function(e)FALSE))
 if(ext=="json")return(tryCatch({jsonlite::fromJSON(path);TRUE},error=function(e)FALSE))
 TRUE
}
wf_hash_outputs <- function(paths)setNames(lapply(paths,function(p)digest::digest(file=file.path(RESULT_ROOT,p),algo="sha256")),paths)
wf_parameters <- function(module)list(geneA=GENE_A,
 geneB=if(module %in% c("targeted_dependency","lineage_dependency","adjusted_dependency","cn_threshold_sensitivity","genomewide_dependency"))GENE_B else NULL,
 geneB_provided=if(module=="genomewide_dependency")GENE_B_PROVIDED else NULL,workflow=WORKFLOW,
 bootstrap=if(module=="lineage_dependency")BOOT_R else NULL,min_group_n=MIN_N,low=.585,deep=.35,seed=1234,thresholds=c(.585,.50,.40,.35))
wf_cache_valid <- function(old,key,outputs) {
 if(is.null(old)||!identical(old$key,key)||!identical(names(old$outputs),outputs))return(FALSE)
 if(!all(vapply(file.path(RESULT_ROOT,outputs),wf_validate_artifact,logical(1))))return(FALSE)
 identical(old$outputs,wf_hash_outputs(outputs))
}
wf_run_tcga <- function() {
 path<-wf_table("TCGA_Reference_Samples.csv")
 status<-system2("C:/Python312/python.exe",c(shQuote(file.path(PROJECT_ROOT,"scripts/utils/export_workflow_tcga.py")),
             "--gene",GENE_A,"--output",shQuote(path)))
 if(status!=0)stop("Local processed TCGA extraction failed")
 dat<-fread(path);stopifnot(!anyDuplicated(dat$SampleID),all(dat$TumorNormal=="Tumor"),all(dat$DataLayer=="PanCanAtlas_Xena_reference"))
 order<-cancer_order(dat,TCGA_CANCERS);wf_write(order,"TCGA_Cancer_Order.csv")
 prev<-reference_prevalence(dat,order);wf_write(prev,"TCGA_CNA_Percentage.csv")
 p1<-tcga_prevalence_plot(prev,order);p2<-tcga_landscape_plot(dat,order)
 # ggplot's discrete y after coord_flip is drawn bottom-to-top; reverse levels
 # therefore produce the requested ascending medians from top to bottom.
 stopifnot(identical(levels(p1$data$CancerType),rev(order$CancerType)),
           identical(levels(p2$data$CancerType),levels(p1$data$CancerType)))
 build<-ggplot_build(p2);stopifnot(nrow(build$data[[2]])==sum(is.finite(dat$CopyNumber)))
 workflow_save(p1,file.path(RESULT_ROOT,PATHS["cna"]),width=11,height=11)
 workflow_save(p2,file.path(RESULT_ROOT,PATHS["landscape"]),width=10,height=11)
 stats<-reference_cancer_statistics(dat,order)
 wf_write(stats,paste0("TCGA_",GENE_A,"_CN_mRNA_AllCancer_Statistics.csv"))
 jsonlite::write_json(list(visual_top_to_bottom=order$CancerType,
  CNA_factor_levels=levels(p1$data$CancerType),Landscape_factor_levels=levels(p2$data$CancerType),
  CN_mRNA_statistics_order=stats$CancerType,landscape_points=nrow(build$data[[2]]),
  factor_display_rule="coord_flip: reverse factor levels = table Order top to bottom"),
  wf_provenance("TCGA_Plot_Order.json"),pretty=TRUE,auto_unbox=TRUE)
 cna_counts<-list()
 for(cancer in order$CancerType) {
  s<-stats[CancerType==cancer]
  workflow_save(tcga_expression_plot(dat,cancer,s),file.path(RESULT_ROOT,PATHS["rna"],paste0(cancer,"_",GENE_A,"_CN_mRNA.pdf")),width=9,height=6.5)
  d<-dat[CancerType==cancer & is.finite(CopyNumber)&is.finite(Expression)&is.finite(GISTIC)]
  cna_counts[[cancer]]<-data.table(CancerType=cancer,CNA=CNA_STATES,N=as.integer(table(factor(d$GISTIC,levels=-2:2))))
 }
 wf_write(rbindlist(cna_counts),"TCGA_CN_mRNA_CNA_Counts.csv")
}
wf_run_targeted <- function() {
 pair<-as.data.table(load_target_pair(GENE_A,GENE_B))
 stopifnot(!anyDuplicated(pair$ModelID))
 s<-targeted_statistics(pair);wf_write(s,paste0(GENE_A,"_",GENE_B,"_Targeted_Statistics.csv"))
 wf_write(pair,paste0(GENE_A,"_",GENE_B,"_CellLines.csv"))
 w<-targeted_plots(pair,s,file.path(RESULT_ROOT,PATHS[c("scatter","waterfall","box")]))
 wf_write(w,paste0(GENE_A,"_",GENE_B,"_Waterfall_Order.csv"))
}
wf_run_screen <- function() {
 res<-screen_core_statistics();wf_write(res,"GenomeWide_Dependency.csv")
 wf_write(dependency_top(res),"Top_Dependency_Candidates.csv")
 workflow_save(volcano_plot(res),file.path(RESULT_ROOT,PATHS["volcano"]),width=11,height=7.5)
}
wf_run_covariation <- function() {
 res<-covariation_statistics(fread(FILES$cn),GENE_A);wf_write(res,paste0(GENE_A,"_CN_Covariation.csv"))
 top<-covariation_top(res,GENE_A);wf_write(top,"Top_CN_Covariation.csv")
 workflow_save(covariation_plot(top),file.path(RESULT_ROOT,PATHS["covariation"]),width=12,height=8)
}
wf_run_supplement <- function(module) {
 if(module=="lineage_dependency") {
  pair<-as.data.table(load_target_pair(GENE_A,GENE_B));model<-as.data.table(load_model())
  d<-merge(pair,model[,.(ModelID,OncotreeLineage)],by="ModelID")
  eligibility<-d[!is.na(OncotreeLineage),.(N_low=sum(CN_binary=="CN-Low"),N_nonlow=sum(CN_binary=="CN-NonLow")),by=OncotreeLineage]
  eligibility[,Eligible:=N_low>=MIN_N & N_nonlow>=MIN_N];wf_write(eligibility,"Lineage_Eligibility.csv")
  run_lineage()
  if(!any(eligibility$Eligible)) {
   wf_write(data.table(Lineage=character(),N_low=integer(),N_nonlow=integer(),Delta_median=numeric(),P=numeric(),FDR=numeric()),"Lineage_Dependency.csv")
   workflow_save(ggplot()+annotate("text",0,0,label="Insufficient group N in every lineage")+theme_void(),
                 file.path(out_dir(module),"Lineage_Dependency_Forest.pdf"))
  }
 }
 if(module=="adjusted_dependency") {
  run_adjusted();s<-fread(wf_table("Continuous_CN_Adjusted.csv"));a<-s[term=="CN_log"]
  p<-ggplot(a,aes(estimate,term))+geom_vline(xintercept=0,linetype="dashed")+
   geom_errorbar(aes(xmin=estimate-1.96*std.error,xmax=estimate+1.96*std.error),orientation="y",width=.2)+
   geom_point()+workflow_theme()+labs(title=paste(GENE_B,"Dependency Adjusted for Lineage"),
   x="CN_log coefficient (95% normal interval)",y=NULL,caption="Chronos ~ log2(relative CN + 1) + OncotreeLineage; association only.")
  workflow_save(p,file.path(out_dir(module),"Adjusted_CN_Coefficient.pdf"),height=4)
 }
 if(module=="cn_threshold_sensitivity")run_threshold_sensitivity()
}
wf_summary <- function() {
 old<-options(scipen=0);on.exit(options(old))
 separator<-paste(rep("=",50),collapse="")
 section<-function(title)c("",separator,title,separator)
 lines<-c(separator,"Copy Number–Driven Dependency Analysis",separator,"",paste("Gene A:",GENE_A),
          paste("Gene B:",if(GENE_B_PROVIDED)GENE_B else "未提供（发现性筛选）"),
          "TCGA: TCGA Pan-Cancer","DepMap: Public 26Q1",paste("Run date:",format(Sys.time(),tz="Asia/Shanghai",usetz=TRUE)),
          "DATA_MODE: local（仅读取本地 processed 数据）",
          "TCGA 主图统一使用 reference 层；当前 DR46 可用性与真实来源见 Provenance。",
          "DepMap 为用户提供的 26Q1 导出文件；相对官方完整 release 的完整性未验证。",
          "CN-low: log2(relative CN + 1) < 0.585；分析定义，不是五级 GISTIC。")
 descriptions<-c(cna=paste0("展示 ",GENE_A," 的五级 CNA 比例：Deep Deletion、Shallow Deletion、Diploid、Gain、Amplification；观察缺失比例。"),
 landscape="每个灰点为一个 TCGA tumor sample；箱线图保留全部有限 CN；标签 n 为实际有效样本数；从上到下按 median CN 递增。",
 rna=paste0("逐癌种检验 ",GENE_A," CN 与自身 mRNA 的关系；每张图含 N、Pearson、Spearman、黑色回归线与五级 CNA 百分比。N < 20 不报告相关统计。"),
 scatter=paste0(GENE_A," relative CN 与 ",GENE_B," dependency 的相关性；右上相关统计使用 relative CN。表中另保留历史 CN_log 统计。"),
 waterfall="匹配细胞按 Chronos 从高到低排序，即弱依赖到强依赖；黑色为 CN-Normal，红色为 CN-Low；参考线为 -1。",
 box="比较 CN-Normal 与 CN-Low 的 Chronos；展示全部点、各组 N、median、Delta median、Wilcoxon P 与显著性星号。",
 volcano="完整本地 Chronos genome 筛选；Chronos 越负依赖越强。Delta median < 0 表示 CN-low 组依赖更强；用户排名为 Eligible_Rank。图标记每侧 Top10；候选表保留每侧 Top20。",
 covariation="CN 共变图展示正相关 Top20 与负相关 Top20；正相关表示 CN 往往同向变化，负相关表示往往反向变化；这是相关性，不是因果关系。")
 for(key in names(PATHS))lines<-c(lines,section(paste(basename(PATHS[key]))),paste("文件/目录：",PATHS[key]),"数据：",if(key %in% c("cna","landscape","rna"))"TCGA" else "DepMap","用途/怎么看：",descriptions[key])
 if(GENE_B_PROVIDED)lines<-c(lines,section("Supplementary"),
  "01_Lineage：各 OncotreeLineage 内 CN-low vs nonlow 的 dependency；小组不够 N 时明确不检验。",
  "02_Adjusted：Chronos ~ CN_log + Lineage，并保留二元 CN-low 调整表。",
  "03_CN_Threshold_Sensitivity：固定预设 cutoff 的稳健性分析；不为预期结果调整阈值。")
 prev<-fread(wf_table("TCGA_CNA_Percentage.csv"));loss<-prev[CNA %in% CNA_STATES[1:2],.(Loss_percentage=sum(Percentage)),by=CancerType];setorder(loss,-Loss_percentage)
 stats<-fread(wf_table(paste0("TCGA_",GENE_A,"_CN_mRNA_AllCancer_Statistics.csv")))
 findings<-c(paste0("- ",GENE_A," deletion 比例最高的癌种：",paste(paste0(head(loss$CancerType,3)," ",round(head(loss$Loss_percentage,3),1),"%"),collapse="；")),
 paste0("- CN-mRNA：",nrow(stats)," 癌种图；",sum(stats$Pearson_FDR<.05,na.rm=TRUE)," 癌种 Pearson BH-FDR < 0.05。"))
 if(GENE_B_PROVIDED) {
  target<-fread(wf_table(paste0(GENE_A,"_",GENE_B,"_Targeted_Statistics.csv")))
  screen<-fread(wf_table("GenomeWide_Dependency.csv"))[Gene==GENE_B]
  findings<-c(findings,paste0("- ",GENE_B," Eligible_Rank：",if(nrow(screen))screen$Eligible_Rank else "不在本地 Chronos 导出中"),
   paste0("- CN-low N=",target$N_low,"；CN-normal N=",target$N_nonlow,"。"),
   paste0("- median Chronos：low ",signif(target$Median_low,5),"；normal ",signif(target$Median_nonlow,5),"；Delta ",signif(target$Delta_median,5),"。"),
   paste0("- Wilcoxon P ",p_text(target$Wilcoxon_P),"；CN_log Pearson r ",signif(target$Pearson_r,4),"；Spearman rho ",signif(target$Spearman_rho,4),"。"))
  lin<-fread(wf_table("Lineage_Dependency.csv"))
  supported<-lin[is.finite(FDR)&FDR<.05&Delta_median<0,Lineage]
  findings<-c(findings,paste0("- 同方向且 FDR < 0.05 的 lineage：",if(length(supported))paste(supported,collapse="、") else "无；详见各 lineage 的样本量和效应。"))
 }
 findings<-c(findings,"- 以上为观察性相关与依赖差异，不能解释为 CN 缺失导致 dependency 或其他基因 CN 改变。")
 lines<-c(lines,section("Key Findings"),head(findings,15),section("File Index"))
 files<-sort(list.files(RESULT_ROOT,recursive=TRUE))
 files<-setdiff(files,c("00_Analysis_Summary.txt","Provenance/File_Index.csv"))
 description<-function(f) {
  key<-names(PATHS)[match(f,PATHS)]
  if(!is.na(key))return(unname(descriptions[key]))
  if(startsWith(f,paste0(PATHS["rna"],"/")))return(paste(GENE_A,"在该癌种 CN 与自身 mRNA 的散点和统计"))
  if(startsWith(f,"Tables/"))return("统计数值、样本匹配或绘图输入；用于复核与重新制图")
  if(startsWith(f,"Supplementary/"))return("lineage、调整模型或固定 CN cutoff 的补充结果")
  "数据真实来源、版本、哈希、运行记录或验证信息"
 }
 for(i in seq_along(files))lines<-c(lines,paste0(sprintf("%02d",i)," ",files[i]),paste("作用：",description(files[i])))
 # Index includes itself and the summary; exact on-disk correspondence is validated.
 lines<-c(lines,"00_Analysis_Summary.txt","作用：中文结果阅读指南、关键统计与文件索引。",
          "Provenance/File_Index.csv","作用：本目录完整文件清单（包括本清单与 Summary）。")
 writeLines(enc2utf8(lines),file.path(RESULT_ROOT,"00_Analysis_Summary.txt"),useBytes=TRUE)
 files<-sort(unique(c(files,"00_Analysis_Summary.txt","Provenance/File_Index.csv")))
 data.table::fwrite(data.table(Path=files),wf_provenance("File_Index.csv"))
}
run_final_workflow <- function() {
 PATHS<<-workflow_paths(WORKFLOW,GENE_A,GENE_B)
 TCGA_CANCERS<<-sort(unique(unname(unlist(jsonlite::fromJSON(file.path(PROJECT_ROOT,"config/tcga_cancer_types.json"))$mapping))))
 for(p in c("Main_Results","Tables","Provenance"))dir.create(file.path(RESULT_ROOT,p),recursive=TRUE,showWarnings=FALSE)
 OUTPUT_DIRS[c("lineage_dependency","adjusted_dependency","cn_threshold_sensitivity")]<<-
  c("Supplementary/01_Lineage","Supplementary/02_Adjusted","Supplementary/03_CN_Threshold_Sensitivity")
 wf_route_supplements()
 if(GENE_B_PROVIDED)for(p in OUTPUT_DIRS[c("lineage_dependency","adjusted_dependency","cn_threshold_sensitivity")])dir.create(file.path(RESULT_ROOT,p),recursive=TRUE,showWarnings=FALSE)
 # Legacy supplement functions write their tables through this scoped adapter.
 fwrite<<-function(x,file,...) {
  if(grepl("/Supplementary/",gsub("\\\\","/",file),fixed=TRUE)&&grepl("\\.csv$",file))file<-wf_table(basename(file))
  data.table::fwrite(x,file,...,scipen=0)
 }
 fun<-list(tcga=wf_run_tcga,targeted_dependency=wf_run_targeted,genomewide_dependency=wf_run_screen,cn_covariation=wf_run_covariation)
 runs<-list();inputs<-list()
 for(module in workflow_modules(WORKFLOW)) {
  start<-Sys.time();outputs<-unname(wf_expected(module));versions<-wf_input_versions(wf_input_paths(module));inputs[[module]]<-versions
  code<-wf_code_signature(module)
  params<-wf_parameters(module)
  key<-digest::digest(list(inputs=versions,params=params,code=code$hash),algo="sha256")
  cache_path<-wf_provenance(paste0("Cache_",module,".json"))
  old<-if(file.exists(cache_path))tryCatch(jsonlite::fromJSON(cache_path,simplifyVector=FALSE),error=function(e)NULL) else NULL
  state<-"computed"
  if(!opt$force&&wf_cache_valid(old,key,outputs)) {state<-"cached";msg(module," SKIP recomputation: verified cache")}
  else {
   msg(module," running once")
   if(module %in% names(fun))fun[[module]]() else wf_run_supplement(module)
   if(!all(vapply(file.path(RESULT_ROOT,outputs),wf_validate_artifact,logical(1))))stop("Output validation failed: ",module)
   jsonlite::write_json(list(module=module,key=key,parameters=params,parameter_hash=digest::digest(params,algo="sha256"),
    code_hash=code$hash,code_functions=code$functions,data_versions=versions,outputs=wf_hash_outputs(outputs),validated=TRUE),
    cache_path,pretty=TRUE,auto_unbox=TRUE)
  }
  runs[[module]]<-data.table(Module=module,Status=state,Seconds=as.numeric(difftime(Sys.time(),start,units="secs")))
  gc()
 }
 data.table::fwrite(rbindlist(runs),wf_provenance("Module_Runs.csv"))
 metadata<-list(geneA=GENE_A,geneB=if(GENE_B_PROVIDED)GENE_B else NULL,workflow=WORKFLOW,DATA_MODE=DATA_MODE,
  R_version=as.character(getRversion()),run_date=format(Sys.time(),tz="Asia/Shanghai",usetz=TRUE),
  TCGA_main_layer="PanCanAtlas_Xena_reference",TCGA_current_source="NCI GDC DR46; not used for main figures",
  DepMap_release="26Q1",DepMap_full_release_completeness="unverified; all supplied Chronos gene columns tested",
  input_versions=inputs,parameters=list(low=.585,deep=.35,min_group_n=MIN_N,bootstrap=BOOT_R,seed=1234),
  code_files=setNames(lapply(list.files("scripts/R",pattern="\\.R$",full.names=TRUE),function(f)digest::digest(file=f,algo="sha256")),list.files("scripts/R",pattern="\\.R$",full.names=TRUE)))
 jsonlite::write_json(metadata,wf_provenance("Run_Metadata.json"),pretty=TRUE,auto_unbox=TRUE)
 writeLines(trimws(capture.output(sessionInfo()),which="right"),wf_provenance("SessionInfo.txt"))
 jsonlite::write_json(inputs,wf_provenance("Input_Versions.json"),pretty=TRUE,auto_unbox=TRUE)
 writeLines(c("TCGA current source: NCI GDC", "Release: DR46 / 46.0", "Current processed data retained; not used in default main figures.",
  "TCGA main source: PanCanAtlas reference / UCSC Xena", "Continuous CN: GISTIC2 all_data_by_genes; source scale, not absolute CN.",
  "Five states: thresholded GISTIC2 (-2,-1,0,1,2).", "Expression: matched TCGA Toil RSEM, log2(norm_count+1); not TPM.",
  "Reference joins: exact SampleID; unique IDs in each layer; TumorNormal=Tumor; finite CN/mRNA/GISTIC for scatter correlations.",
  "Landscape/order: all finite reference CN tumor samples; no outlier exclusion.",
  "Source paths, original URLs, source hashes and processor versions: Input_Versions.json.",
  "DepMap: Public 26Q1, user-supplied portal exports; full-release completeness unverified.",
  "Dependency grouping and continuous core: CN_log=log2(relative CN+1), low<0.585.",
  "Targeted scatter statistics: relative CN; historical CN_log correlation retained in table.",
  "Screen delta: median(low)-median(nonlow); BH over the complete supplied screen, including untested NA rows as in supplied code.",
  "Eligible_Rank excludes undefined FDR; original Rank retained for provenance.",
  "All associations are observational; no causal inference."),wf_provenance("Source_Metadata.txt"))
 # A small readable version file accompanies machine-readable source hashes.
 writeLines(c("DATA_MODE=local","TCGA reference=PanCanAtlas/Xena (existing snapshot)","TCGA current=NCI GDC DR46 (retained separately)","DepMap=Public 26Q1",paste("R=",getRversion()),"Input SHA256 and code hashes: Input_Versions.json / Run_Metadata.json / Cache_*.json"),wf_provenance("Input_Versions.txt"))
 wf_summary()
 source(file.path(PROJECT_ROOT,"scripts/R/validate_final_workflow.R"),encoding="UTF-8")
 validate_final_workflow()
 wf_summary() # Include the newly saved validation report in the exact file index.
 msg("Analysis completed: ",RESULT_ROOT)
}

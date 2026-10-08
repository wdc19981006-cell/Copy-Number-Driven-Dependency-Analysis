# Orchestration, module caches and the final user-facing result contract.
workflow_paths <- function(workflow,geneA,geneB=NULL) {
 pair<-identical(workflow,"geneA_geneB")
 c(cna=paste0("Main_Results/01_TCGA_CNA_Percentage/01_TCGA_",geneA,"_CNA_Percentage.pdf"),
 landscape=paste0("Main_Results/02_TCGA_CopyNumber_Landscape/02_TCGA_",geneA,"_CopyNumber_Landscape.pdf"),
 rna="Main_Results/03_TCGA_CN_mRNA",
 if(pair)c(scatter=paste0("Main_Results/04_DepMap_CN_vs_Dependency/04_DepMap_",geneA,"_CN_vs_",geneB,"_Dependency.pdf"),
 waterfall=paste0("Main_Results/05_DepMap_Dependency_Waterfall/05_DepMap_",geneB,"_Dependency_Waterfall.pdf"),
 box=paste0("Main_Results/06_DepMap_Dependency_CNlow_vs_Normal/06_DepMap_",geneB,"_Dependency_CNlow_vs_Normal.pdf")),
 volcano=paste0("Main_Results/",if(pair)"07" else "04","_DepMap_GenomeWide_Dependency/",if(pair)"07" else "04","_DepMap_",geneA,"_GenomeWide_Dependency_Volcano.pdf"),
 covariation=paste0("Main_Results/",if(pair)"08" else "05","_DepMap_CN_Covariation/",if(pair)"08" else "05","_DepMap_",geneA,"_CN_Covariation.pdf"))
}
wf_table_relative <- function(name) {
 key<-switch(name,TCGA_CNA_Percentage.csv="cna",TCGA_Cancer_Order.csv="landscape",
  TCGA_CN_mRNA_CNA_Counts.csv="rna",GenomeWide_Dependency.csv="volcano",Top_Dependency_Candidates.csv="volcano",
  Top_CN_Covariation.csv="covariation",NULL)
 if(name %in% c("TCGA_Current_Samples.csv","TCGA_RNA_Representative_Selection.csv","TCGA_Sample_Baselines.csv"))
  return(file.path("Provenance/Data_Audit",name))
 if(name==paste0("TCGA_",GENE_A,"_CN_mRNA_AllCancer_Statistics.csv"))key<-"rna"
 if(name==paste0(GENE_A,"_CN_Covariation.csv"))key<-"covariation"
 pair<-paste0(GENE_A,"_",GENE_B)
 if(name %in% paste0(pair,c("_Targeted_Statistics.csv","_CellLines.csv")))key<-"scatter"
 if(name==paste0(pair,"_Waterfall_Order.csv"))key<-"waterfall"
 if(name==paste0(pair,"_CNlow_vs_Normal_Statistics.csv"))key<-"box"
 supplements<-c(Lineage_Dependency.csv="01_Lineage",Lineage_Eligibility.csv="01_Lineage",
  Continuous_CN_Adjusted.csv="02_Adjusted",CNLow_Adjusted.csv="02_Adjusted",
  CN_Threshold_Sensitivity.csv="03_CN_Threshold_Sensitivity",CN_Threshold_Continuous_Statistics.csv="03_CN_Threshold_Sensitivity")
 if(name %in% names(supplements))return(file.path("Supplementary",supplements[[name]],name))
 if(is.null(key)||!key %in% names(PATHS))stop("No output module assigned to table: ",name)
 file.path(if(key=="rna")PATHS[[key]] else dirname(PATHS[[key]]),name)
}
wf_table <- function(name)file.path(RESULT_ROOT,wf_table_relative(name))
wf_provenance <- function(name)file.path(RESULT_ROOT,"Provenance",name)
wf_write <- function(x,name) data.table::fwrite(x,wf_table(name),scipen=0)
wf_create_directories <- function() {
 folders<-c("Provenance/Data_Audit",PATHS["rna"],dirname(PATHS[names(PATHS)!="rna"]))
 if(GENE_B_PROVIDED)folders<-c(folders,"Supplementary/01_Lineage","Supplementary/02_Adjusted","Supplementary/03_CN_Threshold_Sensitivity")
 for(p in unique(folders))dir.create(file.path(RESULT_ROOT,p),recursive=TRUE,showWarnings=FALSE)
}
wf_group_statistics <- function(s) {
 # Project the already computed targeted row; never run another statistical test.
 fields<-c("geneA","geneB","N_low","N_nonlow","Median_low","Median_nonlow","Delta_median","Wilcoxon_P")
 s[,intersect(fields,names(s)),with=FALSE]
}
wf_screen_status_path <- function()file.path(dirname(PATHS[["volcano"]]),paste0(if(GENE_B_PROVIDED)"07" else "04","_Analysis_Status.txt"))
wf_write_screen_status <- function(status) {
 lines<-c(paste0("分析状态：",status$Status),paste0("CN-Low n=",status$N_low,"；CN-Normal n=",status$N_nonlow),
  paste0("当前最低每组样本数要求 n>=",status$min_group_n))
 lines<-c(lines,if(status$Status=="SKIPPED")c("低CN组或正常CN组样本量不足。","没有可靠的组间dependency筛选结果。","PDF为状态说明；结果CSV为空，不提供候选基因。") else
  "筛选已完成；效应、P、FDR与Eligible_Rank见本模块CSV。")
 writeLines(enc2utf8(lines),file.path(RESULT_ROOT,wf_screen_status_path()),useBytes=TRUE)
}
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
 ref<-file.path(PROJECT_ROOT,"data/processed/tcga/gdc_DR46")
 switch(module,tcga=file.path(ref,c("TCGA_GeneLevel_CN.parquet","TCGA_GeneLevel_CN.genes.parquet",
                                   "TCGA_STAR_TPM.parquet","TCGA_STAR_TPM.genes.parquet")),
 genomewide_dependency=c(FILES$cn,FILES$chronos),cn_covariation=FILES$cn,
 c(FILES$cn,FILES$chronos,if(module %in% c("lineage_dependency","adjusted_dependency"))FILES$model))
}
wf_input_versions <- function(paths) {
 lapply(paths,function(p) {
  local_required(p);s<-file.info(p);side<-paste0(p,".provenance.json")
  if(!file.exists(side))side<-sub("\\.parquet$",".provenance.json",p)
  meta<-if(file.exists(side))jsonlite::fromJSON(side,simplifyVector=FALSE) else NULL
  list(path=substring(p,nchar(PROJECT_ROOT)+2),bytes=s$size,mtime=format(s$mtime,"%Y-%m-%dT%H:%M:%OS6%z"),
       prepared_sha256=if(!is.null(meta$SHA256))meta$SHA256 else meta$output_sha256,
       provenance=meta)
 })
}
wf_code_signature <- function(module) {
 functions<-switch(module,
 tcga=c("tcga_relative_cn","cancer_order","ordered_cancer","current_cancer_statistics","current_prevalence",
        "tcga_prevalence_plot","tcga_landscape_plot","tcga_expression_plot","tcga_cna_pie","wf_run_tcga","pair_correlations"),
 targeted_dependency=c("targeted_statistics","targeted_plots","targeted_scatter_plot","waterfall_order","wf_run_targeted","wf_targeted_statistics_key","pair_correlations","group_effect"),
 genomewide_dependency=c("run_genomewide","screen_core_statistics","screen_group_counts","empty_screen_statistics",
                        "eligible_rank","dependency_top","volcano_plot","wf_run_screen"),
 cn_covariation=c("covariation_statistics","covariation_top","covariation_plot","wf_run_covariation"),
 lineage_dependency=c("run_lineage","wf_run_supplement"),adjusted_dependency=c("run_adjusted","wf_run_supplement","coefficient_table"),
 cn_threshold_sensitivity=c("run_threshold_sensitivity","cn_threshold_statistics","group_effect","wf_run_supplement"))
 common<-c("prepare_cn","load_target_pair","read_gene","load_model","fread","find_gene_col","clean_gene",
           "processed_depmap","processed_pair","wf_write","workflow_save","workflow_theme","save_pdf",
           "workflow_paths","wf_table","wf_table_relative","wf_group_statistics","wf_write_screen_status","wf_screen_status_path","wf_validate_artifact","p_text","p_star")
 # Hash source text rather than serialized language objects: runtime/JIT or
 # package expression caches can change attributes of an otherwise identical
 # function body after plotting. Text stays stable across cold and warm calls.
 code<-lapply(unique(c(functions,common)),function(n)list(name=n,
  args=paste(deparse(formals(get(n)),width.cutoff=500L),collapse="\n"),
  body=paste(deparse(body(get(n)),width.cutoff=500L),collapse="\n")))
 py<-"scripts/R/local_data.R"
 if(module=="tcga")py<-c(py,"config/workflows.json","scripts/utils/export_workflow_tcga.py","scripts/utils/tcga_current.py","scripts/data_access.py","config/tcga_cancer_types.json")
 list(hash=digest::digest(list(code=code,files=lapply(py,function(p)digest::digest(file=file.path(PROJECT_ROOT,p),algo="sha256")),
               constants=list(states=CNA_STATES,colors=CNA_COLORS,tcga_min_n=TCGA_MIN_N),
               expected=unname(wf_expected(module))),algo="sha256"),functions=unique(c(functions,common)))
}
wf_expected <- function(module) {
 tables<-function(n)vapply(n,wf_table_relative,character(1),USE.NAMES=FALSE)
 switch(module,
 tcga=c(PATHS[c("cna","landscape")],paste0(PATHS["rna"],"/",TCGA_CANCERS,"_",GENE_A,"_CN_mRNA.pdf"),
        tables(c("TCGA_Cancer_Order.csv","TCGA_CNA_Percentage.csv","TCGA_Current_Samples.csv",
        "TCGA_Sample_Baselines.csv","TCGA_RNA_Representative_Selection.csv")),"Provenance/TCGA_Baseline_Method.json",
        tables(paste0("TCGA_",GENE_A,"_CN_mRNA_AllCancer_Statistics.csv")),"Provenance/TCGA_Plot_Order.json",tables("TCGA_CN_mRNA_CNA_Counts.csv")),
 targeted_dependency=c(PATHS[c("scatter","waterfall","box")],tables(paste0(GENE_A,"_",GENE_B,
          c("_Targeted_Statistics.csv","_CellLines.csv","_Waterfall_Order.csv","_CNlow_vs_Normal_Statistics.csv")))),
 genomewide_dependency=c(PATHS["volcano"],tables(c("GenomeWide_Dependency.csv","Top_Dependency_Candidates.csv")),
                        "Provenance/GenomeWide_Dependency_Status.json",wf_screen_status_path()),
 cn_covariation=c(PATHS["covariation"],tables(c(paste0(GENE_A,"_CN_Covariation.csv"),"Top_CN_Covariation.csv"))),
 lineage_dependency=c("Supplementary/01_Lineage/Lineage_Dependency_Forest.pdf",tables(c("Lineage_Dependency.csv","Lineage_Eligibility.csv"))),
 adjusted_dependency=c("Supplementary/02_Adjusted/Adjusted_CN_Coefficient.pdf",tables(c("Continuous_CN_Adjusted.csv","CNLow_Adjusted.csv"))),
 cn_threshold_sensitivity=c("Supplementary/03_CN_Threshold_Sensitivity/CN_Threshold_Sensitivity.pdf",
                            tables(c("CN_Threshold_Sensitivity.csv","CN_Threshold_Continuous_Statistics.csv"))))
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
 path<-wf_table("TCGA_Current_Samples.csv")
 status<-system2("C:/Python312/python.exe",c(shQuote(file.path(PROJECT_ROOT,"scripts/utils/export_workflow_tcga.py")),
             "--gene",GENE_A,"--output",shQuote(path)))
 if(status!=0)stop("Local processed TCGA extraction failed")
 dat<-fread(path);stopifnot(!anyDuplicated(dat$SampleID),all(dat$TumorNormal=="Tumor"),all(dat$DataLayer=="gdc_current_DR46"))
 dat<-tcga_relative_cn(dat);wf_write(dat,"TCGA_Current_Samples.csv")
 order<-cancer_order(dat,TCGA_CANCERS);wf_write(order,"TCGA_Cancer_Order.csv")
 prev<-current_prevalence(dat,order);wf_write(prev,"TCGA_CNA_Percentage.csv")
 p1<-tcga_prevalence_plot(prev,order);p2<-tcga_landscape_plot(dat,order)
 # ggplot's discrete y after coord_flip is drawn bottom-to-top; reverse levels
 # therefore produce the requested ascending medians from top to bottom.
 stopifnot(identical(levels(p1$data$CancerType),rev(order$CancerType)),
           identical(levels(p2$data$CancerType),levels(p1$data$CancerType)))
 build<-ggplot_build(p2);stopifnot(nrow(build$data[[2]])==sum(is.finite(dat$CopyNumber)))
 workflow_save(p1,file.path(RESULT_ROOT,PATHS["cna"]),width=11,height=11)
 workflow_save(p2,file.path(RESULT_ROOT,PATHS["landscape"]),width=10,height=11)
 stats<-current_cancer_statistics(dat,order)
 wf_write(stats,paste0("TCGA_",GENE_A,"_CN_mRNA_AllCancer_Statistics.csv"))
 jsonlite::write_json(list(visual_top_to_bottom=order$CancerType,
  CNA_factor_levels=levels(p1$data$CancerType),Landscape_factor_levels=levels(p2$data$CancerType),
  CN_mRNA_statistics_order=stats$CancerType,landscape_points=nrow(build$data[[2]]),
  CN_metric="Relative_CN_Change",formula="CopyNumber / BaselineCN - 1",
  landscape_zero_reference=build$data[[3]]$yintercept,
  factor_display_rule="coord_flip: reverse factor levels = table Order top to bottom"),
  wf_provenance("TCGA_Plot_Order.json"),pretty=TRUE,auto_unbox=TRUE)
 cna_counts<-list()
 for(cancer in order$CancerType) {
  s<-stats[CancerType==cancer]
  workflow_save(tcga_expression_plot(dat,cancer,s),file.path(RESULT_ROOT,PATHS["rna"],paste0(cancer,"_",GENE_A,"_CN_mRNA.pdf")),width=9,height=6.5)
  d<-dat[CancerType==cancer & is.finite(CopyNumber)&is.finite(Expression)&is.finite(CNAState)]
  cna_counts[[cancer]]<-data.table(CancerType=cancer,CNA=CNA_STATES,N=as.integer(table(factor(d$CNAState,levels=-2:2))))
 }
 wf_write(rbindlist(cna_counts),"TCGA_CN_mRNA_CNA_Counts.csv")
}
wf_targeted_statistics_key <- function() {
 functions<-c("targeted_statistics","pair_correlations","group_effect","prepare_cn","load_target_pair",
              "read_gene","fread","find_gene_col","clean_gene","targeted_plots","waterfall_order")
 digest::digest(list(inputs=wf_input_versions(wf_input_paths("targeted_dependency")),params=wf_parameters("targeted_dependency"),
  code=lapply(functions,function(n)list(formals(get(n)),paste(deparse(body(get(n)),width.cutoff=500L),collapse="\n")))),algo="sha256")
}
wf_run_targeted <- function() {
 cache<-wf_provenance("Cache_targeted_dependency.json")
 old<-if(file.exists(cache))jsonlite::fromJSON(cache,simplifyVector=FALSE) else NULL
 preserved<-setdiff(unname(wf_expected("targeted_dependency")),unname(PATHS["scatter"]))
 if(!opt$force&&!is.null(old$statistics_key)&&identical(old$statistics_key,wf_targeted_statistics_key())&&
    all(vapply(file.path(RESULT_ROOT,preserved),wf_validate_artifact,logical(1)))&&
    identical(old$outputs[preserved],wf_hash_outputs(preserved))) {
  pair<-fread(wf_table(paste0(GENE_A,"_",GENE_B,"_CellLines.csv")))
  s<-fread(wf_table(paste0(GENE_A,"_",GENE_B,"_Targeted_Statistics.csv")))
  wf_write(wf_group_statistics(s),paste0(GENE_A,"_",GENE_B,"_CNlow_vs_Normal_Statistics.csv"))
  workflow_save(targeted_scatter_plot(pair,s),file.path(RESULT_ROOT,PATHS["scatter"]))
  return("plot_regenerated")
 }
 pair<-as.data.table(load_target_pair(GENE_A,GENE_B))
 stopifnot(!anyDuplicated(pair$ModelID))
 s<-targeted_statistics(pair);wf_write(s,paste0(GENE_A,"_",GENE_B,"_Targeted_Statistics.csv"))
 wf_write(wf_group_statistics(s),paste0(GENE_A,"_",GENE_B,"_CNlow_vs_Normal_Statistics.csv"))
 wf_write(pair,paste0(GENE_A,"_",GENE_B,"_CellLines.csv"))
 w<-targeted_plots(pair,s,file.path(RESULT_ROOT,PATHS[c("scatter","waterfall","box")]))
 wf_write(w,paste0(GENE_A,"_",GENE_B,"_Waterfall_Order.csv"))
 invisible(NULL)
}
wf_run_screen <- function() {
 counts<-screen_group_counts()
 skipped<-counts$N_low<MIN_N||counts$N_nonlow<MIN_N
 detail<-if(skipped)paste0("Insufficient CN group N: CN-Low n=",counts$N_low,
  "; CN-Normal n=",counts$N_nonlow,"; each group requires n >= ",MIN_N,". CN_log threshold remains 0.585.") else ""
 status<-list(Status=if(skipped)"SKIPPED" else "COMPLETED",Reason=detail,
  N_low=counts$N_low,N_nonlow=counts$N_nonlow,min_group_n=MIN_N,CN_log_threshold=.585)
 jsonlite::write_json(status,wf_provenance("GenomeWide_Dependency_Status.json"),pretty=TRUE,auto_unbox=TRUE)
 wf_write_screen_status(status)
 if(skipped) {
  res<-empty_screen_statistics();wf_write(res,"GenomeWide_Dependency.csv")
  wf_write(dependency_top(res),"Top_Dependency_Candidates.csv")
  p<-ggplot()+annotate("text",0,0,label=paste0("SKIPPED: insufficient CN group N\nCN-Low n = ",counts$N_low,
    "; CN-Normal n = ",counts$N_nonlow,"\nEach group requires n >= ",MIN_N),size=5)+theme_void()+
    labs(title=paste(GENE_A,"CN-Low Genome-wide Dependency Screen"),
         caption="CN_log = log2(DepMap relative CN + 1); CN-Low < 0.585. No group comparison or candidates reported.")
  workflow_save(p,file.path(RESULT_ROOT,PATHS["volcano"]),width=11,height=7.5)
  return(list(status="SKIPPED",detail=detail))
 }
 res<-screen_core_statistics();wf_write(res,"GenomeWide_Dependency.csv")
 wf_write(dependency_top(res),"Top_Dependency_Candidates.csv")
 workflow_save(volcano_plot(res),file.path(RESULT_ROOT,PATHS["volcano"]),width=11,height=7.5)
 invisible(NULL)
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
          "TCGA: Current TCGA dataset","DepMap: Public 26Q1",paste("Run date:",if(file.exists(wf_provenance("Run_Metadata.json")))jsonlite::fromJSON(wf_provenance("Run_Metadata.json"))$run_date else format(Sys.time(),tz="Asia/Shanghai",usetz=TRUE)),
          "DATA_MODE: local（仅读取本地 processed 数据）",
          "TCGA 主结果使用本地已验证的 current TCGA 数据；Source / release 信息见 Provenance。",
          "五级 CNA 基于 current gene-level CN 与 sample-specific baseline/ploidy，是 analysis-defined 分类，不是官方 GISTIC 五级值。",
          "baseline 使用每个样本常染色体 gene-level CN 的整数众数估算；不是直接测量的 ploidy。",
          "DepMap 为用户提供的 26Q1 导出文件；相对官方完整 release 的完整性未验证。",
          "CN-low: log2(relative CN + 1) < 0.585；分析定义，不是五级 GISTIC。")
 descriptions<-c(cna=paste0("展示 ",GENE_A," 的五级 CNA 比例：Deep Deletion、Shallow Deletion、Diploid、Gain、Amplification；观察缺失比例。"),
 landscape="本图展示 Gene A 相对于每个肿瘤样本自身 CN baseline 的变化。Relative CN Change = Gene CN / Sample Baseline CN - 1；0：相对未改变；<0：相对拷贝数减少；>0：相对拷贝数增加。每个灰点为一个 TCGA tumor sample；箱线图保留全部有限 CN 样本；标签 n 为实际有效样本数；01/02 从上到下按 median(Relative CN Change) 递增。这是 analysis-derived relative copy-number change，不是 GISTIC、log2 GISTIC 或 DepMap CN_relative。",
 rna=paste0("逐癌种检验 ",GENE_A," CN 与自身 mRNA 的关系；横轴及 Pearson/Spearman 统一使用 Relative CN Change = CN / BaselineCN - 1；Expression = log2(TPM + 1)。每张图含 N、Pearson、Spearman、黑色回归线与匹配 cohort 的五级 CNA pie chart。N < 20 不报告相关统计。"),
 scatter=paste0(GENE_A," relative CN 与 ",GENE_B," dependency 的相关性；右上相关统计使用 relative CN。表中另保留历史 CN_log 统计。"),
 waterfall="匹配细胞按 Chronos 从高到低排序，即弱依赖到强依赖；黑色为 CN-Normal，红色为 CN-Low；参考线为 -1。",
 box="比较 CN-Normal 与 CN-Low 的 Chronos；展示全部点、各组 N、median、Delta median、Wilcoxon P 与显著性星号。",
 volcano="完整本地 Chronos genome 筛选；Chronos 越负依赖越强。Delta median < 0 表示 CN-low 组依赖更强；用户排名为 Eligible_Rank。图标记每侧 Top10；候选表保留每侧 Top20。",
 covariation="CN 共变图展示正相关 Top20 与负相关 Top20。Pearson r 表示在 DepMap 细胞系中 Gene A 与另一个基因 CN 的线性相关程度；r > 0 表示往往同向变化，r < 0 表示往往反向变化；|r| 越接近 1 相关越强，r 接近 0 线性相关越弱。相关不代表因果。")
 screen_status<-jsonlite::fromJSON(wf_provenance("GenomeWide_Dependency_Status.json"))
 if(screen_status$Status=="SKIPPED")descriptions["volcano"]<-paste0("SKIPPED：CN 分组样本不足，未进行组间 dependency 筛选；CN-Low n=",
  screen_status$N_low,"，CN-Normal n=",screen_status$N_nonlow,"，每组至少 ",MIN_N,"。阈值仍为 0.585；此 PDF 为跳过说明，候选表为空。")
 for(key in names(PATHS)) {
  folder<-if(key=="rna")PATHS[[key]] else dirname(PATHS[[key]])
  csv<-sort(list.files(file.path(RESULT_ROOT,folder),pattern="\\.csv$"))
  lines<-c(lines,section(basename(folder)),"分析目的：",descriptions[key],
   paste("数据库：",if(key %in% c("cna","landscape","rna"))"TCGA（本地 current DR46）" else "DepMap（本地 26Q1）"),
   paste0("PDF位置： ",PATHS[key],if(key=="rna")" （33癌种，各一个PDF）" else ""),
   "统计CSV位置：",file.path(folder,csv),"如何解释：",descriptions[key])
 }
 if(GENE_B_PROVIDED)lines<-c(lines,section("Supplementary"),
  "01_Lineage：各 OncotreeLineage 内 CN-low vs nonlow 的 dependency；小组不够 N 时明确不检验。",
  "02_Adjusted：Chronos ~ CN_log + Lineage，并保留二元 CN-low 调整表。",
  "03_CN_Threshold_Sensitivity：固定预设 cutoff 的稳健性分析；不为预期结果调整阈值。")
 prev<-fread(wf_table("TCGA_CNA_Percentage.csv"));loss<-prev[CNA %in% CNA_STATES[1:2],.(Loss_percentage=sum(Percentage)),by=CancerType];setorder(loss,-Loss_percentage)
 stats<-fread(wf_table(paste0("TCGA_",GENE_A,"_CN_mRNA_AllCancer_Statistics.csv")))
 findings<-c(paste0("- ",GENE_A," deletion 比例最高的癌种：",paste(paste0(head(loss$CancerType,3)," ",round(head(loss$Loss_percentage,3),1),"%"),collapse="；")),
 paste0("- CN-mRNA：",nrow(stats)," 癌种图；",sum(stats$Pearson_FDR<.05,na.rm=TRUE)," 癌种 Pearson BH-FDR < 0.05。"))
 if(screen_status$Status=="SKIPPED")findings<-c(findings,paste0("- Genome-wide dependency SKIPPED：CN-Low n=",
  screen_status$N_low,"；CN-Normal n=",screen_status$N_nonlow,"；样本量不足，不能据此解释为没有依赖关联。"))
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
 # Keep the exact path index, without presenting technical files as analyses.
 lines<-c(lines,files[!startsWith(files,"Provenance/")],"00_Analysis_Summary.txt",
  "技术审计、版本和验证信息保存在Provenance",files[startsWith(files,"Provenance/")])
 # Index includes itself and the summary; exact on-disk correspondence is validated.
 lines<-c(lines,"Provenance/File_Index.csv")
 writeLines(enc2utf8(lines),file.path(RESULT_ROOT,"00_Analysis_Summary.txt"),useBytes=TRUE)
 files<-sort(unique(c(files,"00_Analysis_Summary.txt","Provenance/File_Index.csv")))
 data.table::fwrite(data.table(Path=files),wf_provenance("File_Index.csv"))
}
run_final_workflow <- function() {
 PATHS<<-workflow_paths(WORKFLOW,GENE_A,GENE_B)
 TCGA_CANCERS<<-sort(unique(unname(unlist(jsonlite::fromJSON(file.path(PROJECT_ROOT,"config/tcga_cancer_types.json"))$mapping))))
 wf_create_directories()
 if(dir.exists(file.path(RESULT_ROOT,"Tables")))stop("Legacy Tables layout: migrate existing results before running this workflow")
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
  state<-"computed";execution<-"computed";detail<-""
  if(!opt$force&&wf_cache_valid(old,key,outputs)) {
   state<-if(identical(old$analysis_status,"SKIPPED"))"SKIPPED" else "cached"
   execution<-"cached";detail<-if(is.null(old$detail))"" else old$detail
   msg(module," SKIP recomputation: verified cache")
  }
  else {
   msg(module," running once")
   result<-if(module %in% names(fun))fun[[module]]() else wf_run_supplement(module)
   if(identical(result,"plot_regenerated"))state<-result
   if(is.list(result)&&identical(result$status,"SKIPPED")) {state<-"SKIPPED";detail<-result$detail;msg(module," SKIPPED: ",detail)}
   if(!all(vapply(file.path(RESULT_ROOT,outputs),wf_validate_artifact,logical(1))))stop("Output validation failed: ",module)
   jsonlite::write_json(list(module=module,key=key,parameters=params,parameter_hash=digest::digest(params,algo="sha256"),
    code_hash=code$hash,code_functions=code$functions,data_versions=versions,outputs=wf_hash_outputs(outputs),validated=TRUE,
    analysis_status=state,detail=detail,
    statistics_key=if(module=="targeted_dependency")wf_targeted_statistics_key() else NULL),
    cache_path,pretty=TRUE,auto_unbox=TRUE)
  }
  runs[[module]]<-data.table(Module=module,Status=state,Execution=execution,Detail=detail,Seconds=as.numeric(difftime(Sys.time(),start,units="secs")))
  gc()
 }
 data.table::fwrite(rbindlist(runs),wf_provenance("Module_Runs.csv"))
 metadata<-list(geneA=GENE_A,geneB=if(GENE_B_PROVIDED)GENE_B else NULL,workflow=WORKFLOW,DATA_MODE=DATA_MODE,
  R_version=as.character(getRversion()),run_date=format(Sys.time(),tz="Asia/Shanghai",usetz=TRUE),
  TCGA_main_layer="gdc_current_DR46",TCGA_current_source="NCI GDC DR46",TCGA_database="TCGA",TCGA_release="46.0",
  DepMap_release="26Q1",DepMap_full_release_completeness="unverified; all supplied Chronos gene columns tested",
  input_versions=inputs,parameters=list(low=.585,deep=.35,min_group_n=MIN_N,bootstrap=BOOT_R,seed=1234),
  code_files=setNames(lapply(list.files("scripts/R",pattern="\\.R$",full.names=TRUE),function(f)digest::digest(file=f,algo="sha256")),list.files("scripts/R",pattern="\\.R$",full.names=TRUE)))
 jsonlite::write_json(metadata,wf_provenance("Run_Metadata.json"),pretty=TRUE,auto_unbox=TRUE)
 writeLines(trimws(capture.output(sessionInfo()),which="right"),wf_provenance("SessionInfo.txt"))
 jsonlite::write_json(inputs,wf_provenance("Input_Versions.json"),pretty=TRUE,auto_unbox=TRUE)
 writeLines(c("Database = TCGA", "Source = NCI GDC", "Release = DR46 / 46.0",
  "TCGA user-facing database: TCGA", "Underlying source: NCI GDC", "Main data layer: gdc_current_DR46",
  "RNA: STAR TPM; analysis: log2(TPM + 1).", "Gene-level CN: current TCGA absolute Gene-Level Copy Number.",
  "TCGA_Relative_CN_Change = CopyNumber / existing BaselineCN - 1; analysis-derived, not GISTIC, log2 GISTIC or DepMap CN_relative.",
  "TCGA 01/02 order: ascending median relative CN change; 02 and 03 x-axis and Pearson/Spearman use this same variable.",
  "Current TCGA CN Five-State Classification: analysis-defined from current GDC DR46 absolute CN and sample-specific baseline/ploidy.",
  "No verified measured ploidy available. Baseline estimate: modal integer autosomal gene-level CN; ties choose smallest mode; never fixed CN=2.",
  "Deep Deletion: CN=0; Shallow Deletion: 0<CN<baseline; Diploid: CN=baseline; Gain: baseline<CN<2*baseline; Amplification: CN>=2*baseline.",
  "These categories are not PanCanAtlas GISTIC -2/-1/0/+1/+2. CNAState numeric codes are internal analysis labels only.",
  "PanCanAtlas/Xena: available locally as historical/reference layer; not used in default 01-03 main results.",
  "Exact SampleID joins, with CaseID/ProjectID consistency checks. RNA representative: matching CN aliquot first, then lexical FileID; audit in Provenance/Data_Audit.",
  "Scatter/pie denominator: current cancer finite matched CN/STAR expression/five-state cohort. Pearson and Spearman BH: separate eligible cancer tests.",
  "Landscape/order: all finite current CN tumor samples; no outlier exclusion.",
  "Source paths, original URLs, source hashes and processor versions: Input_Versions.json.",
  "DepMap: Public 26Q1, user-supplied portal exports; full-release completeness unverified.",
  "Dependency grouping and continuous core: CN_log=log2(relative CN+1), low<0.585.",
  "Targeted scatter statistics: relative CN; historical CN_log correlation retained in table.",
  "Screen delta: median(low)-median(nonlow); BH over the complete supplied screen, including untested NA rows as in supplied code.",
  "Eligible_Rank excludes undefined FDR; original Rank retained for provenance.",
  "All associations are observational; no causal inference."),wf_provenance("Source_Metadata.txt"))
 # A small readable version file accompanies machine-readable source hashes.
 writeLines(c("DATA_MODE=local","TCGA main=NCI GDC DR46 / 46.0","TCGA RNA=STAR log2(TPM + 1)","TCGA reference=historical only; unused by default","DepMap=Public 26Q1",paste("R=",getRversion()),"Input SHA256 and code hashes: Input_Versions.json / Run_Metadata.json / Cache_*.json"),wf_provenance("Input_Versions.txt"))
 wf_summary()
 source(file.path(PROJECT_ROOT,"scripts/R/validate_final_workflow.R"),encoding="UTF-8")
 validate_final_workflow()
 wf_summary() # Include the newly saved validation report in the exact file index.
 msg("Analysis completed: ",RESULT_ROOT)
}

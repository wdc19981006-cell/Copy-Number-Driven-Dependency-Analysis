# Independent recomputation from the exported matched cohort; no workflow reruns.
validate_final_workflow <- function() {
 checks<-list()
 check<-function(name,value) {
  checks[[name]]<<-isTRUE(value)
  if(!isTRUE(value))stop("Final validation failed: ",name)
 }
 equal<-function(a,b) {
  a<-as.numeric(a);b<-as.numeric(b)
  if(length(a)!=length(b)||!identical(is.na(a),is.na(b)))return(FALSE)
  ok<-!is.na(a)
  # Elementwise relative comparison also tests tiny P values: a nonzero 1e-22
  # must never pass merely because an absolute tolerance accepts zero.
  all(abs(a[ok]-b[ok])<=1e-9*pmax(abs(a[ok]),abs(b[ok]),.Machine$double.xmin))
 }
 order<-data.table::fread(wf_table("TCGA_Cancer_Order.csv"));d<-data.table::fread(wf_table("TCGA_Current_Samples.csv"))
 check("R_4_5_0",as.character(getRversion())=="4.5.0")
 check("33_cancer_types",nrow(order)==33L&&setequal(order$CancerType,TCGA_CANCERS))
 check("current_unique_tumors",!anyDuplicated(d$SampleID)&&all(d$TumorNormal=="Tumor")&&all(d$DataLayer=="gdc_current_DR46"))
 check("current_sample_baseline",all(is.finite(d$BaselineCN)&d$BaselineCN>0&d$BaselineCN==floor(d$BaselineCN)))
 expected_state<-with(d,ifelse(!is.finite(CopyNumber),NA_integer_,ifelse(CopyNumber==0,-2L,
   ifelse(CopyNumber<BaselineCN,-1L,ifelse(CopyNumber==BaselineCN,0L,ifelse(CopyNumber<2*BaselineCN,1L,2L))))))
 check("analysis_defined_categories_match_CN_baseline",identical(is.na(d$CNAState),is.na(expected_state))&&all(d$CNAState==expected_state,na.rm=TRUE))
 check("STAR_log2_TPM",equal(d$Expression,log2(d$RNA_TPM+1)))
 versions<-wf_input_paths("tcga")
 check("all_TCGA_inputs_current_processed",all(grepl("data/processed/tcga/gdc_DR46/",versions,fixed=TRUE))&&!any(grepl("PanCanAtlas|Xena|Toil",versions,ignore.case=TRUE)))
 direct<-d[is.finite(CopyNumber),.(N=.N,Median_CN=median(CopyNumber)),by=CancerType]
 setorder(direct,Median_CN,CancerType)
 check("ascending_median_order_and_N",identical(order$CancerType,direct$CancerType)&&identical(order$N,direct$N)&&equal(order$Median_CN,direct$Median_CN))
 plot_order<-jsonlite::fromJSON(wf_provenance("TCGA_Plot_Order.json"))
 check("plot_factors_match_order_table",identical(plot_order$CNA_factor_levels,rev(order$CancerType))&&
       identical(plot_order$Landscape_factor_levels,rev(order$CancerType))&&identical(plot_order$visual_top_to_bottom,order$CancerType))
 check("landscape_all_sample_points",plot_order$landscape_points==sum(is.finite(d$CopyNumber)))
 prev<-data.table::fread(wf_table("TCGA_CNA_Percentage.csv"))
 counts<-d[is.finite(CNAState),.(N=.N),by=.(CancerType,CNAState)]
 observed<-copy(prev);observed[,CNAState:=match(CNA,CNA_STATES)-3L]
 joined<-merge(observed,counts,by=c("CancerType","CNAState"),all.x=TRUE,suffixes=c("_plot","_source"));joined[is.na(N_source),N_source:=0L]
 check("five_state_CNA_counts",nrow(prev)==33L*5L&&all(joined$N_plot==joined$N_source)&&setequal(prev$CNA,CNA_STATES))
 check("CNA_percentage_denominator",all(prev[,sum(N)==unique(denominator)&&abs(sum(Percentage)-100)<1e-8,by=CancerType]$V1))
 stats<-data.table::fread(wf_table(paste0("TCGA_",GENE_A,"_CN_mRNA_AllCancer_Statistics.csv")))
 check("CN_mRNA_summary_order",identical(stats$CancerType,order$CancerType))
 p_pe<-p_sp<-rep(NA_real_,nrow(stats))
 for(i in seq_len(nrow(stats))) {
  z<-d[CancerType==stats$CancerType[i]&is.finite(CopyNumber)&is.finite(Expression)&is.finite(CNAState)]
  check(paste0("RNA_matching_N_",stats$CancerType[i]),stats$N[i]==nrow(z))
  if(nrow(z)>=TCGA_MIN_N&&sd(z$CopyNumber)>0&&sd(z$Expression)>0) {
   pe<-cor.test(z$CopyNumber,z$Expression);sp<-cor.test(z$CopyNumber,z$Expression,method="spearman",exact=FALSE)
   check(paste0("RNA_independent_correlations_",stats$CancerType[i]),
    equal(c(stats$Pearson_r[i],stats$Pearson_P[i],stats$Spearman_rho[i],stats$Spearman_P[i]),c(pe$estimate,pe$p.value,sp$estimate,sp$p.value)))
   p_pe[i]<-pe$p.value;p_sp[i]<-sp$p.value
  }else check(paste0("RNA_insufficient_not_fabricated_",stats$CancerType[i]),all(is.na(unlist(stats[i,.(Pearson_r,Pearson_P,Spearman_rho,Spearman_P)]))))
 }
 check("RNA_independent_BH",equal(stats$Pearson_FDR,p.adjust(p_pe,"BH"))&&equal(stats$Spearman_FDR,p.adjust(p_sp,"BH")))
 pie_counts<-data.table::fread(wf_table("TCGA_CN_mRNA_CNA_Counts.csv"))
 check("pie_uses_matched_cancer_denominator",all(pie_counts[,.(N=sum(N)),by=CancerType][match(stats$CancerType,CancerType),N]==stats$N))
 screen<-data.table::fread(wf_table("GenomeWide_Dependency.csv"))
 check("genomewide_complete_BH",equal(screen$Wilcoxon_FDR,p.adjust(screen$Wilcoxon_P,"BH"))&&equal(screen$Pearson_FDR,p.adjust(screen$Pearson_P,"BH")))
 eligible<-screen[is.finite(Wilcoxon_FDR)&is.finite(Delta_median)][order(Wilcoxon_FDR,Delta_median,Rank)]
 check("eligible_rank",identical(eligible$Eligible_Rank,seq_len(nrow(eligible))))
 # Compare every common numeric field with preserved pre-migration screen if available.
 baseline<-c(file.path(PROJECT_ROOT,"results/VPS4B_VPS4A/02_GenomeWide_Dependency 1/GenomeWide_Dependency.csv"),
             file.path(PROJECT_ROOT,".runtime/prior_case/02_GenomeWide_Dependency/GenomeWide_Dependency.csv"))
 existing<-baseline[file.exists(baseline)]
 if(GENE_A=="VPS4B"&&GENE_B=="VPS4A"&&length(existing)) {
  old<-data.table::fread(existing[1]);s<-screen[match(old$Gene,Gene)]
  numeric_cols<-intersect(names(old)[vapply(old,is.numeric,logical(1))],names(s))
  numeric_cols<-setdiff(numeric_cols,c("Eligible_Rank","Eligible_N"))
  check("original_genomewide_numeric_fields_preserved",nrow(old)==nrow(s)&&all(vapply(numeric_cols,function(n)equal(old[[n]],s[[n]]),logical(1))))
 }
 if(GENE_B_PROVIDED) {
  pair<-data.table::fread(wf_table(paste0(GENE_A,"_",GENE_B,"_CellLines.csv")))
  s<-data.table::fread(wf_table(paste0(GENE_A,"_",GENE_B,"_Targeted_Statistics.csv")))
  check("targeted_unique_finite_cohort",!anyDuplicated(pair$ModelID)&&all(is.finite(pair$CN_relative)&is.finite(pair$Chronos)))
  check("targeted_CN_groups",all((pair$CN_binary=="CN-Low")== (log2(pair$CN_relative+1)<.585)))
  pb<-ggplot_build(targeted_scatter_plot(pair,s))
  check("targeted_three_reference_lines",identical(pb$data[[2]]$yintercept,0)&&
   equal(pb$data[[3]]$xintercept,2^.585-1)&&equal(pb$data[[4]]$xintercept,2^.35-1))
  check("targeted_two_visual_groups",length(unique(pb$data[[1]]$colour))==length(unique(pair$CN_binary)))
  pe<-cor.test(pair$CN_log,pair$Chronos);sp<-cor.test(pair$CN_log,pair$Chronos,method="spearman",exact=FALSE)
  per<-cor.test(pair$CN_relative,pair$Chronos);spr<-cor.test(pair$CN_relative,pair$Chronos,method="spearman",exact=FALSE)
  check("targeted_independent_correlations",equal(c(s$Pearson_r,s$Pearson_P,s$Spearman_rho,s$Spearman_P,
    s$Relative_CN_Pearson_r,s$Relative_CN_Pearson_P,s$Relative_CN_Spearman_rho,s$Relative_CN_Spearman_P),
    c(pe$estimate,pe$p.value,sp$estimate,sp$p.value,per$estimate,per$p.value,spr$estimate,spr$p.value)))
  low<-pair[CN_binary=="CN-Low",Chronos];non<-pair[CN_binary=="CN-NonLow",Chronos]
  if(length(low)>=MIN_N&&length(non)>=MIN_N)check("targeted_independent_Wilcoxon",equal(s$Wilcoxon_P,wilcox.test(low,non,exact=FALSE)$p.value))
  check("targeted_group_medians_and_N",s$N_low==length(low)&&s$N_nonlow==length(non)&&equal(c(s$Median_low,s$Median_nonlow,s$Delta_median),c(median(low),median(non),median(low)-median(non))))
  w<-data.table::fread(wf_table(paste0(GENE_A,"_",GENE_B,"_Waterfall_Order.csv")))
  sorted<-copy(pair);setorder(sorted,-Chronos,ModelID)
  check("waterfall_weak_to_strong",identical(w$ModelID,sorted$ModelID)&&all(diff(w$Chronos)<=0))
  candidate<-screen[Gene==GENE_B]
  if(nrow(candidate))check("screen_targeted_agree",equal(c(candidate$Wilcoxon_P,candidate$Delta_median),c(s$Wilcoxon_P,s$Delta_median)))
 }
 top<-data.table::fread(wf_table("Top_CN_Covariation.csv"))
 check("covariation_direction",all(top[Direction=="Positive CN correlation",Pearson_r]>0)&&all(top[Direction=="Negative CN correlation",Pearson_r]<0)&&!GENE_A %in% top$Gene)
 files<-sort(list.files(RESULT_ROOT,recursive=TRUE))
 expected<-unname(unlist(lapply(workflow_modules(WORKFLOW),wf_expected)))
 check("all_expected_files",all(vapply(file.path(RESULT_ROOT,expected),wf_validate_artifact,logical(1))))
 check("main_PDF_count",length(list.files(file.path(RESULT_ROOT,"Main_Results"),pattern="\\.pdf$",recursive=TRUE))==if(GENE_B_PROVIDED)40L else 37L)
 check("no_default_hidden_modules",!any(grepl("reverse|mutation|expression_dependency|genomewide_adjusted",files,ignore.case=TRUE)))
 check("summary_TXT",file.exists(file.path(RESULT_ROOT,"00_Analysis_Summary.txt"))&&!any(basename(files)=="Summary.md"))
 index<-data.table::fread(wf_provenance("File_Index.csv"))
 check("file_index_matches_actual_files",identical(sort(index$Path),files))
 report<-if(file.exists(wf_provenance("Validation.json")))jsonlite::fromJSON(wf_provenance("Validation.json"),simplifyVector=FALSE) else list()
 report$status<-"PASS";report$checks<-checks;report$number_of_checks<-length(checks)
 jsonlite::write_json(report,wf_provenance("Validation.json"),pretty=TRUE,auto_unbox=TRUE)
 msg("PASS: ",length(checks)," final workflow checks")
 invisible(checks)
}

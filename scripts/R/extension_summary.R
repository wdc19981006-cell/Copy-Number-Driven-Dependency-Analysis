extend_case_summary <- function() {
 prefix<-paste0(GENE_A,"_",GENE_B);folder<-file.path(RESULT_ROOT,"Summary")
 path<-file.path(folder,paste0(prefix,"_Summary.md"));key<-file.path(folder,paste0(prefix,"_Key_Statistics.csv"))
 lines<-if(file.exists(path))readLines(path,warn=FALSE) else character();metrics<-if(file.exists(key))fread(key,colClasses="character") else data.table(module=character(),statistic=character(),value=character())
 sections<-character();read_if<-function(dir,name){p<-file.path(out_dir(dir),name);if(file.exists(p))as.data.frame(readr::read_csv(p,show_col_types=FALSE)) else NULL}
 add<-function(module,row,fields=names(row))for(n in fields)metrics<<-rbind(metrics,data.table(module=module,statistic=n,value=as.character(row[[n]][1])),fill=TRUE)
 fmt<-function(x)format(signif(x,6),scientific=TRUE,trim=TRUE)
 tcga<-read_if("tcga","TCGA_CN_Expression_Adjusted.csv");std<-read_if("tcga","TCGA_CN_Expression_Adjusted_Standardized.csv")
 if(!is.null(tcga)){
  cn<-tcga[tcga$term=="CN",,drop=FALSE];add("TCGA cancer-adjusted CN",cn)
  sections<-c(sections,"","Pan-cancer overall correlation may be influenced by between-cancer differences.",paste("GDC DR46 cancer-adjusted Expression ~ CN + CancerType: CN beta",fmt(cn$estimate),"; P",fmt(cn$p.value),"; N",cn$N,". Association only."))
  if(!is.null(std)){cnz<-std[std$term=="CN",,drop=FALSE];add("TCGA cancer-adjusted standardized CN",cnz);sections<-c(sections,paste("Standardized cancer-adjusted CN beta",fmt(cnz$estimate),"; P",fmt(cnz$p.value)))}
 }
 by<-read_if("tcga","TCGA_CN_Expression_ByCancer.csv")
 if(!is.null(by)){
  eligible<-sum(by$eligible);sig<-sum(by$eligible&!is.na(by$Pearson_FDR)&by$Pearson_FDR<.05)
  sigsp<-sum(by$eligible&!is.na(by$Spearman_FDR)&by$Spearman_FDR<.05)
  add("TCGA per-cancer",data.frame(eligible_cancers=eligible,Pearson_significant_cancers=sig,Spearman_significant_cancers=sigsp))
  sections<-c(sections,paste("TCGA per-cancer N >= 20: eligible",eligible,"; Pearson BH-FDR < 0.05",sig,"; Spearman BH-FDR < 0.05",sigsp,". See ByCancer CSV/forest for all cancers."))
 }
 threshold<-read_if("cn_threshold_sensitivity","CN_Threshold_Sensitivity.csv")
 if(!is.null(threshold)){
  sections<-c(sections,"","CN threshold sensitivity (analysis-defined; not official DepMap GISTIC):")
  for(i in seq_len(nrow(threshold))){a<-threshold[i,,drop=FALSE];add(paste("CN threshold",a$Threshold),a);sections<-c(sections,paste("- Threshold",a$Threshold,": N_low/nonlow",a$N_low,"/",a$N_nonlow,"; delta median",fmt(a$Delta_median),"; P",fmt(a$Wilcoxon_P),"; FDR",fmt(a$Wilcoxon_FDR),"; eligible",a$eligible))}
  continuous<-read_if("cn_threshold_sensitivity","CN_Threshold_Continuous_Statistics.csv");if(!is.null(continuous))add("CN sensitivity continuous (once)",continuous)
 }
 expression<-read_if("expression_dependency","Expression_Dependency.csv");continuous<-read_if("expression_dependency","Expression_Dependency_Continuous.csv")
 if(!is.null(expression)){
  sections<-c(sections,"","Expression-defined dependency (independent of CN-loss and mutation):")
  if(!is.null(continuous)){add("Expression continuous",continuous);sections<-c(sections,paste("Continuous expression -> Chronos: N",continuous$N,"; Pearson",fmt(continuous$Pearson_r),"P",fmt(continuous$Pearson_P),"; Spearman",fmt(continuous$Spearman_rho),"P",fmt(continuous$Spearman_P)))}
  for(i in seq_len(nrow(expression))){a<-expression[i,,drop=FALSE];add(paste("Expression quantile",a$threshold_quantile),a);sections<-c(sections,paste("- Bottom",100*a$threshold_quantile,"%: cutoff",fmt(a$expression_cutoff),"; N_low/nonlow",a$N_low,"/",a$N_nonlow,"; delta median",fmt(a$Delta_median),"; P",fmt(a$Wilcoxon_P),"; FDR",fmt(a$Wilcoxon_FDR),"; eligible",a$eligible))}
 }
 adjusted<-read_if("genomewide_adjusted_dependency",paste0(GENE_B,"_Adjusted_Candidate.csv"))
 if(!is.null(adjusted)&&nrow(adjusted)){
  add("Genomewide lineage-adjusted candidate",adjusted)
  sections<-c(sections,"",paste(GENE_B,"lineage-adjusted genome-wide: Eligible_Rank_adjusted",adjusted$Eligible_Rank_adjusted,"of",adjusted$Eligible_N_adjusted,"; Beta_CN",fmt(adjusted$Beta_CN),"; FDR_CN",fmt(adjusted$FDR_CN),"; Beta_CNLow",fmt(adjusted$Beta_CNLow),"; FDR_CNLow",fmt(adjusted$FDR_CNLow)))
 }
 writeLines(c(lines,sections),path,useBytes=TRUE);fwrite(metrics,key)
}

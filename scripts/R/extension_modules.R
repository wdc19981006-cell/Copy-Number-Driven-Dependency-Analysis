# Modules use existing Parquet, leaving the supplied raw-CSV statistical core intact.
processed_depmap <- function(dataset,gene=NULL,value_name=NULL) {
 path<-file.path(PROJECT_ROOT,"data/processed/depmap/26Q1",paste0(dataset,".parquet"))
 if(!file.exists(path))stop("Existing processed DepMap file missing: ",path)
 reader<-arrow::ParquetFileReader$create(path);nms<-reader$GetSchema()$names
 if(!is.null(gene)){
  column<-find_gene_col(nms,gene)
  dat<-as.data.frame(arrow::read_parquet(path,col_select=all_of(c("ModelID",column))))
  names(dat)<-c("ModelID",value_name)
 }else dat<-as.data.frame(arrow::read_parquet(path))
 if(anyDuplicated(dat$ModelID))stop("Duplicate processed ModelID in ",dataset)
 dat
}
processed_pair <- function(kind="cn") {
 x<-processed_depmap(if(kind=="cn")"OmicsCNGeneWGS" else "ExpressionProteinCoding",GENE_A,if(kind=="cn")"CN_relative" else "Expression")
 y<-processed_depmap("CRISPRGeneEffect",GENE_B,"Chronos")
 dat<-merge(x,y,by="ModelID",sort=FALSE)
 if(kind=="cn")dat$CN_log<-log2(dat$CN_relative+1)
 dat
}
run_threshold_sensitivity <- function() {
 d<-processed_pair();f<-out_dir("cn_threshold_sensitivity")
 res<-cn_threshold_statistics(d$CN_log,d$Chronos)
 res$Definition<-"analysis-defined threshold; not DepMap official GISTIC classification"
 fwrite(res,file.path(f,"CN_Threshold_Sensitivity.csv"))
 fwrite(pair_correlations(d$CN_log,d$Chronos),file.path(f,"CN_Threshold_Continuous_Statistics.csv"))
 plotted<-res[is.finite(res$Delta_median),,drop=FALSE]
 skipped<-res[!res$eligible,,drop=FALSE]
 caption<-if(nrow(skipped))paste("Not tested:",paste(paste0("threshold ",skipped$Threshold," (N_low=",skipped$N_low,", N_nonlow=",skipped$N_nonlow,")"),collapse="; ")) else NULL
 p<-ggplot(plotted,aes(Threshold,Delta_median))+geom_hline(yintercept=0,linetype="dashed")+geom_line()+geom_point(aes(shape=eligible),size=3)+
  geom_text(aes(label=paste0("N_low=",N_low,ifelse(!is.na(Wilcoxon_FDR)&Wilcoxon_FDR<.05," *","")),hjust=ifelse(Threshold==max(res$Threshold),1,.5)),vjust=-1,size=3)+
  scale_x_continuous(breaks=sort(res$Threshold),limits=range(res$Threshold))+theme_classic()+labs(title=paste(GENE_A,"->",GENE_B,"CN threshold sensitivity"),subtitle="Analysis-defined CN_log thresholds; * BH-FDR < 0.05",caption=paste(strwrap(caption,width=88),collapse="\n"),y="Delta median Chronos (low - non-low)",x="CN_log threshold")
 save_pdf(p,f,"CN_Threshold_Sensitivity",width=8)
}
run_expression_dependency <- function() {
 d<-processed_pair("expression");good<-is.finite(d$Expression)&is.finite(d$Chronos);d<-d[good,,drop=FALSE]
 f<-out_dir("expression_dependency");res<-expression_dependency_statistics(d$Expression,d$Chronos)
 fwrite(res$groups,file.path(f,"Expression_Dependency.csv"));fwrite(res$continuous,file.path(f,"Expression_Dependency_Continuous.csv"))
 fwrite(d,file.path(f,"Expression_Dependency_Samples.csv"))
 plots<-list(ggplot(d,aes(Expression,Chronos))+geom_point(alpha=.3,size=.8)+geom_smooth(method="lm",se=FALSE)+theme_classic()+
  labs(title=paste("A:",GENE_A,"expression ->",GENE_B),subtitle=paste0("N=",res$continuous$N,"; Pearson r=",signif(res$continuous$Pearson_r,3),"; Spearman rho=",signif(res$continuous$Spearman_rho,3)),x="Supplied expression scale",y="Chronos"))
 for(i in seq_len(nrow(res$groups))){
  row<-res$groups[i,];if(!row$eligible)next
  d$Group<-factor(ifelse(d$Expression<=row$expression_cutoff,"Expression-Low","Expression-NonLow"),levels=c("Expression-Low","Expression-NonLow"))
  plots[[length(plots)+1L]]<-ggplot(d,aes(Group,Chronos,fill=Group))+geom_boxplot(outlier.shape=NA)+geom_jitter(width=.12,alpha=.2,size=.7)+theme_classic()+theme(legend.position="none")+
   labs(title=paste0(LETTERS[i+1],": bottom ",100*row$threshold_quantile,"% expression"),subtitle=paste0("N=",row$N_low,"/",row$N_nonlow,"; P=",format(row$Wilcoxon_P,digits=3,scientific=TRUE),"; FDR=",format(row$Wilcoxon_FDR,digits=3,scientific=TRUE)),x=NULL,y="Chronos")
 }
 save_pdf(patchwork::wrap_plots(plots,nrow=1),f,"Expression_Dependency",width=5*length(plots),height=5.5)
}
run_genomewide_adjusted <- function() {
 cn<-processed_depmap("OmicsCNGeneWGS",GENE_A,"CN_relative");cn$CN_log<-log2(cn$CN_relative+1)
 model<-processed_depmap("Model");meta<-merge(cn,model[,c("ModelID","OncotreeLineage")],by="ModelID",sort=FALSE)
 matrix<-processed_depmap("CRISPRGeneEffect");names(matrix)[-1]<-clean_gene(names(matrix)[-1])
 if(anyDuplicated(names(matrix)))stop("Ambiguous Chronos gene symbols")
 res<-genomewide_adjusted_statistics(matrix,meta);f<-out_dir("genomewide_adjusted_dependency")
 fwrite(res,file.path(f,"GenomeWide_Adjusted_Dependency.csv"))
 candidate<-res[res$Gene==GENE_B,,drop=FALSE]
 if(GENE_B_PROVIDED)fwrite(candidate,file.path(f,paste0(GENE_B,"_Adjusted_Candidate.csv")))
 eligible<-res[is.finite(res$FDR_CN)&is.finite(res$Beta_CN),,drop=FALSE];eligible$logFDR<--log10(pmax(eligible$FDR_CN,.Machine$double.xmin))
 p<-ggplot(eligible,aes(Beta_CN,logFDR))+geom_point(alpha=.3,size=.8)+geom_hline(yintercept=-log10(.05),linetype="dashed")+theme_classic()+labs(title=paste(GENE_A,"lineage-adjusted genome-wide dependency"),x="CN_log beta (lineage adjusted)",y="-log10(BH FDR_CN)")
 if(GENE_B_PROVIDED&&nrow(candidate)&&is.finite(candidate$FDR_CN))p<-p+geom_point(data=subset(eligible,Gene==GENE_B),color="#C7473B",size=3)+ggrepel::geom_text_repel(data=subset(eligible,Gene==GENE_B),aes(label=Gene),color="#C7473B")
 save_pdf(p,f,"GenomeWide_Adjusted_Dependency_Volcano",width=8)
}
run_genomewide_expression <- function() {
 ex<-processed_depmap("ExpressionProteinCoding",GENE_A,"Expression");matrix<-processed_depmap("CRISPRGeneEffect")
 names(matrix)[-1]<-clean_gene(names(matrix)[-1]);if(anyDuplicated(names(matrix)))stop("Ambiguous Chronos gene symbols")
 res<-genomewide_expression_statistics(matrix,ex)
 fwrite(res,file.path(out_dir("genomewide_expression_dependency"),"GenomeWide_Expression_Dependency.csv"))
}
run_tcga_expression_extensions <- function() {
 path<-file.path(out_dir("tcga"),"TCGA_CN_Expression_Samples.csv");d<-fread(path)
 if(anyDuplicated(d$SampleID))stop("Duplicate current sample UUID")
 if(!all(d$GDCRelease_CN=="DR46"&d$GDCRelease_RNA=="DR46"))stop("Require GDC current DR46")
 if(!all(d$CancerType_CN==d$CancerType_RNA))stop("CN/RNA cancer mapping disagreement")
 dat<-data.frame(CN=d$Value_CN,Expression=d$Value_RNA,CancerType=d$CancerType_CN)
 f<-out_dir("tcga");res<-tcga_cancer_statistics(dat);res$GDCRelease<-"DR46";res$DataLayer<-"gdc_DR46"
 fwrite(res,file.path(f,"TCGA_CN_Expression_ByCancer.csv"))
 plotted<-res[res$eligible&is.finite(res$Pearson_r),,drop=FALSE]
 plotted$CancerType<-factor(plotted$CancerType,levels=plotted$CancerType[order(plotted$Pearson_r)])
 p<-ggplot(plotted,aes(Pearson_r,CancerType))+geom_vline(xintercept=0,linetype="dashed")+
  geom_errorbar(aes(xmin=Pearson_CI_low,xmax=Pearson_CI_high),orientation="y",width=.2)+geom_point(color="#3A6EA5")+theme_classic()+
  labs(title=paste(GENE_A,"GDC DR46 CN-expression by cancer"),subtitle="N >= 20; Fisher-z 95% CI",x="Pearson r",y=NULL)
 save_pdf(p,f,"TCGA_CN_Expression_ByCancer_Forest",height=9,width=8)
 for(standardized in c(FALSE,TRUE)){
  result<-tcga_adjusted_statistics(dat,standardized)
  if(!result$eligible)stop(result$reason)
  tab<-result$coefficients;tab$N<-result$N;tab$GDCRelease<-"DR46";tab$DataLayer<-"gdc_DR46"
  fwrite(tab,file.path(f,if(standardized)"TCGA_CN_Expression_Adjusted_Standardized.csv" else "TCGA_CN_Expression_Adjusted.csv"))
 }
 writeLines(c("GDC current DR46. Expression = log2(TPM + 1); CN = source total gene-level CN.",
  "Pan-cancer overall correlation may be influenced by between-cancer differences.",
  "Report overall correlations alongside cancer-adjusted CN beta; association does not establish causality.",
  "Per-cancer correlations require N >= 20 and nonconstant CN/expression. BH families are separate for Pearson and Spearman.",
  "Fisher-z CI = tanh(atanh(r) +/- qnorm(0.975)/sqrt(N-3)).",
  "Adjusted regressions require known cancer type; standardization uses global sample SD on that same complete regression cohort."),file.path(f,"TCGA_CN_Expression_Methods.md"))
}
# Extend the existing module without replacing its overall statistics or mapping.
run_tcga_expression_original <- run_tcga_expression
run_tcga_expression <- function() {
 run_tcga_expression_original()
 if(is.null(MODULE_SKIP))run_tcga_expression_extensions()
}

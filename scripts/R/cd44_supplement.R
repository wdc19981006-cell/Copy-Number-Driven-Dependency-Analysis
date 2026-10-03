#!/usr/bin/env Rscript
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
if (.Platform$OS.type == "windows") Sys.setenv(PROCESSOR_ARCHITECTURE="AMD64")
stopifnot(as.character(getRversion()) == "4.5.0")
source(".Rprofile")
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(arrow); library(dplyr)})
data.table::setDTthreads(4)
source("scripts/R/extensions_statistics.R")
out <- "results/CD44_PanCancer/13_CD44_Supplement"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
write <- function(x, name) fwrite(x, file.path(out, paste0(name, ".csv")))
plot_save <- function(p, name, w=10, h=7) {
  ggsave(file.path(out, paste0(name,".png")), p, width=w, height=h, dpi=180, bg="white")
}
states <- c("Deep deletion", "Shallow deletion", "Diploid", "Gain", "Amplification")
colors <- setNames(c("#2166ac", "#67a9cf", "#cccccc", "#ef8a62", "#b2182b"), states)
ref <- fread(file.path(out, "TCGA_Reference_CNA_RNA_Samples.csv"))
ref[,CNA := factor(GISTIC, levels=-2:2, labels=states)]
stopifnot(!anyDuplicated(ref$SampleID), all(ref$TumorNormal=="Tumor"))
group_summary <- function(d) d[,.(N=.N, RNA_median=median(Expression), RNA_Q1=quantile(Expression,.25),
                                RNA_Q3=quantile(Expression,.75)), by=CNA]
overall <- group_summary(ref)
write(overall, "Reference_CNA_RNA_Overall_Groups")
by_groups <- ref[,group_summary(.SD), by=CancerType]
write(by_groups, "Reference_CNA_RNA_ByCancer_Groups")
fit <- lm(Expression ~ relevel(CNA,ref="Diploid") + CancerType, data=ref)
tab <- as.data.table(coefficient_table(fit)); tab[,N:=nrow(ref)]
tab[,P_CNA_BH := NA_real_]
tab[grepl("^relevel",term),P_CNA_BH := p.adjust(p.value,"BH")]
write(tab, "Reference_CNA_RNA_CancerAdjusted")
contrasts <- list()
for (cancer in sort(unique(ref$CancerType))) {
  d <- ref[CancerType==cancer]
  control <- d[GISTIC==0, Expression]
  for (code in c(-2,-1,1,2)) {
    a <- d[GISTIC==code, Expression]
    eligible <- length(a)>=3 && length(control)>=3
    contrasts[[length(contrasts)+1L]] <- data.table(CancerType=cancer, CNA=states[code+3],
      N_CNA=length(a), N_diploid=length(control), eligible=eligible,
      Delta_median=if(length(a)&&length(control))median(a)-median(control) else NA_real_,
      Wilcoxon_P=if(eligible)wilcox.test(a,control,exact=FALSE)$p.value else NA_real_)
  }
}
contrasts <- rbindlist(contrasts); contrasts[,Wilcoxon_FDR:=p.adjust(Wilcoxon_P,"BH")]
write(contrasts, "Reference_CNA_RNA_ByCancer_vsDiploid")
cont <- data.frame(CN=ref$CopyNumber, Expression=ref$Expression, CancerType=ref$CancerType)
write(pair_correlations(cont$CN,cont$Expression), "Reference_Continuous_CN_RNA")
write(tcga_cancer_statistics(cont), "Reference_Continuous_CN_RNA_ByCancer")
labels <- setNames(paste0(states,"\nN=",overall$N[match(states,as.character(overall$CNA))]),states)
p <- ggplot(ref,aes(CNA,Expression,fill=CNA))+geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=colors)+scale_x_discrete(labels=labels)+theme_classic()+theme(legend.position="none")+
  labs(title="CD44 copy-number state and RNA expression",subtitle="TCGA reference tumors; pooled descriptive groups",x=NULL,
       y="CD44 log2(norm_count + 1)",caption="Five GISTIC states; not absolute copy counts. See cancer-adjusted and within-cancer tests.")
plot_save(p,"CD44_CNA_RNA_Overall",10,6)
p <- ggplot(ref,aes(CNA,Expression,fill=CNA))+geom_boxplot(outlier.shape=NA)+facet_wrap(~CancerType,ncol=6)+
  scale_fill_manual(values=colors,drop=FALSE)+theme_bw(base_size=9)+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom")+
  labs(title="CD44 RNA by GISTIC state in each TCGA cancer",x=NULL,y="log2(norm_count + 1)",fill=NULL)
plot_save(p,"CD44_CNA_RNA_ByCancer",16,13)

# Frequencies use all profiled tumors, including those without matched RNA.
prev <- fread("results/CD44_PanCancer/08_TCGA/TCGA_GISTIC_CNA_Prevalence.csv")
grid <- CJ(CancerType=unique(prev$CancerType),CNA=states,unique=TRUE)
prev <- merge(grid,prev,by=c("CancerType","CNA"),all.x=TRUE)
prev[is.na(N),N:=0L]; prev[,denominator:=sum(N),by=CancerType];prev[,prevalence:=N/denominator]
prev[,CNA:=factor(CNA,levels=states)]
write(prev,"Reference_CNA_Prevalence_Complete")
freq <- dcast(prev,CancerType+denominator~CNA,value.var="N")
freq[,Loss_pct:=100*(`Deep deletion`+`Shallow deletion`)/denominator]
freq[,GainAmplification_pct:=100*(Gain+Amplification)/denominator]
setorder(freq,-Loss_pct);write(freq,"Reference_CNA_Prevalence_ByCancer")
global <- prev[,.(N=sum(N)),by=CNA];global[,prevalence:=N/sum(N)]
write(global,"Reference_CNA_Prevalence_Overall")
prev[,CancerType:=factor(CancerType,levels=rev(freq$CancerType))]
p<-ggplot(prev,aes(CancerType,prevalence,fill=CNA))+geom_col()+coord_flip()+
 scale_fill_manual(values=colors,drop=FALSE)+scale_y_continuous(labels=scales::percent)+theme_classic()+
 labs(title="CD44 copy-number landscape",subtitle=paste(sum(global$N),"TCGA reference tumor samples; ordered by deletion frequency"),x=NULL,y="Fraction of tumors",fill=NULL)
plot_save(p,"CD44_CNA_Prevalence",11,10)

# Current continuous-CN analysis is restricted explicitly to tumor specimens.
current <- fread("results/CD44_PanCancer/08_TCGA/TCGA_CN_Expression_Samples.csv")
types <- current[,.(N=.N),by=SampleType_CN];write(types,"Current_SampleTypes_Audit")
tumor <- current[grepl("Primary|Tumor|Cancer|Metastatic",SampleType_CN,ignore.case=TRUE)]
stopifnot(!anyDuplicated(tumor$SampleID),all(tumor$CancerType_CN==tumor$CancerType_RNA))
write(tumor,"Current_Tumor_CN_RNA_Samples")
d <- data.frame(CN=tumor$Value_CN,Expression=tumor$Value_RNA,CancerType=tumor$CancerType_CN)
write(pair_correlations(d$CN,d$Expression),"Current_Tumor_CN_RNA_Overall")
by <- tcga_cancer_statistics(d);write(by,"Current_Tumor_CN_RNA_ByCancer")
for (z in c(FALSE,TRUE)) {
  result<-tcga_adjusted_statistics(d,z);stopifnot(result$eligible)
  t<-result$coefficients;t$N<-result$N
  write(t,if(z)"Current_Tumor_CN_RNA_Adjusted_Standardized" else "Current_Tumor_CN_RNA_Adjusted")
}
# One sample per case is an explicit sensitivity analysis; tumor priority is audited.
tumor[,priority:=ifelse(grepl("Primary",SampleType_CN,ignore.case=TRUE),0L,1L)]
setorder(tumor,CaseID_CN,priority,SampleID);tumor[,SelectedCaseRepresentative:=!duplicated(CaseID_CN)]
write(tumor[,.(SampleID,CaseID_CN,SampleType_CN,SelectedCaseRepresentative)],"Current_Case_Representative_Audit")
case <- tumor[SelectedCaseRepresentative==TRUE]
cd<-data.frame(CN=case$Value_CN,Expression=case$Value_RNA,CancerType=case$CancerType_CN)
write(pair_correlations(cd$CN,cd$Expression),"Current_OneTumorPerCase_CN_RNA")
ct<-tcga_adjusted_statistics(cd)$coefficients;ct$N<-nrow(cd);write(ct,"Current_OneTumorPerCase_Adjusted")
p<-ggplot(by[by$eligible,],aes(Pearson_r,reorder(CancerType,Pearson_r)))+
 geom_vline(xintercept=0,linetype="dashed")+geom_errorbar(aes(xmin=Pearson_CI_low,xmax=Pearson_CI_high),orientation="y",width=.2)+
 geom_point(aes(color=Pearson_FDR<.05),size=2)+theme_classic()+labs(title="CD44 copy number and RNA: current TCGA tumors",x="Pearson r with 95% CI",y=NULL,color="BH FDR < 0.05")
plot_save(p,"CD44_Current_Tumor_CN_RNA_ByCancer",9,10)

# Select candidates transparently from the complete screening families.
raw <- fread("results/CD44_PanCancer/02_GenomeWide_Dependency/GenomeWide_Dependency.csv")
adj <- fread("results/CD44_PanCancer/12_GenomeWide_Adjusted_Dependency/GenomeWide_Adjusted_Dependency.csv")
merged <- merge(raw,adj,by="Gene",suffixes=c("_raw","_adjusted"))
merged[,SupportsLowCN := Delta_median<0 & Beta_CNLow<0]
merged[,JointBinaryFDR05 := SupportsLowCN & Wilcoxon_FDR<.05 & FDR_CNLow<.05]
merged[,ContinuousLowCNFDR05 := Beta_CN>0 & FDR_CN<.05]
setorder(merged,FDR_CNLow,Delta_median,na.last=TRUE)
write(merged,"CD44_AllGenes_CombinedScreen")
write(merged[JointBinaryFDR05==TRUE],"CD44_LowCN_JointSignificant_Candidates")
binary <- merged[SupportsLowCN==TRUE];setorder(binary,FDR_CNLow,Wilcoxon_FDR,Delta_median)
write(head(binary,30),"CD44_LowCN_Exploratory_Top30")
continuous<-merged[Beta_CN>0];setorder(continuous,FDR_CN)
write(head(continuous,30),"CD44_Continuous_LowCN_Top30")
genes <- unique(c(head(binary$Gene,10),head(continuous$Gene,10)))
read_gene <- function(dataset,gene,value) {
 path<-file.path("data/processed/depmap/26Q1",paste0(dataset,".parquet"))
 nms<-ParquetFileReader$create(path)$GetSchema()$names
 hit<-nms[sub("\\s*\\([^)]*\\)\\s*$","",nms)==gene];stopifnot(length(hit)==1L)
 x<-as.data.table(read_parquet(path,col_select=all_of(c("ModelID",hit))))
 setnames(x,c("ModelID",value));x
}
cn<-read_gene("OmicsCNGeneWGS","CD44","CN_relative");cn[,CN_log:=log2(CN_relative+1)]
ex<-read_gene("ExpressionProteinCoding","CD44","CD44_RNA")
model<-as.data.table(read_parquet("data/processed/depmap/26Q1/Model.parquet"))
meta<-merge(merge(cn,ex,by="ModelID",all.x=TRUE),model[,.(ModelID,CellLineName,OncotreeLineage)],by="ModelID")
details<-list();sensitivity<-list();lineage<-list();loo<-list()
for(g in genes) {
 x<-merge(meta,read_gene("CRISPRGeneEffect",g,"Chronos"),by="ModelID")
 x<-x[is.finite(CN_log)&is.finite(Chronos)&!is.na(OncotreeLineage)]
 x[,CNLow:=CN_log<.585];x[,Gene:=g];details[[g]]<-copy(x)
 r<-cn_threshold_statistics(x$CN_log,x$Chronos);r$Gene<-g;write(r,paste0(g,"_ThresholdSensitivity"))
 for(l in sort(unique(x$OncotreeLineage))) {
   a<-x[OncotreeLineage==l];row<-group_effect(a$Chronos,a$CNLow)
   row$Gene<-g;row$Lineage<-l;lineage[[length(lineage)+1L]]<-row
 }
 for(id in x[CNLow==TRUE,ModelID]) {
   a<-x[ModelID!=id];row<-group_effect(a$Chronos,a$CNLow)
   row$Gene<-g;row$OmittedLowModel<-id
   fit<-tryCatch(lm(Chronos~CNLow+OncotreeLineage,data=a),error=function(e)NULL)
   row$AdjustedBeta<-if(is.null(fit))NA_real_ else coef(fit)["CNLowTRUE"]
   loo[[length(loo)+1L]]<-row
 }
 # Own target CN and CD44 RNA are candidate sensitivity covariates, not screening families.
 tcn<-tryCatch(read_gene("OmicsCNGeneWGS",g,"TargetCN"),error=function(e)NULL)
 if(!is.null(tcn))x<-merge(x,tcn,by="ModelID") else x[,TargetCN:=NA_real_]
 x[,TargetCNLog:=log2(TargetCN+1)]
 forms<-list(base=Chronos~CN_log+OncotreeLineage,
   exclude_top1pct=Chronos~CN_log+OncotreeLineage,
   plus_targetCN=Chronos~CN_log+TargetCNLog+OncotreeLineage,
   plus_CD44_RNA=Chronos~CN_log+CD44_RNA+OncotreeLineage)
 for(name in names(forms)) {
   a<-if(name=="exclude_top1pct")x[CN_log<=quantile(CN_log,.99)] else x
   variables<-all.vars(forms[[name]])
   a<-a[complete.cases(a[,..variables])]
   eligible<-nrow(a)>=20L && uniqueN(a$OncotreeLineage)>=2L && sd(a$CN_log)>0 && sd(a$Chronos)>0
   fit<-if(eligible)tryCatch(lm(forms[[name]],data=a),error=function(e)NULL) else NULL
   row<-data.frame(term="CN_log",estimate=NA_real_,std.error=NA_real_,statistic=NA_real_,p.value=NA_real_)
   if(!is.null(fit)) {
     t<-coefficient_table(fit);hit<-t[t$term=="CN_log",]
     if(nrow(hit))row<-hit
   }
   row$Gene<-g;row$Model<-name;row$N<-nrow(a);row$eligible<-!is.null(fit)&&is.finite(row$estimate)
   row$reason<-if(row$eligible)"" else "Insufficient complete data, constant outcome, or aliased CN term"
   sensitivity[[length(sensitivity)+1L]]<-row
 }
}
details<-rbindlist(details);write(details,"Candidate_Dependency_CellLines")
write(rbindlist(sensitivity),"Candidate_Continuous_Sensitivity")
write(rbindlist(lineage,fill=TRUE),"Candidate_WithinLineage_Effects")
write(rbindlist(loo,fill=TRUE),"Candidate_LeaveOneLowModelOut")
selected<-head(binary$Gene,6)
p<-ggplot(details[Gene %in% selected],aes(factor(CNLow,levels=c(FALSE,TRUE)),Chronos,fill=CNLow))+
 geom_boxplot(outlier.shape=NA)+geom_jitter(width=.10,alpha=.18,size=.6)+facet_wrap(~Gene,ncol=3,scales="free_y")+
 scale_x_discrete(labels=c("Non-low","Low"))+theme_classic()+theme(legend.position="none")+
 labs(title="CD44 low-copy-number exploratory dependency candidates",subtitle="Candidates ranked by lineage-adjusted binary FDR; inspect full FDR table",x="CD44 copy-number group",y="Chronos gene effect (more negative = stronger dependency)")
plot_save(p,"CD44_Dependency_Candidates",12,8)
checks<-list(reference_unique=TRUE,current_unique=TRUE,reference_tumors=nrow(ref),current_tumors=nrow(d),
  current_cases=nrow(cd),screened_genes=nrow(merged),joint_significant_binary=sum(merged$JointBinaryFDR05,na.rm=TRUE),
  continuous_low_CN_significant=sum(merged$ContinuousLowCNFDR05,na.rm=TRUE),
  raw_BH_verified=isTRUE(all.equal(raw$Wilcoxon_FDR,p.adjust(raw$Wilcoxon_P,"BH"),tolerance=1e-10)),
  adjusted_BH_verified=isTRUE(all.equal(adj$FDR_CNLow,p.adjust(adj$P_CNLow,"BH"),tolerance=1e-10)))
stopifnot(checks$raw_BH_verified,checks$adjusted_BH_verified,nrow(raw)==18531,nrow(adj)==18531)
jsonlite::write_json(checks,file.path(out,"Supplement_Validation.json"),pretty=TRUE,auto_unbox=TRUE)
print(checks)

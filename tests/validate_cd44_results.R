invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
if(.Platform$OS.type=="windows")Sys.setenv(PROCESSOR_ARCHITECTURE="AMD64")
source(".Rprofile")
suppressPackageStartupMessages({library(data.table);library(arrow);library(dplyr)})
stopifnot(as.character(getRversion())=="4.5.0")
base<-"results/CD44_PanCancer";out<-file.path(base,"13_CD44_Supplement")
read<-function(name)fread(file.path(out,paste0(name,".csv")))
near<-function(a,b)stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=1e-10)))
near_p<-function(a,b){a<-as.numeric(a);b<-as.numeric(b);stopifnot(identical(is.na(a),is.na(b)));ok<-is.finite(a)&is.finite(b)&a>0&b>0;stopifnot(all(abs(log(a[ok])-log(b[ok]))<1e-8))}
raw<-fread(file.path(base,"02_GenomeWide_Dependency/GenomeWide_Dependency.csv"))
adj<-fread(file.path(base,"12_GenomeWide_Adjusted_Dependency/GenomeWide_Adjusted_Dependency.csv"))
stopifnot(nrow(raw)==18531L,nrow(adj)==18531L,!anyDuplicated(raw$Gene),!anyDuplicated(adj$Gene))
near(raw$Wilcoxon_FDR,p.adjust(raw$Wilcoxon_P,"BH"));near(raw$Pearson_FDR,p.adjust(raw$Pearson_P,"BH"))
near(adj$FDR_CN,p.adjust(adj$P_CN,"BH"));near(adj$FDR_CNLow,p.adjust(adj$P_CNLow,"BH"))
eligible<-raw[!is.na(Wilcoxon_FDR)];setorder(eligible,Wilcoxon_FDR,Delta_median,Rank)
stopifnot(all(eligible$Eligible_Rank==seq_len(nrow(eligible))),all(eligible$Eligible_N==nrow(eligible)))
gene<-function(dataset,g,value){
 path<-file.path("data/processed/depmap/26Q1",paste0(dataset,".parquet"))
 nms<-ParquetFileReader$create(path)$GetSchema()$names
 col<-nms[sub("\\s*\\([^)]*\\)\\s*$","",nms)==g];stopifnot(length(col)==1L)
 d<-as.data.table(read_parquet(path,col_select=all_of(c("ModelID",col))));setnames(d,c("ModelID",value));d
}
cn<-gene("OmicsCNGeneWGS","CD44","CN_relative");cn[,CN_log:=log2(CN_relative+1)]
model<-as.data.table(read_parquet("data/processed/depmap/26Q1/Model.parquet"))
meta<-merge(cn,model[,.(ModelID,OncotreeLineage)],by="ModelID")
for(g in c("AOC1","VWA5B2","PSMD10","SNAP23","RPL6","SMIM11")){
 d<-merge(meta,gene("CRISPRGeneEffect",g,"Chronos"),by="ModelID")
 d<-d[is.finite(CN_log)&is.finite(Chronos)&!is.na(OncotreeLineage)]
 low<-d$CN_log<.585;r<-raw[Gene==g];a<-adj[Gene==g]
 stopifnot(nrow(d)==r$N,nrow(d)==a$N,sum(low)==r$N_low,sum(!low)==r$N_nonlow)
 if(sum(low)>=3L&&sum(!low)>=3L){
  near(r$Delta_median,median(d$Chronos[low])-median(d$Chronos[!low]))
  near_p(r$Wilcoxon_P,wilcox.test(d$Chronos[low],d$Chronos[!low],exact=FALSE)$p.value)
 }
 f<-summary(lm(Chronos~CN_log+OncotreeLineage,data=d))$coefficients["CN_log",]
 near(c(a$Beta_CN,a$SE_CN,a$T_CN),f[1:3]);near_p(a$P_CN,f[4])
 f<-summary(lm(Chronos~I(CN_log<.585)+OncotreeLineage,data=d))$coefficients["I(CN_log < 0.585)TRUE",]
 near(c(a$Beta_CNLow,a$SE_CNLow,a$T_CNLow),f[1:3]);near_p(a$P_CNLow,f[4])
}
ref<-read("TCGA_Reference_CNA_RNA_Samples");stopifnot(!anyDuplicated(ref$SampleID),all(ref$TumorNormal=="Tumor"),nrow(ref)==9492L)
prev<-read("Reference_CNA_Prevalence_Complete");stopifnot(all(prev[,sum(N)==unique(denominator),by=CancerType]$V1),sum(prev$N)==10845L)
contrast<-read("Reference_CNA_RNA_ByCancer_vsDiploid")
near(contrast$Wilcoxon_FDR,p.adjust(contrast$Wilcoxon_P,"BH"))
for(i in which(contrast$eligible)){
 row<-contrast[i];code<-match(row$CNA,c("Deep deletion","Shallow deletion","Diploid","Gain","Amplification"))-3L
 d<-ref[CancerType==row$CancerType];a<-d[GISTIC==code,Expression];b<-d[GISTIC==0,Expression]
 stopifnot(length(a)==row$N_CNA,length(b)==row$N_diploid)
 near(row$Delta_median,median(a)-median(b));near_p(row$Wilcoxon_P,wilcox.test(a,b,exact=FALSE)$p.value)
}
current<-read("Current_Tumor_CN_RNA_Samples");stopifnot(nrow(current)==10549L,!anyDuplicated(current$SampleID),!any(grepl("Normal",current$SampleType_CN)))
by<-read("Current_Tumor_CN_RNA_ByCancer")
for(i in which(by$eligible)){
 d<-current[CancerType_CN==by$CancerType[i]];p<-cor.test(d$Value_CN,d$Value_RNA)
 stopifnot(nrow(d)==by$N[i]);near(by$Pearson_r[i],p$estimate);near_p(by$Pearson_P[i],p$p.value)
}
near(by$Pearson_FDR,p.adjust(by$Pearson_P,"BH"))
adjusted<-read("Current_Tumor_CN_RNA_Adjusted")
f<-summary(lm(Value_RNA~Value_CN+CancerType_CN,data=current))$coefficients["Value_CN",]
near(adjusted[term=="CN",.(estimate,std.error,statistic)],f[1:3]);near_p(adjusted[term=="CN",p.value],f[4])
checks<-c("18531 genes and all four complete-screen BH families","eligible ranking","six candidates recomputed from source Parquet",
 "reference tumors, five-state denominators and every eligible within-cancer RNA comparison","current tumor-only sample identity, per-cancer correlations and adjusted coefficient")
jsonlite::write_json(list(all_passed=TRUE,checks=checks,R_version=as.character(getRversion())),file.path(out,"Independent_Validation.json"),pretty=TRUE,auto_unbox=TRUE)
cat("PASS: independent CD44 source alignment, candidate statistics, FDR, TCGA identities and RNA contrasts.\n")

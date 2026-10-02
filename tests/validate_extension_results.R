invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"));source(".Rprofile");library(data.table)
stopifnot(as.character(getRversion())=="4.5.0")
root<-"results/VPS4B_VPS4A"
read<-function(dir,name)as.data.frame(readr::read_csv(file.path(root,dir,name),show_col_types=FALSE))
close<-function(x,y)stopifnot(isTRUE(all.equal(as.numeric(x),as.numeric(y),tolerance=1e-12)))
close_p<-function(x,y){x<-as.numeric(x);y<-as.numeric(y);ok<-is.finite(x)&is.finite(y);stopifnot(identical(is.na(x),is.na(y)),all((x[ok]==0)==(y[ok]==0)));positive<-ok&x>0&y>0;stopifnot(all(abs(log(x[positive])-log(y[positive]))<1e-9))}
gene<-function(dataset,name,value){
 p<-file.path("data/processed/depmap/26Q1",paste0(dataset,".parquet"));schema<-arrow::ParquetFileReader$create(p)$GetSchema()$names
 col<-schema[sub("\\s*\\([^)]*\\)\\s*$","",schema)==name];stopifnot(length(col)==1)
 d<-as.data.frame(arrow::read_parquet(p,col_select=all_of(c("ModelID",col))));names(d)<-c("ModelID",value);d
}
for(case in c("VPS4B_VPS4A","ENO1_ENO2_Method_Test")){
 root<-file.path("results",case);a<-if(grepl("VPS",case))"VPS4B" else "ENO1";b<-if(grepl("VPS",case))"VPS4A" else "ENO2"
 cn<-gene("OmicsCNGeneWGS",a,"cn");dep<-gene("CRISPRGeneEffect",b,"dep");d<-merge(cn,dep,by="ModelID");d$x<-log2(d$cn+1);d<-d[is.finite(d$x)&is.finite(d$dep),]
 t<-read("10_CN_Threshold_Sensitivity","CN_Threshold_Sensitivity.csv")
 for(i in seq_len(nrow(t))){low<-d$x<t$Threshold[i];stopifnot(t$N_low[i]==sum(low),t$N_nonlow[i]==sum(!low));close(t$Delta_median[i],median(d$dep[low])-median(d$dep[!low]));if(t$eligible[i])close_p(t$Wilcoxon_P[i],wilcox.test(d$dep[low],d$dep[!low],exact=FALSE)$p.value) else stopifnot(is.na(t$Wilcoxon_P[i]))}
 close(t$Wilcoxon_FDR,p.adjust(t$Wilcoxon_P,"BH"))
 continuous<-read("10_CN_Threshold_Sensitivity","CN_Threshold_Continuous_Statistics.csv");close(continuous$Pearson_r,cor(d$x,d$dep));close(continuous$Spearman_rho,cor(d$x,d$dep,method="spearman"))
 ex<-gene("ExpressionProteinCoding",a,"expression");d<-merge(ex,dep,by="ModelID");d<-d[is.finite(d$expression)&is.finite(d$dep),]
 r<-read("11_Expression_Dependency","Expression_Dependency.csv");close(r$expression_cutoff,quantile(d$expression,c(.1,.2),type=7,names=FALSE))
 for(i in 1:2){low<-d$expression<=r$expression_cutoff[i];stopifnot(r$N_low[i]==sum(low),r$N_nonlow[i]==sum(!low));close(r$Delta_median[i],median(d$dep[low])-median(d$dep[!low]));close_p(r$Wilcoxon_P[i],wilcox.test(d$dep[low],d$dep[!low],exact=FALSE)$p.value)}
 close(r$Wilcoxon_FDR,p.adjust(r$Wilcoxon_P,"BH"))
 cat("PASS:",case,"actual Parquet-aligned threshold and expression-defined statistics.\n")
}
root<-"results/VPS4B_VPS4A"
pairs<-read("08_TCGA","TCGA_CN_Expression_Samples.csv");by<-read("08_TCGA","TCGA_CN_Expression_ByCancer.csv")
stopifnot(!anyDuplicated(pairs$SampleID),all(by$GDCRelease=="DR46"))
for(i in seq_len(nrow(by))){
 d<-pairs[pairs$CancerType_CN==by$CancerType[i],];stopifnot(nrow(d)==by$N[i])
 if(!by$eligible[i])next
 pe<-cor.test(d$Value_CN,d$Value_RNA);sp<-cor.test(d$Value_CN,d$Value_RNA,method="spearman",exact=FALSE)
 close(by$Pearson_r[i],pe$estimate);close_p(by$Pearson_P[i],pe$p.value);close(by$Spearman_rho[i],sp$estimate);close_p(by$Spearman_P[i],sp$p.value)
 close(c(by$Pearson_CI_low[i],by$Pearson_CI_high[i]),tanh(atanh(pe$estimate)+c(-1,1)*qnorm(.975)/sqrt(nrow(d)-3)))
}
close(by$Pearson_FDR,p.adjust(by$Pearson_P,"BH"));close(by$Spearman_FDR,p.adjust(by$Spearman_P,"BH"))
d<-data.frame(CN=pairs$Value_CN,Expression=pairs$Value_RNA,CancerType=pairs$CancerType_CN)
for(standardized in c(FALSE,TRUE)){
 name<-if(standardized)"TCGA_CN_Expression_Adjusted_Standardized.csv" else "TCGA_CN_Expression_Adjusted.csv"
 r<-read("08_TCGA",name);fit<-if(standardized)lm(scale(Expression)~scale(CN)+CancerType,data=d) else lm(Expression~CN+CancerType,data=d)
 tab<-summary(fit)$coefficients;close(r$estimate,tab[,1]);close(r$std.error,tab[,2]);close_p(r$p.value,tab[,4])
}
screen<-read("12_GenomeWide_Adjusted_Dependency","GenomeWide_Adjusted_Dependency.csv")
stopifnot(nrow(screen)==18531,!anyDuplicated(screen$Gene))
close(screen$FDR_CN,p.adjust(screen$P_CN,"BH"));close(screen$FDR_CNLow,p.adjust(screen$P_CNLow,"BH"))
eligible<-which(!is.na(screen$FDR_CN));ordered<-eligible[order(screen$FDR_CN[eligible],screen$Beta_CN[eligible],screen$Gene[eligible])]
stopifnot(all(screen$Eligible_Rank_adjusted[ordered]==seq_along(ordered)))
cn<-gene("OmicsCNGeneWGS","VPS4B","cn");model<-as.data.frame(arrow::read_parquet("data/processed/depmap/26Q1/Model.parquet"))
meta<-merge(cn,model[,c("ModelID","OncotreeLineage")],by="ModelID");meta$CN_log<-log2(meta$cn+1)
for(g in c("VPS4A","ENO2","RPL5")){
 d<-merge(meta,gene("CRISPRGeneEffect",g,"Chronos"),by="ModelID");d<-d[is.finite(d$Chronos)&is.finite(d$CN_log)&!is.na(d$OncotreeLineage)&nzchar(d$OncotreeLineage),]
 r<-screen[screen$Gene==g,];stopifnot(nrow(d)==r$N)
 tab<-summary(lm(Chronos~CN_log+OncotreeLineage,data=d))$coefficients["CN_log",]
 close(c(r$Beta_CN,r$SE_CN,r$T_CN),tab[1:3]);close_p(r$P_CN,tab[4])
 tab<-summary(lm(Chronos~I(CN_log<.585)+OncotreeLineage,data=d))$coefficients["I(CN_log < 0.585)TRUE",]
 close(c(r$Beta_CNLow,r$SE_CNLow,r$T_CNLow),tab[1:3]);close_p(r$P_CNLow,tab[4])
}
cat("PASS: actual TCGA per-cancer/BH/CI/regression and full adjusted families/candidate/source numeric checks.\n")
jsonlite::write_json(list(checked_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),R_version=as.character(getRversion()),all_passed=TRUE,
 checks=c("VPS4B/VPS4A processed ModelID alignment and threshold/expression statistics","ENO1/ENO2 processed ModelID alignment and threshold/expression statistics","All current DR46 cancer correlations, BH and Fisher-z CI","Raw/standardized cancer-adjusted coefficients and P","Full adjusted continuous/binary BH families and eligible rank","VPS4A, ENO2 and RPL5 independent processed-source adjusted coefficients and P")),
 "data/manifests/extension_results_validation.json",pretty=TRUE,auto_unbox=TRUE)

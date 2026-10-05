# Final workflow statistics. Thresholds and the screen's numeric core are fixed.
CNA_STATES <- c("Deep Deletion","Shallow Deletion","Diploid","Gain","Amplification")
CNA_COLORS <- setNames(c("#2166ac","#67a9cf","#cccccc","#ef8a62","#b2182b"),CNA_STATES)
TCGA_MIN_N <- 20L
tcga_relative_cn <- function(dat) {
 d<-copy(as.data.table(dat))
 if(!all(c("CopyNumber","BaselineCN","CNAState") %in% names(d)))
  stop("TCGA relative CN requires the existing CopyNumber, BaselineCN and CNAState fields")
 if(any(!is.finite(d$BaselineCN)|d$BaselineCN<=0))stop("Invalid existing TCGA sample baseline")
 d[,TCGA_Relative_CN_Change:=CopyNumber/BaselineCN-1]
 d[!is.finite(CopyNumber),TCGA_Relative_CN_Change:=NA_real_]
 d[,Relative_CN_Change:=TCGA_Relative_CN_Change]
 z<-d[is.finite(CopyNumber)]
 valid<-with(z,(CNAState==-2 & CopyNumber==0 & Relative_CN_Change==-1) |
  (CNAState==-1 & CopyNumber>0 & CopyNumber<BaselineCN & Relative_CN_Change> -1 & Relative_CN_Change<0) |
  (CNAState==0 & CopyNumber==BaselineCN & Relative_CN_Change==0) |
  (CNAState==1 & CopyNumber>BaselineCN & CopyNumber<2*BaselineCN & Relative_CN_Change>0 & Relative_CN_Change<1) |
  (CNAState==2 & CopyNumber>=2*BaselineCN & Relative_CN_Change>=1))
 if(anyNA(valid)||!all(valid)||any(is.finite(d$CNAState)&!is.finite(d$CopyNumber)))
  stop("TCGA CNAState and relative copy-number change are inconsistent")
 d
}
workflow_modules <- function(workflow) {
 switch(workflow,geneA_screen=c("tcga","genomewide_dependency","cn_covariation"),
 geneA_geneB=c("tcga","targeted_dependency","lineage_dependency","adjusted_dependency",
              "cn_threshold_sensitivity","genomewide_dependency","cn_covariation"),
 stop("Unsupported workflow: ",workflow))
}
cancer_order <- function(dat,cancers) {
 z<-as.data.table(dat)[is.finite(CopyNumber),.(N=.N,Median_CN=median(CopyNumber),
  Median_Absolute_CN=median(CopyNumber),Median_Relative_CN_Change=median(TCGA_Relative_CN_Change)),by=CancerType]
 z<-merge(data.table(CancerType=cancers),z,by="CancerType",all.x=TRUE)
 z[is.na(N),N:=0L];setorder(z,Median_Relative_CN_Change,CancerType,na.last=TRUE)
 z[,Order:=.I];setcolorder(z,c("Order","CancerType","N","Median_CN"));z
}
ordered_cancer <- function(x,order)factor(x,levels=rev(order$CancerType))
current_cancer_statistics <- function(dat,order,min_n=TCGA_MIN_N) {
 res<-rbindlist(lapply(order$CancerType,function(cancer) {
  d<-dat[CancerType==cancer & is.finite(CopyNumber)&is.finite(Expression)&is.finite(CNAState)]
  r<-as.data.table(pair_correlations(d$TCGA_Relative_CN_Change,d$Expression,min_n))
  r[,`:=`(CN_metric="Relative_CN_Change",Median_Absolute_CN=median(d$CopyNumber),
           Median_Relative_CN_Change=median(d$TCGA_Relative_CN_Change))]
  r[,`:=`(CancerType=cancer,Status=if(.N&&N<min_n)"Insufficient N for correlation" else
           if(!is.finite(Pearson_r))"Constant CN or mRNA; correlation unavailable" else "Eligible")]
  r
 }))
 res[,`:=`(Pearson_FDR=p.adjust(Pearson_P,"BH"),Spearman_FDR=p.adjust(Spearman_P,"BH"))]
 setcolorder(res,c("CancerType",setdiff(names(res),"CancerType")));res
}
current_prevalence <- function(dat,order) {
 d<-copy(dat[is.finite(CNAState)])
 if(!all(d$CNAState %in% -2:2))stop("Invalid five-state CNAState code")
 d[,CNA:=CNA_STATES[as.integer(CNAState)+3L]]
 counts<-d[,.(N=.N),by=.(CancerType,CNA)]
 z<-merge(CJ(CancerType=order$CancerType,CNA=CNA_STATES),counts,all.x=TRUE,
          by=c("CancerType","CNA"))
 z[is.na(N),N:=0L];z[,denominator:=sum(N),by=CancerType]
 z[,Percentage:=ifelse(denominator>0,100*N/denominator,NA_real_)]
 z[,CNA:=factor(CNA,levels=CNA_STATES)]
 z[,Order:=match(CancerType,order$CancerType)];setorder(z,Order,CNA);z
}
targeted_statistics <- function(pair,min_n=MIN_N) {
 # Supplied core reports continuous CN_log statistics; the new scatter uses
 # relative CN, so both scales are explicitly reported rather than conflated.
 core<-pair_correlations(pair$CN_log,pair$Chronos,min_n=3L)
 relative<-pair_correlations(pair$CN_relative,pair$Chronos,min_n=3L)
 effect<-group_effect(pair$Chronos,pair$CN_binary=="CN-Low",min_n)
 cbind(data.table(geneA=GENE_A,geneB=GENE_B),as.data.table(core),
       as.data.table(relative)[,.(Relative_CN_Pearson_r=Pearson_r,Relative_CN_Pearson_P=Pearson_P,
                               Relative_CN_Spearman_rho=Spearman_rho,Relative_CN_Spearman_P=Spearman_P)],
       as.data.table(effect))
}
eligible_rank <- function(res) {
 res<-copy(res);res[,Eligible_Rank:=NA_integer_]
 ok<-which(is.finite(res$Wilcoxon_FDR)&is.finite(res$Delta_median))
 ix<-ok[order(res$Wilcoxon_FDR[ok],res$Delta_median[ok],res$Rank[ok])]
 res$Eligible_Rank[ix]<-seq_along(ix);res[,Eligible_N:=length(ix)];res
}
dependency_top <- function(res,n=20L) {
 a<-res[is.finite(Eligible_Rank)&Delta_median<0];setorder(a,Eligible_Rank);a<-head(a,n)
 b<-res[is.finite(Eligible_Rank)&Delta_median>0];setorder(b,Eligible_Rank);b<-head(b,n)
 a[,Direction:="Stronger dependency in CN-Low"];b[,Direction:="Weaker dependency in CN-Low"]
 rbindlist(list(a,b))
}
screen_core_statistics <- function() {
 # Execute only the supplied statistical prefix once. The plotting/export tail
 # is replaced by the final output specification; every numeric expression stays intact.
 expressions<-as.list(body(run_genomewide))[-1L]
 end<-which(vapply(expressions,function(e)identical(e,quote(result[,Rank := .I])),logical(1)))
 stopifnot(length(end)==1L)
 f<-run_genomewide
 body(f)<-as.call(c(list(as.name("{")),expressions[3L:end],list(quote(result))))
 eligible_rank(f())
}
screen_group_counts <- function() {
 # Read only the target CN column and Chronos IDs before deciding whether the
 # group comparison is possible. Do not change the supplied screen's statistics.
 cn<-as.data.table(prepare_cn(read_gene(FILES$cn,GENE_A,"CN_relative")))
 id<-identify_id_column(names(fread(FILES$chronos,nrows=0)))
 ids<-as.character(fread(FILES$chronos,select=id)[[id]])
 matched<-cn[match(ids,ModelID)][is.finite(CN_log)]
 list(N_low=sum(matched$CN_binary=="CN-Low"),N_nonlow=sum(matched$CN_binary=="CN-NonLow"))
}
empty_screen_statistics <- function() {
 data.table(Gene=character(),N=integer(),Pearson_r=numeric(),Pearson_P=numeric(),
  N_low=integer(),N_nonlow=integer(),Mean_low=numeric(),Mean_nonlow=numeric(),
  Median_low=numeric(),Median_nonlow=numeric(),Delta_mean=numeric(),Delta_median=numeric(),
  Wilcoxon_P=numeric(),Pearson_FDR=numeric(),Wilcoxon_FDR=numeric(),Rank=integer(),
  Eligible_Rank=integer(),Eligible_N=integer())
}
covariation_statistics <- function(dt,gene) {
 id<-identify_id_column(names(dt));cols<-setdiff(names(dt),id)
 a<-find_gene_col(cols,gene);x<-as.numeric(dt[[a]])
 res<-rbindlist(lapply(cols,function(g) {
  y<-as.numeric(dt[[g]]);ok<-is.finite(x)&is.finite(y);n<-sum(ok)
  r<-if(n>=10)suppressWarnings(cor(x[ok],y[ok])) else NA_real_
  p<-if(is.finite(r)&&abs(r)<1)2*pt(-abs(r*sqrt((n-2)/(1-r^2))),df=n-2) else NA_real_
  data.table(Gene=clean_gene(g),N=n,Pearson_r=r,P=p)
 }))
 res[,FDR:=p.adjust(P,"BH")];setorder(res,-Pearson_r);res
}
covariation_top <- function(res,gene,n=20L) {
 a<-head(res[Gene!=gene & is.finite(Pearson_r)&Pearson_r>0][order(-Pearson_r,Gene)],n)
 b<-head(res[Gene!=gene & is.finite(Pearson_r)&Pearson_r<0][order(Pearson_r,Gene)],n)
 a[,Direction:="Positive CN correlation"];b[,Direction:="Negative CN correlation"]
 rbindlist(list(a,b))
}
waterfall_order <- function(pair) {
 d<-as.data.table(pair);d<-copy(d);setorder(d,-Chronos,ModelID);d[,Waterfall_Order:=.I];d
}

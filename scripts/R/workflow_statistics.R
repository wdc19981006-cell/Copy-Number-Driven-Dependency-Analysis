# Final workflow statistics. Thresholds and the screen's numeric core are fixed.
CNA_STATES <- c("Deep Deletion","Shallow Deletion","Diploid","Gain","Amplification")
CNA_COLORS <- setNames(c("#2166ac","#67a9cf","#cccccc","#ef8a62","#b2182b"),CNA_STATES)
TCGA_MIN_N <- 20L
workflow_modules <- function(workflow) {
 switch(workflow,geneA_screen=c("tcga","genomewide_dependency","cn_covariation"),
 geneA_geneB=c("tcga","targeted_dependency","lineage_dependency","adjusted_dependency",
              "cn_threshold_sensitivity","genomewide_dependency","cn_covariation"),
 stop("Unsupported workflow: ",workflow))
}
cancer_order <- function(dat,cancers) {
 z<-as.data.table(dat)[is.finite(CopyNumber),.(N=.N,Median_CN=median(CopyNumber)),by=CancerType]
 z<-merge(data.table(CancerType=cancers),z,by="CancerType",all.x=TRUE)
 z[is.na(N),N:=0L];setorder(z,Median_CN,CancerType,na.last=TRUE)
 z[,Order:=.I];setcolorder(z,c("Order","CancerType","N","Median_CN"));z
}
ordered_cancer <- function(x,order)factor(x,levels=rev(order$CancerType))
reference_cancer_statistics <- function(dat,order,min_n=TCGA_MIN_N) {
 res<-rbindlist(lapply(order$CancerType,function(cancer) {
  d<-dat[CancerType==cancer & is.finite(CopyNumber)&is.finite(Expression)&is.finite(GISTIC)]
  r<-as.data.table(pair_correlations(d$CopyNumber,d$Expression,min_n))
  r[,`:=`(CancerType=cancer,Status=if(.N&&N<min_n)"Insufficient N for correlation" else
           if(!is.finite(Pearson_r))"Constant CN or mRNA; correlation unavailable" else "Eligible")]
  r
 }))
 res[,`:=`(Pearson_FDR=p.adjust(Pearson_P,"BH"),Spearman_FDR=p.adjust(Spearman_P,"BH"))]
 setcolorder(res,c("CancerType",setdiff(names(res),"CancerType")));res
}
reference_prevalence <- function(dat,order) {
 d<-copy(dat[is.finite(GISTIC)])
 if(!all(d$GISTIC %in% -2:2))stop("Invalid five-state GISTIC code")
 d[,CNA:=CNA_STATES[as.integer(GISTIC)+3L]]
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

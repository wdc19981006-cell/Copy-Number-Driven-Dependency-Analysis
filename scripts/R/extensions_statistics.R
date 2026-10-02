# Pure statistical extensions. No files, downloads or original core functions.
finite_pair <- function(x,y) {ok<-is.finite(x)&is.finite(y);list(x=x[ok],y=y[ok])}
pair_correlations <- function(x,y,min_n=3L) {
 d<-finite_pair(x,y);x<-d$x;y<-d$y
 out<-data.frame(N=length(x),Pearson_r=NA_real_,Pearson_P=NA_real_,Spearman_rho=NA_real_,Spearman_P=NA_real_)
 if(length(x)<min_n||sd(x)==0||sd(y)==0)return(out)
 pe<-cor.test(x,y);sp<-cor.test(x,y,method="spearman",exact=FALSE)
 out[1,2:5]<-c(unname(pe$estimate),pe$p.value,unname(sp$estimate),sp$p.value);out
}
group_effect <- function(y,low,min_n=3L) {
 stopifnot(length(y)==length(low));ok<-is.finite(y)&!is.na(low)
 a<-y[ok&!is.na(low)&low];b<-y[ok&!is.na(low)&!low]
 eligible<-length(a)>=min_n&&length(b)>=min_n
 safe<-function(f,x)if(length(x))f(x) else NA_real_
 data.frame(N_low=length(a),N_nonlow=length(b),Mean_low=safe(mean,a),Mean_nonlow=safe(mean,b),Median_low=safe(median,a),Median_nonlow=safe(median,b),
 Delta_mean=safe(mean,a)-safe(mean,b),Delta_median=safe(median,a)-safe(median,b),Wilcoxon_P=if(eligible)wilcox.test(a,b,exact=FALSE)$p.value else NA_real_,eligible=eligible)
}
cn_threshold_statistics <- function(cn_log,chronos,thresholds=c(.585,.50,.40,.35),min_n=3L) {
 d<-finite_pair(cn_log,chronos)
 out<-do.call(rbind,lapply(thresholds,function(t)cbind(Threshold=t,group_effect(d$y,d$x<t,min_n))))
 out$Wilcoxon_FDR<-p.adjust(out$Wilcoxon_P,"BH");out
}
expression_dependency_statistics <- function(expression,chronos,quantiles=c(.10,.20),min_n=3L,cutoffs=NULL) {
 d<-finite_pair(expression,chronos)
 if(is.null(cutoffs))cutoffs<-if(length(d$x))quantile(d$x,quantiles,type=7,names=FALSE) else rep(NA_real_,length(quantiles))
 out<-do.call(rbind,lapply(seq_along(quantiles),function(i)cbind(threshold_quantile=quantiles[i],expression_cutoff=cutoffs[i],group_effect(d$y,d$x<=cutoffs[i],min_n))))
 out$Wilcoxon_FDR<-p.adjust(out$Wilcoxon_P,"BH")
 list(continuous=pair_correlations(d$x,d$y),groups=out)
}
coefficient_table <- function(fit) {
 x<-as.data.frame(summary(fit)$coefficients);names(x)<-c("estimate","std.error","statistic","p.value")
 cbind(term=rownames(x),x,row.names=NULL)
}
tcga_cancer_statistics <- function(dat,min_n=20L) {
 d<-dat[is.finite(dat$CN)&is.finite(dat$Expression)&!is.na(dat$CancerType)&nzchar(dat$CancerType),,drop=FALSE]
 cancers<-sort(unique(as.character(d$CancerType)))
 if(!length(cancers))return(data.frame())
 out<-do.call(rbind,lapply(cancers,function(cancer){
  a<-d[d$CancerType==cancer,,drop=FALSE];r<-pair_correlations(a$CN,a$Expression,min_n)
  eligible<-nrow(a)>=min_n&&is.finite(r$Pearson_r)
  ci<-c(NA_real_,NA_real_)
  if(eligible)ci<-if(abs(r$Pearson_r)>=1)rep(sign(r$Pearson_r),2) else tanh(atanh(r$Pearson_r)+c(-1,1)*qnorm(.975)/sqrt(r$N-3))
  cbind(CancerType=cancer,r,CN_median=median(a$CN),Expression_median=median(a$Expression),Pearson_CI_low=ci[1],Pearson_CI_high=ci[2],eligible=eligible)
 }))
 out$Pearson_FDR<-p.adjust(out$Pearson_P,"BH");out$Spearman_FDR<-p.adjust(out$Spearman_P,"BH");out
}
tcga_adjusted_statistics <- function(dat,standardized=FALSE) {
 d<-dat[is.finite(dat$CN)&is.finite(dat$Expression)&!is.na(dat$CancerType)&nzchar(dat$CancerType),,drop=FALSE]
 d$CancerType<-droplevels(factor(d$CancerType))
 if(nrow(d)<20||nlevels(d$CancerType)<2||sd(d$CN)==0||sd(d$Expression)==0)return(list(eligible=FALSE,reason="Need N >= 20, two cancer levels and variable CN/expression",coefficients=data.frame()))
 if(standardized){d$CN<-as.numeric(scale(d$CN));d$Expression<-as.numeric(scale(d$Expression))}
 fit<-lm(Expression~CN+CancerType,data=d)
 if(fit$rank<ncol(fit$qr$qr))return(list(eligible=FALSE,reason="CN is aliased with cancer effects",coefficients=data.frame()))
 list(eligible=TRUE,N=nrow(d),coefficients=coefficient_table(fit),fit=fit)
}
adjusted_gene <- function(y,cn_log,lineage,min_n=20L) {
 ok<-is.finite(y)&is.finite(cn_log)&!is.na(lineage)&nzchar(as.character(lineage))
 d<-data.frame(Chronos=y[ok],CN_log=cn_log[ok],Lineage=droplevels(factor(lineage[ok])))
 out<-data.frame(N=nrow(d),Beta_CN=NA_real_,SE_CN=NA_real_,T_CN=NA_real_,P_CN=NA_real_,Beta_CNLow=NA_real_,SE_CNLow=NA_real_,T_CNLow=NA_real_,P_CNLow=NA_real_,continuous_status="skipped",binary_status="skipped",skip_reason="")
 if(nrow(d)<min_n||nlevels(d$Lineage)<2){out$skip_reason<-"Need N >= 20 and at least two lineage levels";return(out)}
 if(sd(d$Chronos)==0){out$skip_reason<-"Constant Chronos outcome";return(out)}
 d$CNLow<-d$CN_log<.585
 for(binary in c(FALSE,TRUE)){
  fit<-if(binary)lm(Chronos~CNLow+Lineage,data=d) else lm(Chronos~CN_log+Lineage,data=d)
  if(fit$rank<ncol(fit$qr$qr)){out$skip_reason<-"CN term aliased with lineage effects";next}
  tab<-summary(fit)$coefficients;term<-if(binary)"CNLowTRUE" else "CN_log"
  if(term %in% rownames(tab)&&all(is.finite(tab[term,]))){
   cols<-if(binary)c("Beta_CNLow","SE_CNLow","T_CNLow","P_CNLow") else c("Beta_CN","SE_CN","T_CN","P_CN")
   out[1,cols]<-as.numeric(tab[term,]);out[[if(binary)"binary_status" else "continuous_status"]]<-"completed"
  }else out$skip_reason<-"CN term not estimable after lineage adjustment"
 }
 out
}
align_screen <- function(matrix,metadata) {
 stopifnot("ModelID" %in% names(matrix),"ModelID" %in% names(metadata),!anyDuplicated(matrix$ModelID),!anyDuplicated(metadata$ModelID))
 indices<-match(as.character(matrix$ModelID),as.character(metadata$ModelID));keep<-which(!is.na(indices))
 list(matrix=matrix[keep,,drop=FALSE],metadata=metadata[indices[keep],,drop=FALSE])
}
genomewide_adjusted_statistics <- function(matrix,metadata,min_n=20L) {
 a<-align_screen(matrix,metadata);genes<-setdiff(names(a$matrix),"ModelID")
 out<-do.call(rbind,lapply(genes,function(g)cbind(Gene=g,adjusted_gene(as.numeric(a$matrix[[g]]),a$metadata$CN_log,a$metadata$OncotreeLineage,min_n))))
 out$FDR_CN<-p.adjust(out$P_CN,"BH");out$FDR_CNLow<-p.adjust(out$P_CNLow,"BH")
 out$Eligible_Rank_adjusted<-NA_integer_;ok<-which(!is.na(out$FDR_CN));order<-ok[order(out$FDR_CN[ok],out$Beta_CN[ok],out$Gene[ok])]
 out$Eligible_Rank_adjusted[order]<-seq_along(order);out$Eligible_N_adjusted<-length(order);out
}
genomewide_expression_statistics <- function(matrix,metadata) {
 a<-align_screen(matrix,metadata);good<-is.finite(a$metadata$Expression)
 a$matrix<-a$matrix[good,,drop=FALSE];x<-a$metadata$Expression[good]
 cutoffs<-if(length(x))quantile(x,c(.1,.2),type=7,names=FALSE) else c(NA_real_,NA_real_)
 genes<-setdiff(names(a$matrix),"ModelID");continuous<-list();groups<-list()
 for(i in seq_along(genes)){
  r<-expression_dependency_statistics(x,as.numeric(a$matrix[[genes[i]]]),cutoffs=cutoffs)
  continuous[[i]]<-cbind(Gene=genes[i],r$continuous)
  groups[[i]]<-cbind(Gene=genes[i],r$groups)
 }
 continuous<-do.call(rbind,continuous);continuous$Pearson_FDR<-p.adjust(continuous$Pearson_P,"BH")
 groups<-do.call(rbind,groups)
 for(q in c(.1,.2)){ok<-groups$threshold_quantile==q;groups$Wilcoxon_FDR[ok]<-p.adjust(groups$Wilcoxon_P[ok],"BH")}
 merge(groups,continuous,by="Gene",sort=FALSE)
}

source("tests/extensions_fixture.R")
set.seed(942)
n<-90;cn<-runif(n);lineage<-rep(c("A","B","C"),each=30)
meta<-data.frame(ModelID=paste0("m",seq_len(n)),CN_log=cn,OncotreeLineage=lineage)
mat<-data.frame(ModelID=meta$ModelID,G1=.7*cn+match(lineage,LETTERS)+rnorm(n),G2=-.4*(cn<.585)+match(lineage,LETTERS)+rnorm(n),G3=rep(NA_real_,n))
mat$G1[c(5,12)]<-NA;mat<-mat[rev(seq_len(n)),];meta<-meta[c(31:90,1:30),]
r<-genomewide_adjusted_statistics(mat,meta)
aligned<-meta[match(mat$ModelID,meta$ModelID),]
for(g in c("G1","G2")){
 d<-data.frame(y=mat[[g]],x=aligned$CN_log,lineage=aligned$OncotreeLineage);a<-r[r$Gene==g,]
 fit<-lm(y~x+lineage,data=d);tab<-summary(fit)$coefficients["x",]
 assert_close(c(a$Beta_CN,a$SE_CN,a$T_CN,a$P_CN),tab)
 binary<-summary(lm(y~I(x<.585)+lineage,data=d))$coefficients["I(x < 0.585)TRUE",]
 assert_close(c(a$Beta_CNLow,a$SE_CNLow,a$T_CNLow,a$P_CNLow),binary)
 stopifnot(a$N==sum(complete.cases(d)))
}
assert_close(r$FDR_CN,p.adjust(r$P_CN,"BH"));assert_close(r$FDR_CNLow,p.adjust(r$P_CNLow,"BH"))
stopifnot(r$continuous_status[r$Gene=="G3"]=="skipped",all(r$Eligible_N_adjusted==2))
stopifnot(adjusted_gene(1:19,1:19,rep(c("A","B"),length.out=19))$continuous_status=="skipped")
stopifnot(adjusted_gene(1:30,1:30,rep("A",30))$continuous_status=="skipped")
stopifnot(inherits(try(align_screen(mat,rbind(meta,meta[1,])),silent=TRUE),"try-error"))
aliased<-adjusted_gene(seq_len(30),rep(c(.1,.8),each=15),rep(c("A","B"),each=15))
stopifnot(aliased$continuous_status=="skipped",aliased$binary_status=="skipped",is.na(aliased$P_CN))
# Across-lineage confounding disappears after adjustment, matching direct lm.
cat("PASS: beta/SE/t/P, separate BH, exact shuffled-ID alignment, lineage adjustment, N/levels skips and duplicate rejection.\n")

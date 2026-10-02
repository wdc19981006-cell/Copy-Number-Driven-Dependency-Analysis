source("tests/extensions_fixture.R")
x<-c(.2,.34,.35,.39,.40,.49,.50,.584,.585,.7,.8,.9,NA,Inf)
y<-c(-2,-1.8,-1.5,-1.4,-1,-.9,-.8,-.5,-.3,0,.1,.2,-10,-10)
r<-cn_threshold_statistics(x,y)
for(i in seq_len(nrow(r))){
 ok<-is.finite(x)&is.finite(y);low<-ok&x<r$Threshold[i];non<-ok&x>=r$Threshold[i]
 stopifnot(r$N_low[i]==sum(low),r$N_nonlow[i]==sum(non),r$eligible[i]==(sum(low)>=3&&sum(non)>=3))
 assert_close(r$Delta_median[i],median(y[low])-median(y[non]))
 if(r$eligible[i])assert_p(r$Wilcoxon_P[i],wilcox.test(y[low],y[non],exact=FALSE)$p.value) else stopifnot(is.na(r$Wilcoxon_P[i]))
}
stopifnot(r$N_low[r$Threshold==.35]==2,!r$eligible[r$Threshold==.35])
assert_close(r$Wilcoxon_FDR,p.adjust(r$Wilcoxon_P,"BH"))
stopifnot(all(!cn_threshold_statistics(c(NA,Inf),c(1,2))$eligible))
cat("PASS: threshold boundary grouping, counts, delta, Wilcoxon, BH and insufficient-N skips.\n")

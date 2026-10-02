source("tests/extensions_fixture.R")
x<-c(1:100,NA,Inf,2);y<-c(sin(1:100/8)-(1:100)/100,-50,-50,NA)
r<-expression_dependency_statistics(x,y);ok<-is.finite(x)&is.finite(y)
assert_close(r$groups$expression_cutoff,quantile(x[ok],c(.1,.2),type=7,names=FALSE))
stopifnot(identical(r$groups$N_low,c(10L,20L)),all(r$groups$N_low+r$groups$N_nonlow==100),r$continuous$N==100)
assert_close(r$continuous$Pearson_r,cor(x[ok],y[ok]));assert_close(r$continuous$Spearman_rho,cor(x[ok],y[ok],method="spearman"))
for(i in 1:2){low<-x[ok]<=r$groups$expression_cutoff[i];assert_p(r$groups$Wilcoxon_P[i],wilcox.test(y[ok][low],y[ok][!low],exact=FALSE)$p.value)}
assert_close(r$groups$Wilcoxon_FDR,p.adjust(r$groups$Wilcoxon_P,"BH"))
stopifnot(!expression_dependency_statistics(1:10,1:10)$groups$eligible[1])
ties<-expression_dependency_statistics(c(rep(1,30),2:71),seq_len(100))
stopifnot(all(ties$groups$expression_cutoff==1),all(ties$groups$N_low==30))
# Shared expression cutoffs, explicit ModelID alignment and separate BH families.
mat<-data.frame(ModelID=paste0("m",100:1),A=rev(y[1:100]),B=rev(cos(1:100)))
meta<-data.frame(ModelID=paste0("m",1:100),Expression=1:100)
screen<-genomewide_expression_statistics(mat,meta)
for(q in c(.1,.2)){a<-screen[screen$threshold_quantile==q,];assert_close(a$Wilcoxon_FDR,p.adjust(a$Wilcoxon_P,"BH"));assert_close(a$expression_cutoff,rep(quantile(1:100,q),2))}
assert_close(screen$Pearson_r[screen$Gene=="A"],rep(cor(1:100,y[1:100]),2))
cat("PASS: 10/20% type-7 quantiles, ties/NA exclusion, continuous and grouped tests, genome-wide alignment/BH.\n")

Sys.setlocale("LC_CTYPE","English_United States.utf8")
source(".Rprofile")
stopifnot(as.character(getRversion())=="4.5.0")
folder<-"results/VPS4B_VPS4A/08_TCGA"
dat<-readr::read_csv(file.path(folder,"TCGA_CN_Expression_Samples.csv"),show_col_types=FALSE)
stat<-readr::read_csv(file.path(folder,"TCGA_CN_Expression_Statistics.csv"),show_col_types=FALSE)
pe<-cor.test(dat$Value_CN,dat$Value_RNA)
sp<-cor.test(dat$Value_CN,dat$Value_RNA,method="spearman",exact=FALSE)
close<-function(x,y)isTRUE(all.equal(as.numeric(x),as.numeric(y),tolerance=1e-12))
# Absolute tolerance alone would incorrectly accept zero for these tiny P values.
close_p<-function(x,y)all(is.finite(c(x,y))) && x>0 && y>0 && abs(log(x)-log(y))<1e-10
stopifnot(!anyDuplicated(dat$SampleID),stat$N==nrow(dat),close(stat$Pearson_r,pe$estimate),close_p(stat$Pearson_P,pe$p.value),close(stat$Spearman_rho,sp$estimate),close_p(stat$Spearman_P,sp$p.value))
cat("PASS: actual matched sample CSV independently reproduces Pearson/Spearman and P values in R 4.5.0.\n")

Sys.setlocale("LC_CTYPE","English_United States.utf8")
source(".Rprofile")
library(data.table);library(dplyr);library(ggplot2)
source("scripts/R/tcga_analysis.R")
GENE_A<-"TESTGENE"
test_dir<-tempfile("tcga_selection_",tmpdir=".runtime")
dir.create(test_dir)
out_dir<-function(...)test_dir
tcga_ready<-function()TRUE
tcga_gene<-function(dataset){
 if(dataset=="cn")return(data.table(SampleID=paste0("s",1:4),AliquotID=paste0("a",1:4),Value=c(1,2,3,4)))
 data.table(SampleID=c("s1","s1","s2","s2","s3","s4"),AliquotID=c("other","a1","other","other","a3","a4"),FileID=c("0","9","9","0","1","1"),Value=c(99,2,99,3,5,7))
}
save_pdf<-function(...)invisible(NULL)
run_tcga_expression()
audit<-fread(file.path(test_dir,"TCGA_RNA_Representative_Selection.csv"))
selected<-audit[selected_for_sample_analysis==TRUE]
stopifnot(nrow(audit)==6,nrow(selected)==4,selected[SampleID=="s1",FileID]=="9",selected[SampleID=="s2",FileID]=="0")
dat<-fread(file.path(test_dir,"TCGA_CN_Expression_Samples.csv"))
stats<-fread(file.path(test_dir,"TCGA_CN_Expression_Statistics.csv"))
stopifnot(!anyDuplicated(dat$SampleID),identical(dat$Value_RNA,c(2L,3L,5L,7L)),stats$N==4,isTRUE(all.equal(stats$Pearson_r,unname(cor.test(1:4,c(2,3,5,7))$estimate))))
cat("PASS: sample UUID join; prefer matched CN aliquot, then RNA file UUID; alternatives audited; statistics use selected samples. Synthetic check only.\n")

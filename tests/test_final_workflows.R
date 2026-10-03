invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
source(".Rprofile")
library(data.table);library(dplyr);library(ggplot2)
source("scripts/R/data_access.R")
source("scripts/R/dependency_screen.R")
source("scripts/R/lineage_analysis.R")
source("scripts/R/extensions_statistics.R")
source("scripts/R/workflow_statistics.R")
source("scripts/R/workflow_plots.R")
source("scripts/R/final_workflow.R")
GENE_A<-"A";GENE_B<-"B";GENE_B_PROVIDED<-TRUE;MIN_N<-3L
wf_route_supplements()
stopifnot(grepl("Supplementary/01_Lineage",paste(deparse(body(run_lineage)),collapse=" "),fixed=TRUE),
 grepl("Supplementary/02_Adjusted",paste(deparse(body(run_adjusted)),collapse=" "),fixed=TRUE))
stopifnot(identical(workflow_modules("geneA_screen"),c("tcga","genomewide_dependency","cn_covariation")),
 !any(c("reverse_dependency","mutation_dependency","expression_dependency","genomewide_expression_dependency","genomewide_adjusted_dependency") %in% workflow_modules("geneA_geneB")),
 sum(workflow_modules("geneA_geneB")=="genomewide_dependency")==1L)
d<-data.table(SampleID=paste0("s",1:15),CancerType=rep(c("Z","A","C"),each=5),
 CopyNumber=c(-9,-2,-1,0,99,rep(0,5),1:5),Expression=seq_len(15),GISTIC=rep(-2:2,3))
order<-cancer_order(d,c("Z","A","C"));stopifnot(identical(order$CancerType,c("Z","A","C")),all(order$N==5L))
prev<-reference_prevalence(d,order);stopifnot(nrow(prev)==15L,all(prev$N==1L),all(prev$Percentage==20))
p1<-tcga_prevalence_plot(prev,order);p2<-tcga_landscape_plot(d,order)
stopifnot(identical(levels(p1$data$CancerType),rev(order$CancerType)),identical(levels(p1$data$CancerType),levels(p2$data$CancerType)))
b<-ggplot_build(p2);stopifnot(nrow(b$data[[2]])==nrow(d),isTRUE(all.equal(sort(b$data[[2]]$y),sort(d$CopyNumber))))
# Finite paired cohorts only; insufficient cancers still retain a statistics row.
d$Expression[1]<-NA
stats<-reference_cancer_statistics(d,order,min_n=5L)
stopifnot(stats$N[1]==4L,is.na(stats$Pearson_r[1]),stats$Status[1]=="Insufficient N for correlation")
for(i in 2:3) {
 z<-d[CancerType==order$CancerType[i]]
 if(sd(z$CopyNumber)>0)stopifnot(abs(stats$Pearson_r[i]-cor(z$CopyNumber,z$Expression))<1e-10)
}
pair<-as.data.table(prepare_cn(data.frame(ModelID=paste0("m",1:12),CN_relative=c(rep(.2,6),seq(1,2,length.out=6)))))
pair[,Chronos:=c(-2,-1.8,-1.7,-1.6,-1.5,-1.4,-.1,-.2,-.3,-.4,-.5,-.6)]
s<-targeted_statistics(pair)
stopifnot(s$N_low==6L,s$N_nonlow==6L,s$Delta_median<0,
 isTRUE(all.equal(s$Wilcoxon_P,wilcox.test(pair$Chronos[1:6],pair$Chronos[7:12],exact=FALSE)$p.value)))
w<-waterfall_order(pair);stopifnot(all(diff(w$Chronos)<=0),w$ModelID[1]=="m7",w$ModelID[12]=="m1")
# Execute the original screen's numeric prefix using an in-memory synthetic matrix.
core_fread<-fread;core_read_gene<-read_gene;core_msg<-msg
matrix<-data.table(ModelID=pair$ModelID,`B (1)`=pair$Chronos,`C (2)`=-pair$Chronos,`D (3)`=NA_real_)
fread<-function(...)copy(matrix)
read_gene<-function(...)pair[,.(ModelID,CN_relative)]
FILES<-list(cn="unused",chronos="unused")
msg<-function(...)invisible(NULL)
res<-screen_core_statistics()
stopifnot(res[Gene=="B",Delta_median]<0,res[Gene=="C",Delta_median]>0,
 is.na(res[Gene=="D",Eligible_Rank]),res[Gene=="B",Eligible_Rank]==1L,
 isTRUE(all.equal(res$Wilcoxon_FDR,p.adjust(res$Wilcoxon_P,"BH"))))
top<-dependency_top(res,20);stopifnot(nrow(top)==2L,all(top[Direction=="Stronger dependency in CN-Low",Delta_median]<0))
cov<-covariation_statistics(data.table(ModelID=paste0("m",1:12),A=1:12,POS=2*(1:12),NEG=-(1:12)),"A")
ct<-covariation_top(cov,"A");stopifnot(nrow(ct)==2L,ct[Gene=="POS",Pearson_r]==1,ct[Gene=="NEG",Pearson_r]==-1)
fread<-core_fread;read_gene<-core_read_gene;msg<-core_msg
# Cache only accepts complete, structurally valid files with unchanged content hashes.
RESULT_ROOT<-tempfile("wf_cache_",tmpdir=".runtime");dir.create(file.path(RESULT_ROOT,"Main_Results"),recursive=TRUE)
path<-"Main_Results/test.pdf";writeLines("%PDF-1.4 synthetic",file.path(RESULT_ROOT,path))
old<-list(key="same",outputs=wf_hash_outputs(path))
stopifnot(wf_cache_valid(old,"same",path),!wf_cache_valid(old,"changed-parameter-or-code",path))
jsonlite::write_json(old,file.path(RESULT_ROOT,"cache.json"),auto_unbox=TRUE)
roundtrip<-jsonlite::fromJSON(file.path(RESULT_ROOT,"cache.json"),simplifyVector=FALSE)
stopifnot(wf_cache_valid(roundtrip,"same",path))
writeLines("%PDF-1.4 modified",file.path(RESULT_ROOT,path));stopifnot(!wf_cache_valid(old,"same",path))
unlink(RESULT_ROOT,recursive=TRUE)
source("scripts/R/local_data.R")
e<-tryCatch(local_required("missing.parquet"),error=conditionMessage)
stopifnot(grepl("Local processed dataset unavailable.",e,fixed=TRUE),grepl("Run the separate data-update/preparation workflow first.",e,fixed=TRUE))
e<-tryCatch(fread("data/raw/forbidden.csv"),error=conditionMessage)
stopifnot(grepl("forbids raw data reads",e,fixed=TRUE))
cat("PASS: workflows, shared visual order, all points, CNA counts, paired statistics, eligible ranks, directions, cache corruption and local-only failures.\n")

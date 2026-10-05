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
# Both baseline scales must give the same five-state relative values.
for(base in c(2,4)) {
 z<-data.table(CopyNumber=base*c(0,.5,1,1.5,2),BaselineCN=base,CNAState=-2:2)
 rel<-tcga_relative_cn(z)
 stopifnot(identical(rel$Relative_CN_Change,c(-1,-.5,0,.5,1)),
  identical(rel$TCGA_Relative_CN_Change,rel$Relative_CN_Change),identical(z$CopyNumber,rel$CopyNumber))
 for(i in 1:5) {
  bad<-copy(z);bad$CNAState[i]<-if(i==5)-2L else z$CNAState[i+1]
  stopifnot(inherits(try(tcga_relative_cn(bad),silent=TRUE),"try-error"))
 }
 for(invalid in c(0,-1,NA_real_,Inf)) {
  bad<-copy(z);bad$BaselineCN[1]<-invalid
  stopifnot(inherits(try(tcga_relative_cn(bad),silent=TRUE),"try-error"))
 }
}
d<-tcga_relative_cn(data.table(SampleID=paste0("s",1:15),CancerType=rep(c("Z","A","C"),each=5),
 CopyNumber=c(0,1,1,2,4,0,2,4,6,8,0:4),BaselineCN=rep(c(2,4,2),each=5),
 Expression=seq_len(15),CNAState=c(-2,-1,-1,0,2,rep(-2:2,2))))
order<-cancer_order(d,c("Z","A","C"));stopifnot(identical(order$CancerType,c("Z","A","C")),all(order$N==5L))
stopifnot(!identical(order$CancerType,order[order(Median_Absolute_CN,CancerType),CancerType]))
prev<-current_prevalence(d,order);stopifnot(nrow(prev)==15L,all(prev$denominator==5L),sum(prev$N)==15L)
p1<-tcga_prevalence_plot(prev,order);p2<-tcga_landscape_plot(d,order)
rounded<-copy(prev);rounded[,Percentage:=Percentage+1e-12]
rounded_build<-ggplot_build(tcga_prevalence_plot(rounded,order))$data[[1]]
stopifnot(all(is.finite(rounded_build$ymin)),all(is.finite(rounded_build$ymax)),max(rounded_build$ymax)<=1)
stopifnot(identical(levels(p1$data$CancerType),rev(order$CancerType)),identical(levels(p1$data$CancerType),levels(p2$data$CancerType)))
b<-ggplot_build(p2);stopifnot(nrow(b$data[[2]])==nrow(d),isTRUE(all.equal(sort(b$data[[2]]$y),sort(d$Relative_CN_Change))),
 identical(b$data[[3]]$yintercept,0),identical(b$data[[3]]$linetype,"dashed"),inherits(p2$coordinates,"CoordFlip"))
limits<-p2$scales$get_scales("y")$limits;stopifnot(limits[1]== -limits[2])
# Finite paired cohorts only; insufficient cancers still retain a statistics row.
d$Expression[1]<-NA
stats<-current_cancer_statistics(d,order,min_n=5L)
stopifnot(stats$N[1]==4L,is.na(stats$Pearson_r[1]),stats$Status[1]=="Insufficient N for correlation")
for(i in 2:3) {
 z<-d[CancerType==order$CancerType[i]]
 if(sd(z$Relative_CN_Change)>0)stopifnot(abs(stats$Pearson_r[i]-cor(z$Relative_CN_Change,z$Expression))<1e-10)
}
# Vary baselines within a cancer so absolute and relative correlations differ.
rna<-data.table(CancerType="A",CopyNumber=rep(c(1,2,3,4),10),BaselineCN=rep(c(2,4),each=20))
rna[,CNAState:=ifelse(CopyNumber<BaselineCN,-1L,ifelse(CopyNumber==BaselineCN,0L,ifelse(CopyNumber<2*BaselineCN,1L,2L)))]
rna<-tcga_relative_cn(rna);rna[,Expression:=TCGA_Relative_CN_Change+rep(c(0,.1),20)]
rs<-current_cancer_statistics(rna,cancer_order(rna,"A"))
pe<-cor.test(rna$Relative_CN_Change,rna$Expression)
sp<-cor.test(rna$Relative_CN_Change,rna$Expression,method="spearman",exact=FALSE)
stopifnot(rs$CN_metric=="Relative_CN_Change",abs(rs$Pearson_r-pe$estimate)<1e-12,
 abs(rs$Spearman_rho-sp$estimate)<1e-12,abs(rs$Pearson_P-pe$p.value)<1e-12,
 abs(rs$Spearman_P-sp$p.value)<1e-12,abs(rs$Pearson_r-cor(rna$CopyNumber,rna$Expression))>.1)
rp<-tcga_expression_plot(rna,"A",rs)[[1]]
stopifnot(rlang::as_label(rp$mapping$x)=="TCGA_Relative_CN_Change",
 identical(rp$layers[[2]]$data$xintercept,0))
# The fixed current mapping keeps COAD/READ separate, including empty cancers.
TCGA_CANCERS<-sort(unique(unlist(jsonlite::fromJSON("config/tcga_cancer_types.json")$mapping)))
all_order<-cancer_order(rna,TCGA_CANCERS);all_stats<-current_cancer_statistics(rna,all_order)
PATHS<-workflow_paths("geneA_screen",GENE_A)
rna_files<-wf_expected("tcga")[grepl("03_TCGA_CN_mRNA/",wf_expected("tcga"),fixed=TRUE)]
stopifnot(length(rna_files)==33L,all(c("COAD","READ") %in% TCGA_CANCERS),nrow(all_stats)==33L,
 all(all_stats$Status=="Insufficient N for correlation"))
pair<-as.data.table(prepare_cn(data.frame(ModelID=paste0("m",1:12),CN_relative=c(rep(.2,6),seq(1,2,length.out=6)))))
pair[,Chronos:=c(-2,-1.8,-1.7,-1.6,-1.5,-1.4,-.1,-.2,-.3,-.4,-.5,-.6)]
s<-targeted_statistics(pair)
scatter<-targeted_scatter_plot(pair,s);layers<-ggplot_build(scatter)$data
stopifnot(identical(layers[[2]]$yintercept,0),identical(layers[[3]]$xintercept,2^.585-1),
 identical(layers[[4]]$xintercept,2^.35-1),setequal(layers[[1]]$colour,c("black","red")))
pie<-tcga_cna_pie(data.table(CNA=factor(c(CNA_STATES[1],rep(CNA_STATES[3],9)),levels=CNA_STATES)))
stopifnot(inherits(pie$coordinates,"CoordPolar"),nrow(pie$data)==5L,sum(pie$data$N)==10L,
 identical(pie$scales$scales[[1]]$drop,FALSE),all(pie$data$pct==c(10,0,90,0,0)))
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
cp<-covariation_plot(ct)
stopifnot(all(vapply(c("Pearson r","r > 0","r < 0","|r|","correlation","causation"),
 function(term)grepl(term,cp$labels$caption,fixed=TRUE),logical(1))),
 identical(ggplot_build(cp)$data[[3]]$xintercept,c(0,0)))
# Exercise the complete SKIP wrapper at n=0,1,2 and either insufficient group.
# A guard proves the full Chronos matrix/statistical core is never run here.
RESULT_ROOT<-tempfile("wf_skip_",tmpdir=".runtime")
for(folder in c("Main_Results","Tables","Provenance"))dir.create(file.path(RESULT_ROOT,folder),recursive=TRUE,showWarnings=FALSE)
core_screen<-screen_core_statistics;screen_core_statistics<-function()stop("Screen must not run for insufficient groups")
fread<-function(input,nrows=NULL,select=NULL,...) {
 if(identical(nrows,0))return(matrix[0])
 stopifnot(identical(select,"ModelID"));matrix[,.(ModelID)]
}
for(nlow in c(0L,1L,2L,10L,11L,12L)) {
 fixture<-data.table(ModelID=pair$ModelID,CN_relative=c(rep(.2,nlow),rep(1,12L-nlow)))
 read_gene<-function(...)copy(fixture)
 skipped<-wf_run_screen()
 status<-jsonlite::fromJSON(wf_provenance("GenomeWide_Dependency_Status.json"))
 stopifnot(skipped$status=="SKIPPED",status$N_low==nlow,status$N_nonlow==12L-nlow,
  status$CN_log_threshold==.585,nrow(data.table::fread(wf_table("GenomeWide_Dependency.csv")))==0L,
  nrow(data.table::fread(wf_table("Top_Dependency_Candidates.csv")))==0L,
  all(vapply(file.path(RESULT_ROOT,wf_expected("genomewide_dependency")),wf_validate_artifact,logical(1))))
}
screen_core_statistics<-core_screen
unlink(RESULT_ROOT,recursive=TRUE)
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

# Data/output adapters. The supplied statistical function bodies remain intact.
OUTPUT_DIRS <- c(cn_expression="01_CN_Expression",genomewide_dependency="02_GenomeWide_Dependency",
 targeted_dependency="03_Targeted_Dependency",lineage_dependency="04_Lineage_Dependency",adjusted_dependency="05_Adjusted_Dependency",
 reverse_dependency="06_Reverse_Dependency",cn_covariation="07_CN_Covariation",mutation_dependency="09_Mutation_Dependency",qc="00_QC",tcga="08_TCGA")
MODULE_SKIP <- NULL
# Preserve numerical results while serializing small P/FDR values scientifically.
# Long fixed decimals from global scipen=999 can be misread as zero by CSV readers.
fwrite <- function(x,file,...)data.table::fwrite(x,file,...,scipen=0)
skip_module <- function(reason){MODULE_SKIP <<- reason;msg(reason);invisible(NULL)}
# Redirect folder literals only; no transformations of numeric/statistical expressions.
for(nm in c("run_cn_expression","run_targeted","run_genomewide","run_lineage","run_adjusted")){
 f<-get(nm);code<-paste(deparse(body(f),width.cutoff=500L),collapse="\n")
 for(key in names(OUTPUT_DIRS))code<-gsub(paste0('"',key,'"'),paste0('"',OUTPUT_DIRS[[key]],'"'),code,fixed=TRUE)
 code<-gsub('annotate("text",','annotate("label", fill="white", linewidth=0,',code,fixed=TRUE)
 body(f)<-parse(text=code)[[1]];assign(nm,f)
}
out_dir <- function(key) file.path(RESULT_ROOT,OUTPUT_DIRS[[key]])
move_output <- function(folder,old,new){
 from<-file.path(folder,old);to<-file.path(folder,new)
 if(file.exists(from)&&from!=to){if(file.exists(to))unlink(to);if(!file.rename(from,to))stop("Cannot rename result ",from)}
}
finalize_outputs <- function(mode){
 prefix<-paste0(GENE_A,"_",GENE_B)
 if(mode=="genomewide_dependency"){
  status<-system2("C:/Python312/python.exe",c(shQuote(file.path(PROJECT_ROOT,"scripts/utils/add_eligible_rank.py")),"--case",prefix,"--candidate",GENE_B))
  if(status!=0)stop("Eligible_Rank decoration failed")
 }
 if(mode=="depmap_cn_expression"){
  f<-out_dir("cn_expression");s<-fread(file.path(f,"CN_Expression_Statistics.csv"))
  cn<-read_gene(FILES$cn,GENE_A,"CN_relative");ex<-read_gene(FILES$expression,GENE_A,"Expression_A");dep<-read_gene(FILES$chronos,GENE_B,"Chronos_B")
  s[,N:=c(nrow(inner_join(cn,ex,by="ModelID") %>% filter(is.finite(CN_relative),is.finite(Expression_A))),
   nrow(inner_join(ex,dep,by="ModelID") %>% filter(is.finite(Expression_A),is.finite(Chronos_B))))]
  fwrite(s,file.path(f,"CN_Expression_Statistics.csv"))
  move_output(f,"CN_Expression_Dependency.pdf",paste0(GENE_A,"_CN_vs_Expression.pdf"))
 }
 if(mode=="targeted_dependency"){
  f<-out_dir("targeted_dependency")
  for(n in c("Statistics","Group_Summary","CellLines"))move_output(f,paste0(n,".csv"),paste0(prefix,"_",n,".csv"))
  move_output(f,"Targeted_Dependency.pdf",paste0(prefix,"_Targeted_Dependency.pdf"))
  run_three_groups()
 }
 if(mode=="lineage_dependency"){
  f<-out_dir("lineage_dependency")
  move_output(f,"Lineage_Dependency.csv",paste0(prefix,"_Lineage_Dependency.csv"))
  move_output(f,"Lineage_Dependency_Forest.pdf",paste0(prefix,"_Lineage_Forest.pdf"))
 }
}
run_reverse_safe <- function(){
 oldA<-GENE_A;oldB<-GENE_B;oldRoot<-RESULT_ROOT
 on.exit({GENE_A<<-oldA;GENE_B<<-oldB;RESULT_ROOT<<-oldRoot},add=TRUE)
 f<-out_dir("reverse_dependency")
 dir.create(file.path(f,OUTPUT_DIRS[["targeted_dependency"]]),recursive=TRUE,showWarnings=FALSE)
 GENE_A<<-oldB;GENE_B<<-oldA;RESULT_ROOT<<-f
 pair<-load_target_pair(GENE_A,GENE_B)
 counts<-table(pair$CN_binary)
 if(length(counts)<2L||any(counts<MIN_N))return(skip_module("Reverse dependency: at least one CN group has fewer than 3 models."))
 run_targeted()
 rev<-fread(file.path(f,OUTPUT_DIRS[["targeted_dependency"]],"Statistics.csv"))
 rev[,direction:="reverse"]
 forward_pair<-load_target_pair(oldA,oldB)
 corr<-cor.test(forward_pair$CN_log,forward_pair$Chronos)
 spear<-cor.test(forward_pair$CN_log,forward_pair$Chronos,method="spearman",exact=FALSE)
 wt<-wilcox.test(Chronos~CN_binary,data=forward_pair,exact=FALSE)
 fw<-data.table(geneA=oldA,geneB=oldB,N=nrow(forward_pair),pearson_r=unname(corr$estimate),pearson_p=corr$p.value,
  spearman_rho=unname(spear$estimate),spearman_p=spear$p.value,wilcoxon_p=wt$p.value,direction="forward")
 fw[,delta_median:=with(forward_pair,median(Chronos[CN_binary=="CN-Low"])-median(Chronos[CN_binary=="CN-NonLow"]))]
 rev[,delta_median:=with(pair,median(Chronos[CN_binary=="CN-Low"])-median(Chronos[CN_binary=="CN-NonLow"]))]
 fwrite(rbindlist(list(fw,rev),fill=TRUE),file.path(f,"Reverse_Dependency_Summary.csv"))
 for(n in c("Statistics.csv","Group_Summary.csv","CellLines.csv","Targeted_Dependency.pdf"))
  move_output(file.path(f,OUTPUT_DIRS[["targeted_dependency"]]),n,paste0(oldB,"_",oldA,"_",n))
}
run_three_groups <- function(){
 pair<-load_target_pair(GENE_A,GENE_B);counts<-table(factor(pair$CN_status,levels=c("CN Non-Low","Shallow CN Loss","Deep CN Loss")))
 status<-data.table(group=names(counts),N=as.integer(counts),eligible=as.integer(counts)>=MIN_N)
 fwrite(status,file.path(out_dir("targeted_dependency"),"Three_Group_Eligibility.csv"))
 if(any(counts<MIN_N))return(invisible(NULL))
 pair$CN_status<-factor(pair$CN_status,levels=names(counts))
 tests<-combn(names(counts),2,simplify=FALSE)
 stat<-rbindlist(lapply(seq_along(tests),function(i){g<-tests[[i]];p<-wilcox.test(pair$Chronos[pair$CN_status==g[1]],pair$Chronos[pair$CN_status==g[2]],exact=FALSE)$p.value
  data.table(group1=g[1],group2=g[2],P=p,y.position=max(pair$Chronos)+diff(range(pair$Chronos))*(.15*i))}))
 stat[,FDR:=p.adjust(P,"BH")];stat[,label:=paste0(vapply(P,p_star,character(1)),"\nP ",vapply(P,p_text,character(1)))]
 p<-ggplot(pair,aes(CN_status,Chronos,fill=CN_status))+geom_boxplot(outlier.shape=NA)+geom_jitter(width=.12,alpha=.25,size=1)+
  stat_pvalue_manual(stat,label="label",inherit.aes=FALSE)+theme_classic()+theme(legend.position="none")+
  labs(x=NULL,y=paste(GENE_B,"Chronos"),title=paste(GENE_A,"analysis-defined CN groups"))
 save_pdf(p,out_dir("targeted_dependency"),"Three_Group_Dependency",width=9,height=6)
 fwrite(stat,file.path(out_dir("targeted_dependency"),"Three_Group_Statistics.csv"))
}

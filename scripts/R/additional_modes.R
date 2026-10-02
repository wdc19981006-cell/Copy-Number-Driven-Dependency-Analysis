run_qc <- function(){
 model<-load_model();if(anyDuplicated(model$ModelID))stop("Duplicate ModelID in Model.csv")
 cn<-read_gene(FILES$cn,GENE_A,"CN_relative") %>% prepare_cn()
 ex<-read_gene(FILES$expression,GENE_A,"Expression_A");dep<-read_gene(FILES$chronos,GENE_B,"Chronos")
 dat<-cn %>% left_join(ex,by="ModelID") %>% left_join(dep,by="ModelID") %>% left_join(model,by="ModelID")
 f<-out_dir("qc");fwrite(dat,file.path(f,"QC_Model_Alignment.csv"))
 valid<-cn %>% filter(is.finite(CN_relative),is.finite(CN_log))
 counts<-data.table(CN_status=c("Deep CN Loss","Shallow CN Loss","CN Non-Low"),CN_binary=c("CN-Low","CN-Low","CN-NonLow"))
 counts[,N:=vapply(CN_status,function(g)sum(valid$CN_status==g),integer(1))]
 fwrite(counts,file.path(f,"CN_Group_Counts.csv"))
 p<-ggplot(cn %>% filter(is.finite(CN_relative)),aes(CN_relative))+geom_histogram(bins=50,fill="#3A6EA5",color="white")+
  theme_classic()+labs(title=paste(GENE_A,"relative CN distribution"),x="Relative CN",y="Models")
 q<-ggplot(cn %>% filter(is.finite(CN_log)),aes(CN_log))+geom_histogram(bins=50,fill="#C7473B",color="white")+
  geom_vline(xintercept=c(.35,.585),linetype="dashed")+theme_classic()+labs(x="log2(relative CN + 1)",y="Models")
 save_pdf(p+q,f,"CN_Distribution",width=11,height=5)
 measures<-data.table(metric=c("CN_models","expression_models","Chronos_models","CN_Chronos_overlap","CN_expression_overlap","all_three_overlap","CN_missing","CN_negative"),
  value=c(nrow(cn),nrow(ex),nrow(dep),length(intersect(cn$ModelID,dep$ModelID)),length(intersect(cn$ModelID,ex$ModelID)),
   length(Reduce(intersect,list(cn$ModelID,ex$ModelID,dep$ModelID))),sum(!is.finite(cn$CN_relative)),sum(cn$CN_relative<0,na.rm=TRUE)))
 fwrite(measures,file.path(f,"QC_Statistics.csv"))
 fwrite(data.table(variable=c("CN_relative","CN_log"),
  minimum=c(min(cn$CN_relative,na.rm=TRUE),min(cn$CN_log,na.rm=TRUE)),
  median=c(median(cn$CN_relative,na.rm=TRUE),median(cn$CN_log,na.rm=TRUE)),
  maximum=c(max(cn$CN_relative,na.rm=TRUE),max(cn$CN_log,na.rm=TRUE))),file.path(f,"CN_Distribution_Statistics.csv"))
}
run_cn_covariation <- function(){
 dt<-fread(FILES$cn,check.names=FALSE,showProgress=FALSE);id<-identify_id_column(names(dt))
 cols<-setdiff(names(dt),id);a<-find_gene_col(cols,GENE_A);x<-as.numeric(dt[[a]])
 res<-rbindlist(lapply(cols,function(g){y<-as.numeric(dt[[g]]);ok<-is.finite(x)&is.finite(y);n<-sum(ok)
  r<-if(n>=10)suppressWarnings(cor(x[ok],y[ok])) else NA_real_
  p<-if(is.finite(r)&&abs(r)<1)2*pt(-abs(r*sqrt((n-2)/(1-r^2))),df=n-2) else NA_real_
  data.table(Gene=clean_gene(g),N=n,Pearson_r=r,P=p)}))
 res[,FDR:=p.adjust(P,"BH")];setorder(res,-Pearson_r)
 fwrite(res,file.path(out_dir("cn_covariation"),paste0(GENE_A,"_CN_Covariation.csv")))
 top<-head(res[Gene!=GENE_A & is.finite(Pearson_r)],20)
 p<-ggplot(top,aes(Pearson_r,reorder(Gene,Pearson_r)))+geom_point(color="#3A6EA5")+theme_classic()+
  labs(x="Pearson r (relative CN)",y=NULL,title=paste(GENE_A,"CN covariation; association only"))
 save_pdf(p,out_dir("cn_covariation"),"Top_CN_Covariation",height=7)
 rm(dt);gc()
}
run_mutation <- function(){
 dep<-read_gene(FILES$chronos,GENE_B,"Chronos");stats<-list();plots<-list();skip<-list();cell<-list()
 for(kind in c("damaging","hotspot")){
  header<-names(fread(FILES[[kind]],nrows=0,check.names=FALSE));column<-find_gene_col(header,GENE_A,required=FALSE)
  if(is.null(column)){skip[[kind]]<-paste(kind,"gene absent from supplied export; mutation state unknown");next}
  mu<-read_gene(FILES[[kind]],GENE_A,"MutationValue")
  dat<-inner_join(mu,dep,by="ModelID") %>% filter(is.finite(MutationValue),is.finite(Chronos))
  if(any(dat$MutationValue<0))stop("Negative mutation indicators")
  dat$MutationStatus<-ifelse(dat$MutationValue>0,"Mutant","WT")
  low<-dat$Chronos[dat$MutationStatus=="Mutant"];non<-dat$Chronos[dat$MutationStatus=="WT"]
  if(length(low)<MIN_N||length(non)<MIN_N){skip[[kind]]<-paste(kind,"insufficient mutation groups: mutant",length(low),"WT",length(non));next}
  w<-wilcox.test(low,non,exact=FALSE)
  stats[[kind]]<-data.table(definition=kind,N_mutant=length(low),N_WT=length(non),Median_mutant=median(low),Median_WT=median(non),Delta_median=median(low)-median(non),P=w$p.value)
  dat$definition<-kind;cell[[kind]]<-as.data.table(dat)
  bracket<-data.frame(group1="WT",group2="Mutant",y.position=max(dat$Chronos)+.1*diff(range(dat$Chronos)),label=paste0(p_star(w$p.value),"\nP ",p_text(w$p.value)))
  plots[[kind]]<-ggplot(dat,aes(factor(MutationStatus,levels=c("WT","Mutant")),Chronos,fill=MutationStatus))+
   geom_boxplot(outlier.shape=NA)+geom_jitter(width=.12,alpha=.25,size=1)+stat_pvalue_manual(bracket,inherit.aes=FALSE)+
   theme_classic()+theme(legend.position="none")+labs(x=NULL,title=paste(GENE_A,kind,"mutation ->",GENE_B),y="Chronos")
 }
 f<-out_dir("mutation_dependency")
 if(length(skip))fwrite(data.table(definition=names(skip),reason=unlist(skip)),file.path(f,"Mutation_Skipped.csv"))
 if(!length(stats))return(skip_module(paste(unlist(skip),collapse="; ")))
 res<-rbindlist(stats);res[,FDR:=p.adjust(P,"BH")];fwrite(res,file.path(f,"Mutation_Dependency.csv"))
 fwrite(rbindlist(cell),file.path(f,"Mutation_CellLines.csv"))
 save_pdf(wrap_plots(plots),f,"Mutation_Dependency",width=6*length(plots),height=5)
}
write_case_summary <- function(){
 old_options<-options(scipen=0);on.exit(options(old_options),add=TRUE)
 prefix<-paste0(GENE_A,"_",GENE_B);summary_dir<-file.path(RESULT_ROOT,"Summary")
 metrics<-list();add<-function(module,stat,value){metrics[[length(metrics)+1L]]<<-data.table(module=module,statistic=stat,value=as.character(value))}
 lines<-c(paste0("# ",GENE_A," → ",GENE_B," analysis"),"",paste("R",getRversion(),"; DepMap 26Q1 user-supplied exports."),
  "","CN groups use analysis-defined thresholds on log2(relative CN + 1); these are not official GISTIC states.",
  "Continuous dependency and adjusted CN analyses use CN_log exactly as in the supplied R code; CN-expression and covariation use CN_relative.",
  "BH adjustments are separate for each complete screen/family; rank orders Wilcoxon FDR then delta median, exactly as supplied.",
  "Export completeness relative to official full-release files is unverified. All supplied Chronos columns are tested.","")
 read_if<-function(path)if(file.exists(path))as.data.table(readr::read_csv(path,show_col_types=FALSE,name_repair="minimal")) else NULL
 groups<-read_if(file.path(out_dir("qc"),"CN_Group_Counts.csv"))
 if(!is.null(groups))for(i in seq_len(nrow(groups))){add("QC",groups$CN_status[i],groups$N[i]);lines<-c(lines,paste(groups$CN_status[i],"N =",groups$N[i]))}
 if(!is.null(groups))lines<-c(lines,paste("Raw CN-Low N =",sum(groups$N[groups$CN_binary=="CN-Low"]),"; CN-NonLow N =",sum(groups$N[groups$CN_binary=="CN-NonLow"])))
 paired<-read_if(file.path(out_dir("targeted_dependency"),paste0(prefix,"_Group_Summary.csv")))
 if(!is.null(paired))for(i in seq_len(nrow(paired)))lines<-c(lines,paste("Matched",paired$CN_binary[i],"N =",paired$N[i]))
 candidates<-read_if(file.path(out_dir("genomewide_dependency"),paste0(GENE_B,"_Candidate_Rank.csv")))
 if(!is.null(candidates)&&nrow(candidates)){
  for(n in names(candidates))add("Genomewide",n,candidates[[n]][1])
  lines<-c(lines,"",paste(GENE_B,"rank =",candidates$Rank,"; delta median =",signif(candidates$Delta_median,5),"; Wilcoxon P =",signif(candidates$Wilcoxon_P,5),"; FDR =",signif(candidates$Wilcoxon_FDR,5)))
  complete<-read_if(file.path(out_dir("genomewide_dependency"),"GenomeWide_Dependency.csv"))
  if(!is.null(complete)){
   untested<-sum(is.na(complete$Wilcoxon_FDR));preceding<-sum(is.na(complete$Wilcoxon_FDR)&complete$Rank<candidates$Rank[1])
   add("Genomewide","undefined_Wilcoxon_FDR_rows",untested);add("Genomewide","undefined_rows_preceding_candidate",preceding)
   lines<-c(lines,paste("The supplied setorder() places NA FDR first:",preceding,"untested rows precede this candidate. Original Rank is retained unchanged. The candidate has the minimum eligible-test FDR:",isTRUE(candidates$Wilcoxon_FDR[1]==min(complete$Wilcoxon_FDR,na.rm=TRUE))))
  }
 }
 expr<-read_if(file.path(out_dir("cn_expression"),"CN_Expression_Statistics.csv"))
 if(!is.null(expr))for(i in seq_len(nrow(expr))){for(n in setdiff(names(expr),"analysis"))add(expr$analysis[i],n,expr[[n]][i]);lines<-c(lines,"",paste(expr$analysis[i],"Pearson",signif(expr$pearson[i],5),"P",signif(expr$pearson_p[i],5),"; Spearman",signif(expr$spearman[i],5),"P",signif(expr$spearman_p[i],5),"; N",expr$N[i]))}
 targeted<-read_if(file.path(out_dir("targeted_dependency"),paste0(prefix,"_Statistics.csv")))
 if(!is.null(targeted)){for(n in names(targeted))add("Targeted",n,targeted[[n]][1]);lines<-c(lines,"",paste("Targeted Pearson",signif(targeted$pearson_r,5),"P",signif(targeted$pearson_p,5),"; Spearman",signif(targeted$spearman_rho,5),"P",signif(targeted$spearman_p,5),"; Wilcoxon P",signif(targeted$wilcoxon_p,5),"; N",targeted$N))}
 for(n in c("Continuous_CN_Adjusted.csv","CNLow_Adjusted.csv")){
  adj<-read_if(file.path(out_dir("adjusted_dependency"),n));if(!is.null(adj)){
   a<-adj[grepl("^CN_log$|^I\\(",term)]
   for(i in seq_len(nrow(a))){for(k in setdiff(names(a),"term"))add(paste("Adjusted",n,a$term[i]),k,a[[k]][i]);lines<-c(lines,"",paste("Adjusted",a$term[i],"beta",signif(a$estimate[i],5),"SE",signif(a$std.error[i],5),"t",signif(a$statistic[i],5),"P",signif(a$p.value[i],5)))}
  }
 }
 line<-read_if(file.path(out_dir("lineage_dependency"),paste0(prefix,"_Lineage_Dependency.csv")))
 if(!is.null(line)){lines<-c(lines,"",paste("Eligible lineages:",nrow(line),"; FDR < 0.05:",sum(line$FDR<.05)),"Forest displays delta median with percentile 95% CI from 1,000 bootstrap draws (seed 1234).")
  for(i in seq_len(nrow(line))){for(n in setdiff(names(line),"Lineage"))add(paste("Lineage",line$Lineage[i]),n,line[[n]][i]);lines<-c(lines,paste(line$Lineage[i],": delta",signif(line$Delta_median[i],4),"CI [",signif(line$CI_low[i],4),",",signif(line$CI_high[i],4),"] FDR",signif(line$FDR[i],4)))}}
 rev<-read_if(file.path(out_dir("reverse_dependency"),"Reverse_Dependency_Summary.csv"))
 if(!is.null(rev))for(i in seq_len(nrow(rev))){for(n in names(rev))add(paste("Direction",rev$direction[i]),n,rev[[n]][i]);lines<-c(lines,"",paste(rev$direction[i],rev$geneA[i],"→",rev$geneB[i],"delta median",signif(rev$delta_median[i],5),"Pearson",signif(rev$pearson_r[i],5),"Wilcoxon P",signif(rev$wilcoxon_p[i],5)))}
 cov<-read_if(file.path(out_dir("cn_covariation"),paste0(GENE_A,"_CN_Covariation.csv")))
 if(!is.null(cov)){hit<-cov[Gene==GENE_B];for(n in names(hit))if(nrow(hit))add("CN covariation candidate",n,hit[[n]][1]);lines<-c(lines,"","CN covariation is an association/co-deletion signal; it does not establish physical proximity or causality.")}
 runs<-read_if(file.path(summary_dir,"Module_Runs.csv"));if(!is.null(runs)){lines<-c(lines,"","Module status (latest per mode):")
  latest<-runs[!duplicated(mode,fromLast=TRUE)];for(i in seq_len(nrow(latest)))lines<-c(lines,paste0("- ",latest$mode[i],": ",latest$status[i],if(is.na(latest$detail[i])) "." else paste0(". ",latest$detail[i])))}
 lines<-c(lines,"","These are observational cell-line associations. Lineage adjustment reduces measured lineage confounding; it does not establish a causal synthetic-lethal mechanism or clinical benefit.")
 if(length(metrics))fwrite(rbindlist(metrics),file.path(summary_dir,paste0(prefix,"_Key_Statistics.csv")))
 writeLines(lines,file.path(summary_dir,paste0(prefix,"_Summary.md")),useBytes=TRUE)
}

workflow_save <- function(p,path,width=9,height=6) {
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 ggplot2::ggsave(path,p,width=width,height=height,bg="white",device="pdf")
}
workflow_theme <- function()theme_classic(base_size=11)+theme(plot.title=element_text(face="bold"),
 plot.caption=element_text(hjust=0),legend.position="bottom")
tcga_prevalence_plot <- function(prev,order) {
 d<-copy(prev);d[,CancerType:=ordered_cancer(CancerType,order)]
 ggplot(d,aes(CancerType,Percentage/100,fill=CNA))+geom_col(width=.8)+coord_flip()+
 scale_fill_manual(values=CNA_COLORS,drop=FALSE)+scale_x_discrete(drop=FALSE)+
 # Clamp round-trip CSV precision at 100% instead of dropping the final slice.
 scale_y_continuous(labels=scales::percent,limits=c(0,1),expand=c(0,0),oob=scales::squish)+workflow_theme()+
 theme(plot.margin=margin(5.5,18,5.5,5.5))+
 labs(title=paste(GENE_A,"CNA Percentage Across TCGA Pan-Cancer"),
      x=NULL,y="Percentage of TCGA tumor samples",fill="Copy Number Status",
      caption="Cancer types follow ascending median relative CN change (CN / sample baseline - 1); source details in Provenance.")
}
tcga_landscape_plot <- function(dat,order) {
 d<-copy(dat[is.finite(CopyNumber)]);d[,CancerType:=ordered_cancer(CancerType,order)]
 labels<-setNames(paste0(order$CancerType," (n=",order$N,")"),order$CancerType)
 extent<-max(c(1,abs(d$TCGA_Relative_CN_Change)),na.rm=TRUE)
 ggplot(d,aes(CancerType,TCGA_Relative_CN_Change))+geom_boxplot(outlier.shape=NA,fill="white",width=.55)+
 geom_point(position=position_jitter(width=.17,height=0,seed=1234),color="grey45",alpha=.35,size=.55)+
 # coord_flip turns this into the vertical x=0 reference in the final figure.
 geom_hline(yintercept=0,linetype="dashed",color="grey25")+
 scale_y_continuous(limits=c(-extent,extent))+
 scale_x_discrete(labels=labels,drop=FALSE)+coord_flip()+workflow_theme()+
 labs(title=paste(GENE_A,"Relative Copy Number Change Across TCGA Pan-Cancer"),x=NULL,
      y=paste0(GENE_A," Relative Copy Number Change\n(CN / sample baseline - 1)"),
      caption=paste0("0 = sample-specific copy-number baseline;\nnegative values indicate relative copy-number loss; positive values indicate relative copy-number gain.\n",
                     "Analysis-derived change; each point represents one TCGA tumor sample."))
}
tcga_expression_plot <- function(dat,cancer,stat) {
 d<-copy(dat[CancerType==cancer & is.finite(CopyNumber)&is.finite(Expression)&is.finite(CNAState)])
 d[,CNA:=factor(CNAState,levels=-2:2,labels=CNA_STATES)]
 label<-if(stat$N<TCGA_MIN_N)paste0("N = ",stat$N,"\nInsufficient N for correlation") else
  if(!is.finite(stat$Pearson_r))paste0("N = ",stat$N,"\nConstant CN or mRNA") else
  paste0("N = ",stat$N,"\nPearson r = ",signif(stat$Pearson_r,3),"\nP = ",formatC(stat$Pearson_P,format="g",digits=3),
         "\nSpearman rho = ",signif(stat$Spearman_rho,3),"\nP = ",formatC(stat$Spearman_P,format="g",digits=3))
 p<-ggplot(d,aes(TCGA_Relative_CN_Change,Expression))+geom_point(aes(color=CNA),alpha=.55,size=1,show.legend=TRUE)+
 geom_vline(xintercept=0,linetype="dashed",color="grey25")+
 scale_color_manual(values=CNA_COLORS,drop=FALSE)+workflow_theme()+
 guides(color=guide_legend(override.aes=list(alpha=1,size=2)))+
 scale_y_continuous(expand=expansion(mult=c(.05,.65)))+
 labs(title=paste0(cancer,": ",GENE_A," Copy Number vs mRNA Expression"),
      x=paste0(GENE_A," Relative Copy Number Change\n(CN / sample baseline - 1)"),y=paste0(GENE_A," mRNA expression\nlog2(TPM + 1)"),
      color="Copy Number Status",caption="0 = sample-specific baseline; negative = relative loss; positive = relative gain.\nMatched TCGA tumor samples; correlations describe association.")
 if(nrow(d)>=2L&&length(unique(d$TCGA_Relative_CN_Change))>1L)p<-p+geom_smooth(method="lm",se=FALSE,color="black",linewidth=.65)
 # Draw the opaque statistics label after the regression layer, so no line
 # crosses the text. P values use scientific formatting independent of scipen.
 p<-p+annotate("label",x=Inf,y=Inf,label=label,hjust=1.02,vjust=1.08,size=3.1,linewidth=0,fill="white")
 inset<-tcga_cna_pie(d)
 p+patchwork::inset_element(inset,left=.015,bottom=.66,right=.31,top=.99,align_to="panel",on_top=TRUE)
}
tcga_cna_pie <- function(d) {
 counts<-as.data.table(table(factor(d$CNA,levels=CNA_STATES)));setnames(counts,c("CNA","N"))
 counts[,CNA:=factor(CNA,levels=CNA_STATES)]
 counts[,pct:=if(sum(N))100*N/sum(N) else 0]
 # Reserve an empty band above every observation for the pie and statistics.
 # Tiny/nonexistent slices have no text; the main legend always retains 5 states.
 inset<-ggplot(counts,aes(x="",y=N,fill=CNA))+geom_col(width=1,color="white",linewidth=.2)+
 coord_polar(theta="y")+geom_text(aes(label=ifelse(pct>=8,sprintf("%.1f%%",pct),"")),
 position=position_stack(vjust=.5),size=2.6)+scale_fill_manual(values=CNA_COLORS,drop=FALSE)+
 theme_void(base_size=9)+theme(legend.position="none",plot.title=element_text(hjust=.5),
 plot.margin=margin(0,0,0,0))+labs(title="CNA status")
 if(!nrow(d))inset<-ggplot()+annotate("text",0,0,label="No matched samples",size=3)+theme_void()+labs(title="CNA status")
 inset
}
targeted_scatter_plot <- function(pair,stat) {
 d<-copy(as.data.table(pair));d[,Group:=factor(CN_binary,levels=c("CN-NonLow","CN-Low"),
                        labels=paste(GENE_A,c("CN-Normal","CN-Low")))]
 colors<-setNames(c("black","red"),levels(d$Group))
 annotation<-paste0("Pearson r = ",signif(stat$Relative_CN_Pearson_r,3),
 "\nSpearman rho = ",signif(stat$Relative_CN_Spearman_rho,3),"\nN = ",stat$N)
 p<-ggplot(d,aes(CN_relative,Chronos))+geom_point(aes(color=Group),alpha=.55,size=1.2)+
 geom_hline(yintercept=0,linetype="dotted",color="grey55",linewidth=.5)+
 geom_vline(xintercept=2^.585-1,linetype="dashed",color="grey30",linewidth=.6)+
 geom_vline(xintercept=2^.35-1,linetype="dotted",color="black",linewidth=.6)+
 scale_color_manual(values=colors)+workflow_theme()+
 scale_y_continuous(expand=expansion(mult=c(.16,.08)))+
 annotate("label",x=2^.35-1,y=-Inf,label="Deep CN Loss",hjust=0,vjust=-1.3,size=3,linewidth=0,fill="white")+
 annotate("label",x=2^.585-1,y=-Inf,label="CN-Low threshold",hjust=0,vjust=-.05,size=3,linewidth=0,fill="white")+
 annotate("label",x=Inf,y=Inf,label=annotation,hjust=1.04,vjust=1.1,size=3.3,linewidth=0,fill="white")+
 labs(title=paste(GENE_B,"Dependency in Cells with Low",GENE_A,"Copy Number"),
 x=paste(GENE_A,"Relative Copy Number"),y=paste(GENE_B,"Chronos Gene Effect"),color=NULL,
 caption="More negative Chronos indicates stronger dependency. Deep CN Loss is a visual reference only.")
 p
}
targeted_plots <- function(pair,stat,paths) {
 d<-copy(as.data.table(pair));d[,Group:=factor(CN_binary,levels=c("CN-NonLow","CN-Low"),
                        labels=paste(GENE_A,c("CN-Normal","CN-Low")))]
 colors<-setNames(c("black","red"),levels(d$Group))
 workflow_save(targeted_scatter_plot(pair,stat),paths[1])
 w<-waterfall_order(d)
 p<-ggplot(w,aes(Waterfall_Order,Chronos,color=Group))+geom_segment(aes(xend=Waterfall_Order,yend=0),linewidth=.35)+
 geom_hline(yintercept=-1,linetype="dashed")+scale_color_manual(values=colors)+workflow_theme()+
 labs(title=paste(GENE_B,"Dependency Waterfall"),x="Matched cell lines: weaker to stronger dependency",y=paste(GENE_B,"Chronos Gene Effect"),
 color=NULL,caption="Sorted by descending Chronos; ties use ModelID. Dashed reference: Chronos = -1.")
 workflow_save(p,paths[2],width=10)
 n<-c(stat$N_nonlow,stat$N_low);med<-c(stat$Median_nonlow,stat$Median_low)
 labels<-setNames(paste0(levels(d$Group),"\nn=",n,"; median=",signif(med,3)),levels(d$Group))
 note<-paste0("Delta median (low - normal) = ",signif(stat$Delta_median,4),
              "; Wilcoxon P ",p_text(stat$Wilcoxon_P)," ",p_star(stat$Wilcoxon_P))
 if(!stat$eligible)note<-paste(note,"; Insufficient group N for test")
 p<-ggplot(d,aes(Group,Chronos))+geom_boxplot(aes(fill=Group),width=.55,outlier.shape=NA,alpha=.25)+
 geom_point(aes(color=Group),position=position_jitter(width=.12,height=0,seed=1234),alpha=.5,size=1)+
 scale_fill_manual(values=colors)+scale_color_manual(values=colors)+scale_x_discrete(labels=labels)+
 workflow_theme()+theme(legend.position="none")+labs(title=paste(GENE_B,"Dependency by",GENE_A,"Copy Number"),
 subtitle=note,x=NULL,y=paste(GENE_B,"Chronos Gene Effect"),caption="Each point is one matched cell line; more negative Chronos indicates stronger dependency.")
 workflow_save(p,paths[3]);invisible(w)
}
volcano_plot <- function(res) {
 d<-copy(res[is.finite(Delta_median)&is.finite(Wilcoxon_FDR)])
 d[,LogFDR:=-log10(pmax(Wilcoxon_FDR,1e-300))]
 d[,Direction:=ifelse(Wilcoxon_FDR<.05 & Delta_median<0,"Stronger dependency in CN-Low",
   ifelse(Wilcoxon_FDR<.05 & Delta_median>0,"Weaker dependency in CN-Low","Not significant"))]
 labels<-dependency_top(res,10)$Gene
 if(GENE_B_PROVIDED)labels<-unique(c(labels,GENE_B))
 d[,Label:=ifelse(Gene %in% labels,Gene,NA_character_)]
 p<-ggplot(d,aes(Delta_median,LogFDR))+geom_point(aes(color=Direction),alpha=.55,size=1.1)+
 geom_vline(xintercept=0,linetype="dashed")+geom_hline(yintercept=-log10(.05),linetype="dashed")+
 ggrepel::geom_text_repel(aes(label=Label),na.rm=TRUE,max.overlaps=Inf,size=3,seed=1234)+
 scale_color_manual(values=c("Stronger dependency in CN-Low"="#C7473B","Weaker dependency in CN-Low"="#3A6EA5","Not significant"="grey75"))+
 workflow_theme()+labs(title=paste(GENE_A,"CN-Low Genome-wide Dependency Screen"),
 x="Delta median Chronos (CN-low - CN-nonlow)",y="-log10(Wilcoxon BH FDR)",color=NULL,
 caption="Negative delta: stronger dependency in CN-low. Top 10 per direction labeled; Eligible_Rank in Tables.")
 if(GENE_B_PROVIDED)p<-p+geom_point(data=d[Gene==GENE_B],shape=21,fill="gold",color="black",size=3)
 p
}
covariation_plot <- function(top) {
 d<-copy(top);d[,Gene:=factor(Gene,levels=unique(Gene[order(Pearson_r)]))]
 d[,Direction:=factor(Direction,levels=c("Positive CN correlation","Negative CN correlation"))]
 ggplot(d,aes(Pearson_r,Gene,color=Direction))+geom_segment(aes(x=0,xend=Pearson_r,yend=Gene),linewidth=.6)+
 geom_point(size=2)+geom_vline(xintercept=0,linetype="dashed")+
 facet_wrap(~Direction,scales="free",nrow=1,labeller=as_labeller(c(
  "Positive CN correlation"=paste0("Positive CN correlation\n",GENE_A," low: gene usually low"),
  "Negative CN correlation"=paste0("Negative CN correlation\n",GENE_A," low: gene usually high"))))+
 scale_color_manual(values=c("Positive CN correlation"="#3A6EA5","Negative CN correlation"="#C7473B"))+
 workflow_theme()+theme(legend.position="none")+labs(title=paste0("Genes Positively and Negatively Covarying with\n",GENE_A," Copy Number"),
 x="Pearson r (relative copy number)",y=NULL,
 caption=paste0("Pearson r = linear correlation coefficient between ",GENE_A," copy number and each gene's copy number\n",
  "across DepMap cell lines. Top 20 in each direction.\n",
  "r > 0: same-direction copy-number change; r < 0: opposite-direction change;\n",
  "|r| closer to 1 indicates stronger correlation; r near 0 indicates weaker linear correlation.\n",
  "Correlation indicates covariation, not causation."))
}

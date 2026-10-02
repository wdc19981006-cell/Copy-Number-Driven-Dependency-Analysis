# Python bridge reads only requested genes from Parquet; no R arrow dependency.
tcga_ready <- function(){
 cfg<-jsonlite::fromJSON(file.path(PROJECT_ROOT,"config/tcga_layers.json"))
 live<-file.path(PROJECT_ROOT,"data/raw/tcga/gdc_current_DR46/manifests/live/current_layer_status.json")
 if(file.exists(live))cfg<-jsonlite::fromJSON(live)
 if(!identical(cfg$current_status,"complete"))return(FALSE)
 TRUE
}
tcga_gene <- function(dataset,layer="current"){
 target<-file.path(out_dir("tcga"),paste0(GENE_A,"_",dataset,"_",layer,".csv"))
 status<-system2("C:/Python312/python.exe",c(shQuote(file.path(PROJECT_ROOT,"scripts/utils/export_tcga_gene.py")),
  "--gene",GENE_A,"--dataset",dataset,"--layer",layer,"--output",shQuote(target)))
 if(status!=0)stop("TCGA gene export failed")
 fread(target)
}
run_tcga_landscape <- function(){
 if(!tcga_ready())return(skip_module("TCGA module skipped because current GDC data are not ready."))
 dat<-tcga_gene("cn");fwrite(dat,file.path(out_dir("tcga"),"TCGA_CN_Landscape.csv"))
 if(anyDuplicated(dat$SampleID))stop("Duplicate selected CN sample UUID")
 stat<-dat[,.(N_samples=.N,N_finite=sum(is.finite(Value)),CN_median=median(Value,na.rm=TRUE)),by=.(CancerType,Workflow,GDCRelease,DataLayer)]
 fwrite(stat,file.path(out_dir("tcga"),"TCGA_CN_Landscape_Statistics.csv"))
 p<-ggplot(dat,aes(reorder(CancerType,Value,median,na.rm=TRUE),Value))+geom_boxplot(outlier.shape=NA)+coord_flip()+theme_classic()+labs(x=NULL,y="GDC gene-level copy number",title=paste(GENE_A,"GDC DR46 CN by cancer"))
 save_pdf(p,out_dir("tcga"),"TCGA_CN_Landscape",height=9)
}
run_tcga_prevalence <- function(){
 dat<-tcga_gene("gistic",layer="reference");dat<-dat[TumorNormal=="Tumor" & !is.na(CancerType) & is.finite(Value)]
 if(!nrow(dat))stop("No tumor observations in reference GISTIC; check metadata classification")
 dat[,CNA:=factor(Value,levels=-2:2,labels=c("Deep deletion","Shallow deletion","Diploid","Gain","Amplification"))]
 res<-dat[,.(N=.N),by=.(CancerType,CNA)];res[,denominator:=sum(N),by=CancerType];res[,prevalence:=N/denominator]
 fwrite(res,file.path(out_dir("tcga"),"TCGA_GISTIC_CNA_Prevalence.csv"))
 p<-ggplot(res,aes(CancerType,prevalence,fill=CNA))+geom_col()+coord_flip()+theme_classic()+scale_y_continuous(labels=scales::percent)+labs(x=NULL,y="Prevalence",title=paste(GENE_A,"PanCanAtlas reference GISTIC (not current GDC)"))
 save_pdf(p,out_dir("tcga"),"TCGA_GISTIC_CNA_Prevalence",height=9)
}
run_tcga_expression <- function(){
 if(!tcga_ready())return(skip_module("TCGA module skipped because current GDC data are not ready."))
 cn<-tcga_gene("cn");ex<-tcga_gene("expression")
 # sample UUIDs are the join key; aliquots remain in the exported provenance table.
 if(anyDuplicated(cn$SampleID))stop("CN representative sample UUID is not unique")
 # Retain all raw RNA aliquots; select one reproducible representative for this
 # sample-level analysis. Prefer the selected CN aliquot, then lexical file UUID.
 ex<-as.data.table(ex);ex[,matches_CN_aliquot:=AliquotID==cn$AliquotID[match(SampleID,cn$SampleID)]]
 ex[is.na(matches_CN_aliquot),matches_CN_aliquot:=FALSE]
 setorder(ex,SampleID,-matches_CN_aliquot,FileID)
 ex[,selected_for_sample_analysis:=!duplicated(SampleID)]
 fwrite(ex[,.(SampleID,AliquotID,FileID,matches_CN_aliquot,selected_for_sample_analysis,GDCRelease="DR46",DataLayer="gdc_DR46")],file.path(out_dir("tcga"),"TCGA_RNA_Representative_Selection.csv"))
 ex<-ex[selected_for_sample_analysis==TRUE]
 dat<-inner_join(cn,ex,by="SampleID",suffix=c("_CN","_RNA")) %>% filter(is.finite(Value_CN),is.finite(Value_RNA))
 if(anyDuplicated(dat$SampleID))stop("Sample representative selection failed")
 fwrite(dat,file.path(out_dir("tcga"),"TCGA_CN_Expression_Samples.csv"))
 pe<-cor.test(dat$Value_CN,dat$Value_RNA);sp<-cor.test(dat$Value_CN,dat$Value_RNA,method="spearman",exact=FALSE)
 fwrite(data.table(N=nrow(dat),Pearson_r=unname(pe$estimate),Pearson_P=pe$p.value,Spearman_rho=unname(sp$estimate),Spearman_P=sp$p.value,GDCRelease="DR46",DataLayer="gdc_DR46"),file.path(out_dir("tcga"),"TCGA_CN_Expression_Statistics.csv"))
 p<-ggplot(dat,aes(Value_CN,Value_RNA))+geom_point(alpha=.25,size=.8)+geom_smooth(method="lm",se=FALSE)+theme_classic()+labs(x="GDC gene-level CN",y="log2(TPM + 1)",title=paste(GENE_A,"GDC DR46 CN-expression"))
 save_pdf(p,out_dir("tcga"),"TCGA_CN_Expression")
}

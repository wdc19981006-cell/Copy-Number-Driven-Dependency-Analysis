tests<-c("test_r_core.R","test_cn_threshold_sensitivity.R","test_expression_dependency.R","test_tcga_by_cancer.R","test_genomewide_adjusted.R","test_tcga_sample_selection.R")
reports<-list()
for(name in tests){
 cat("\nRunning",name,"\n")
 status<-system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",file.path("tests",name)))
 reports[[name]]<-data.frame(test=name,status=if(status==0)"PASS" else "FAIL",exit_code=status)
 if(status!=0)stop(name," failed")
}
dir.create("data/manifests",recursive=TRUE,showWarnings=FALSE)
write.csv(do.call(rbind,reports),"data/manifests/synthetic_test_results.csv",row.names=FALSE)

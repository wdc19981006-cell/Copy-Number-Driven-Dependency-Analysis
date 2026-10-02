invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"));source(".Rprofile")
options(repos=c(CRAN="https://cloud.r-project.org"))
results<-list()
for(p in c("arrow","duckdb")){
 local<-list.files(".runtime/R-optional",pattern=paste0("^",p,"_.*\\.zip$"),full.names=TRUE)
 ok<-tryCatch({if(length(local))utils::install.packages(normalizePath(local[1],winslash="/"),repos=NULL,type="win.binary",lib=.libPaths()[1]) else renv::install(p,type="binary",prompt=FALSE);requireNamespace(p,quietly=TRUE)},error=function(e){message(p,": ",conditionMessage(e));FALSE})
 results[[p]]<-ok
}
print(results)
if(all(unlist(results))){
 renv::snapshot(packages=c("data.table","dplyr","tidyr","ggplot2","ggpubr","ggrepel","patchwork","broom","scales","stringr","purrr","readr","jsonlite","digest","optparse","future.apply","renv","arrow","duckdb"),prompt=FALSE)
}

Sys.setlocale("LC_CTYPE","English_United States.utf8")
Sys.setenv(LANG="English_United States.utf8",LC_ALL="English_United States.utf8",LC_CTYPE="English_United States.utf8")
source(".Rprofile")
Sys.setenv(RTOOLS45_HOME="D:/R/rtools45",PATH=paste("D:/R/rtools45/usr/bin","D:/R/rtools45/x86_64-w64-mingw32.static.posix/bin",Sys.getenv("PATH"),sep=";"))
dir.create(".runtime/compatible-src",recursive=TRUE,showWarnings=FALSE)
versions<-c(rlang="1.1.7",cli="3.6.4")
for(p in names(versions)){
 ver<-versions[[p]];dest<-file.path(".runtime/compatible-src",paste0(p,"_",ver,".tar.gz"))
 url<-paste0("https://cloud.r-project.org/src/contrib/Archive/",p,"/",p,"_",ver,".tar.gz")
 if(!file.exists(dest))download.file(url,dest,mode="wb")
 utils::install.packages(normalizePath(dest,winslash="/"),repos=NULL,type="source",lib=.libPaths()[1],INSTALL_opts="--no-test-load")
}

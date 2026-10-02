Sys.setlocale("LC_CTYPE","English_United States.utf8")
Sys.setenv(LANG="English_United States.utf8",LC_ALL="English_United States.utf8",LC_CTYPE="English_United States.utf8")
source(".Rprofile")
Sys.setenv(RTOOLS45_HOME="D:/R/rtools45",PATH=paste("D:/R/rtools45/usr/bin","D:/R/rtools45/x86_64-w64-mingw32.static.posix/bin",Sys.getenv("PATH"),sep=";"))
options(repos=c(CRAN="https://cloud.r-project.org"))
print(Sys.which(c("make","gcc")))
for(p in c("rlang","cli")){
 ver<-as.character(packageVersion(p))
 tar<-file.path(getwd(),".runtime/renv-root/source/repository",p,paste0(p,"_",ver,".tar.gz"))
 stopifnot(file.exists(tar))
 utils::install.packages(tar,repos=NULL,type="source",lib=.libPaths()[1])
}

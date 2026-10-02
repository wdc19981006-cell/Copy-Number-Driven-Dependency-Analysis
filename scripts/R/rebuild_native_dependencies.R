Sys.setlocale("LC_CTYPE","English_United States.utf8")
source(".Rprofile")
stopifnot(as.character(getRversion())=="4.5.0")
Sys.setenv(RTOOLS45_HOME="D:/R/rtools45")
Sys.setenv(PATH=paste("D:/R/rtools45/usr/bin","D:/R/rtools45/x86_64-w64-mingw32.static.posix/bin",Sys.getenv("PATH"),sep=";"))
options(repos=c(CRAN="https://cloud.r-project.org"))
# Rebuild the same locked versions with R 4.5.0; do not change statistical code.
for(p in c("rlang","cli")){
 version<-as.character(packageVersion(p))
 renv::install(paste0(p,"@",version),type="source",rebuild=TRUE,prompt=FALSE)
}

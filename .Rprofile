invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
if(.Platform$OS.type=="windows" && !nzchar(Sys.getenv("PROCESSOR_ARCHITECTURE"))){
 arch<-switch(R.version$arch,x86_64="AMD64",aarch64="ARM64",stop("Unknown Windows runtime architecture"))
 Sys.setenv(PROCESSOR_ARCHITECTURE=arch)
}
Sys.setenv(RENV_PATHS_CACHE=file.path(getwd(),".runtime/renv-cache"),RENV_PATHS_ROOT=file.path(getwd(),".runtime/renv-root"),
 RENV_PATHS_LIBRARY=file.path(getwd(),"renv/library"),RENV_PATHS_SANDBOX=file.path(getwd(),".runtime/renv-sandbox"))
source("renv/activate.R")

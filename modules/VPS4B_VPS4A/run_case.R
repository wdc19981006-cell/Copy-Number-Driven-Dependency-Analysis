#!/usr/bin/env Rscript
invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
stopifnot(as.character(getRversion())=="4.5.0")
args<-commandArgs(trailingOnly=FALSE)
script<-sub("^--file=","",args[grepl("^--file=",args)][1])
root<-normalizePath(file.path(dirname(script),"../.."),winslash="/",mustWork=TRUE)
extra<-commandArgs(trailingOnly=TRUE)
mode<-if(length(extra))extra[1] else "full"
executable<-file.path(R.home("bin"),"Rscript.exe")
if(!file.exists(executable))executable<-file.path(R.home("bin"),"Rscript")
status<-system2(executable,c("--vanilla",shQuote(file.path(root,"scripts/R/run_analysis.R")),
 "--project",shQuote(root),"--mode",mode,"--geneA","VPS4B","--geneB","VPS4A"))
quit(status=status)

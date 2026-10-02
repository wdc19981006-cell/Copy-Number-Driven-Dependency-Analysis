Sys.setlocale("LC_CTYPE","English_United States.utf8");source(".Rprofile");args<-commandArgs(TRUE);for(p in args)library(p,character.only=TRUE);cat("PROBE FINISHED\n")

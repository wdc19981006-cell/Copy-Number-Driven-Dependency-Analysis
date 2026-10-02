invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"));source(".Rprofile");library(data.table)
for(root in c("results/ENO1_ENO2","results/VPS4B_VPS4A")){
 for(file in list.files(root,pattern="\\.csv$",recursive=TRUE,full.names=TRUE)){
  dat<-as.data.table(readr::read_csv(file,show_col_types=FALSE,name_repair="minimal"))
  if(any(grepl("(^P$|_P$|_FDR$|_p$|^p.value$|^FDR$)",names(dat))))fwrite(dat,file,scipen=0)
 }
}
cat("Statistical CSV values preserved; P/FDR notation serialized with scipen=0.\n")

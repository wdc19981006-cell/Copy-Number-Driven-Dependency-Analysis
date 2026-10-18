# Shared Gene A project bookkeeping. No data access or statistical execution.
wf_geneBs <- function()sort(basename(list.dirs(file.path(RESULT_ROOT,"Main_Results/06_GeneB_Dependency"),recursive=FALSE)))
wf_project_expected <- function() {
 saved<-list(PATHS=PATHS,GENE_B=GENE_B,GENE_B_PROVIDED=GENE_B_PROVIDED)
 on.exit(list2env(saved,envir=.GlobalEnv))
 PATHS<<-workflow_paths("geneA_screen",GENE_A)
 expected<-unname(unlist(lapply(wf_base_modules(),wf_expected)))
 for(b in wf_geneBs()) {
  GENE_B<<-b;GENE_B_PROVIDED<<-TRUE;PATHS<<-workflow_paths("geneA_geneB",GENE_A,b)
  expected<-c(expected,unname(unlist(lapply(wf_geneB_modules(),wf_expected))),
   file.path("Main_Results/06_GeneB_Dependency",b,"00_GeneB_Summary.txt"))
 }
 expected
}
wf_geneB_summary <- function() {
 s<-data.table::fread(wf_table(paste0(GENE_A,"_",GENE_B,"_Targeted_Statistics.csv")))
 paths<-unname(unlist(lapply(wf_geneB_modules(),wf_expected)))
 lines<-c(paste0(GENE_A," → ",GENE_B," Dependency"),
  paste0("分析目的：检验 ",GENE_A," relative CN 与 ",GENE_B," Chronos dependency 的关联。"),
  paste0("N=",s$N,"；CN-Low N=",s$N_low,"；CN-Normal N=",s$N_nonlow),
  paste0("Relative CN Pearson r=",format(s$Relative_CN_Pearson_r,digits=16),
         "；Spearman rho=",format(s$Relative_CN_Spearman_rho,digits=16)),
  paste0("Median Chronos：low=",format(s$Median_low,digits=16),"；normal=",format(s$Median_nonlow,digits=16)),
  paste0("Delta median=",format(s$Delta_median,digits=16),"；Wilcoxon P=",format(s$Wilcoxon_P,digits=16)),
  "分组仍使用原始 CN_log < 0.585；统计仍使用原始观测值；相关不代表因果。",
  "三张主图：",unname(PATHS[c("scatter","waterfall","box")]),
  "统计CSV位置：",paths[grepl("\\.csv$",paths)],
  paste0("Supplementary位置：Supplementary/",GENE_B,"/"),
  "补充分析：Lineage、Adjusted、CN Threshold Sensitivity。")
 path<-file.path(RESULT_ROOT,"Main_Results/06_GeneB_Dependency",GENE_B,"00_GeneB_Summary.txt")
 writeLines(enc2utf8(lines),path,useBytes=TRUE)
 invisible(lines)
}
wf_summary <- function() {
 # Retain the original Gene A narrative byte-for-byte between appends.
 base_path<-wf_provenance("GeneA_Summary.txt")
 if(!file.exists(base_path)) {
  saved<-list(PATHS=PATHS,GENE_B=GENE_B,GENE_B_PROVIDED=GENE_B_PROVIDED)
  on.exit(list2env(saved,envir=.GlobalEnv),add=TRUE)
  PATHS<<-workflow_paths("geneA_screen",GENE_A);GENE_B<<-"";GENE_B_PROVIDED<<-FALSE
  path<-file.path(RESULT_ROOT,"00_Analysis_Summary.txt")
  if(!file.exists(path))wf_summary_base()
  lines<-readLines(path,encoding="UTF-8",warn=FALSE)
  i<-match("File Index",lines)
  if(is.na(i))stop("Cannot preserve existing Gene A summary: missing File Index")
  writeLines(enc2utf8(lines[seq_len(i-3L)]),base_path,useBytes=TRUE)
  list2env(saved,envir=.GlobalEnv)
 }
 if(GENE_B_PROVIDED)wf_geneB_summary()
 lines<-readLines(base_path,encoding="UTF-8",warn=FALSE)
 for(b in wf_geneBs()) {
  path<-file.path(RESULT_ROOT,"Main_Results/06_GeneB_Dependency",b,"00_GeneB_Summary.txt")
  lines<-c(lines,"",paste0("06. ",GENE_A," → ",b," Dependency"),readLines(path,encoding="UTF-8",warn=FALSE))
 }
 files<-sort(unique(c(list.files(RESULT_ROOT,recursive=TRUE),"00_Analysis_Summary.txt","Provenance/File_Index.csv")))
 lines<-c(lines,"","File Index",files[!startsWith(files,"Provenance/")],
          "技术审计、版本和验证信息保存在Provenance",files[startsWith(files,"Provenance/")])
 writeLines(enc2utf8(lines),file.path(RESULT_ROOT,"00_Analysis_Summary.txt"),useBytes=TRUE)
 data.table::fwrite(data.table(Path=files),wf_provenance("File_Index.csv"))
 invisible(TRUE)
}

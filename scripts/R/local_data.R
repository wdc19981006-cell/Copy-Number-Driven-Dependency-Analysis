# Runtime I/O adapter: preserve supplied statistical function bodies while
# directing every ordinary mode to existing processed Parquet only.
DATA_MODE <- "local"
LOCAL_UNAVAILABLE <- paste("Local processed dataset unavailable.",
 "Run the separate data-update/preparation workflow first.",sep="\n")
local_required <- function(path) {
 if(!file.exists(path))stop(LOCAL_UNAVAILABLE,"\nMissing: ",path,call.=FALSE)
 invisible(path)
}
depmap_files <- function(root) {
 datasets<-c(model="Model",condition="ModelCondition",profiles="OmicsProfiles",cn="OmicsCNGeneWGS",
 expression="ExpressionProteinCoding",chronos="CRISPRGeneEffect",dependency="CRISPRGeneDependency",
 damaging="MutationDamaging",hotspot="MutationHotspot",signatures="GlobalSignatures",subtype="MolecularSubtypes")
 as.list(setNames(file.path(root,"data/processed/depmap/26Q1",paste0(datasets,".parquet")),names(datasets)))
}
fread <- function(input,...,nrows=NULL,select=NULL,check.names=FALSE,showProgress=FALSE) {
 if(is.character(input)&&length(input)==1L&&grepl("\\.parquet$",input)) {
  local_required(input)
  nms<-arrow::ParquetFileReader$create(input)$GetSchema()$names
  if(identical(nrows,0)||identical(nrows,0L))return(setNames(as.data.table(rep(list(character()),length(nms))),nms))
  if(!is.null(select)) {
   if(is.numeric(select))select<-nms[select]
   x<-arrow::read_parquet(input,col_select=tidyselect::all_of(select))
  }else x<-arrow::read_parquet(input)
  x<-as.data.table(x)
  if(!is.null(nrows))x<-head(x,nrows)
  return(x)
 }
 if(is.character(input)&&any(grepl("(^|[/\\\\])data[/\\\\]raw[/\\\\]",input)))
  stop("DATA_MODE=local forbids raw data reads in ordinary analysis",call.=FALSE)
 data.table::fread(input,...,nrows=if(is.null(nrows))Inf else nrows,select=select,
                   check.names=check.names,showProgress=showProgress)
}
tcga_ready <- function() {
 local_required(file.path(PROJECT_ROOT,"data/processed/tcga/gdc_DR46/TCGA_GeneLevel_CN.parquet"))
 local_required(file.path(PROJECT_ROOT,"data/processed/tcga/gdc_DR46/TCGA_STAR_log2TPMplus1.parquet"))
 TRUE
}

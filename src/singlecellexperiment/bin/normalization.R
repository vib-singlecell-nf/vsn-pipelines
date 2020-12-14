suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(scran))

option_list = list(
  make_option(
    "--sceObj",
    default=NA,
    type='character',
    help="File path to the rds file containing a SCE object."
  ),
  make_option(
    "--useSumFactors",
    default=FALSE,
    action = "store_true",
    help="compute sum factors before normalization"
  ),
  make_option(
    "--quickCluster",
    default=FALSE,
    action = "store_true",
    help="use the quick.cluster function before computing sum factors"
  ),
  make_option(
    "--output",
    default=NA,
    type='character',
    help="output file name"
  )
)
opt <- parse_args(OptionParser(option_list=option_list))

if(!is.na(opt$sceObj)){
  sce <- readRDS(opt$sceObj)
  if(opt$useSumFactors){
	  if(opt$quickCluster){
	    set.seed(123)
	    q.clust<- quickCluster(sce)
	    sce <- scran::computeSumFactors(sce, cluster=q.clust)
	  }else{
	    sce <- scran::computeSumFactors(sce)
	  }
  }

  sce <- scater::logNormCounts(sce)

  sce@metadata$tools$diagnostics[['dimBeforeSeuratObj']]<-paste0(nrow(sce)," genes - ",ncol(sce)," cells")

  if(is.na(opt$output)){
    filename = paste0("SINGLE_CELL_EXPERIMENT__NORMALIZATION.rds")
    saveRDS(sce,file=filename, compress = T)
  } else {
    saveRDS(sce,file=as.character(opt$output), compress = T)
  }
} else {
  stop("Please provide a path to a .rds file containing a SCE object")
}

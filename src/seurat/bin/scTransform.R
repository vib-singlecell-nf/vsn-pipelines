suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))

option_list = list(
  make_option(
    "--seuratObj",
    default=NA,
    type='character',
    help="File path to the rds file containing a Seurat object."
  ),
  make_option(
    "--output",
    default=NA,
    type='character',
    help="output file name"
  ),
  make_option(
    "--regressSubsamples",
    default=FALSE,
    action = "store_true",
    help="Add subsamples as variable to regress in SCTransform"
  )
)
opt <- parse_args(OptionParser(option_list=option_list))

seuratObj <- readRDS(file = opt$seuratObj)

################################
########## SCT transformation
################################

# does Normalizedata, scaledata and findvariablefeatures in one

if(opt$regressSubsamples){
  # subsamples could be added as batch to regress out
  seuratObj <- SCTransform(seuratObj, verbose = TRUE,  vars.to.regress = "subsamples", new.assay.name = "SCT", return.only.var.genes = F)
}else {
  seuratObj <- SCTransform(seuratObj, verbose = TRUE, new.assay.name = "SCT", return.only.var.genes = F)
}

# Outputs writing into rds file
if(is.na(opt$output)){
  filename = paste0(seuratObj@project.name,".SEURAT__SCTRANSFORM.rds")
  saveRDS(seuratObj,file=filename, compress = T)
} else {
  saveRDS(seuratObj,file=as.character(opt$output), compress = T)
}

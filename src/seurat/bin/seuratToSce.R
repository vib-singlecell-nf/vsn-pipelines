suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SingleCellExperiment))

option_list = list(
  make_option(
    "--seuratObj",
    default=NA,
    type='character',
    help="File path to the rds file containing the seurat object."
  ),
  make_option(
    "--output",
    default=NA,
    type='character',
    help="output file name"
  )
)

opt <- parse_args(OptionParser(option_list=option_list)) 

seuratObj <- readRDS(file = opt$seuratObj)

sceObj <- as.SingleCellExperiment(seuratObj)

sceObj@metadata$project.name <- seuratObj@project.name
sceObj@metadata$misc <- seuratObj@misc
sceObj@metadata$tools <- seuratObj@tools

if(is.na(opt$output)){
  filename = paste0(seuratObj@project.name,"_SEURAT__SEURAT_TO_SCE.rds")
  saveRDS(sceObj,file=filename, compress = T)
} else {
  saveRDS(sceObj,file=as.character(opt$output), compress = T)
}
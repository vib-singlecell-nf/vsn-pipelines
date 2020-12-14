suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(gridExtra))

option_list = list(
  make_option(
    "--sceObj",
    default=NA,
    type='character',
    help="File path to the rds file containing a SCE object."
  ),
  make_option(
    "--seuratObj",
    default=NA,
    type='character',
    help="To merge the SCE object, file path to a rds file containing a Seurat object."
  ),
  make_option(
    "--output",
    default=NA,
    type='character',
    help="output file name"
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

sce <- readRDS(opt$sceObj)
seuratObj_tmp <- as.Seurat(sce, assay = "RNA", counts = "counts", data = "logcounts")

if(!is.na(opt$seuratObj)){
  seuratObj <- readRDS(file = opt$seuratObj)
  cells.to.keep <- colnames(sce)
  seuratObj<- subset(seuratObj, cells=cells.to.keep)
  seuratObj[["RNA"]] <- seuratObj_tmp[["RNA"]]
  diffVariable <- !colnames(seuratObj_tmp@meta.data) %in% colnames(seuratObj@meta.data)
  meta_temp <- seuratObj_tmp@meta.data[,diffVariable]
  metaD <- merge(seuratObj@meta.data,meta_temp,by=0)
  rownames(metaD) <- metaD$Row.names
  metaD$Row.names <- NULL
  seuratObj@meta.data <- metaD
}else {
  seuratObj <- seuratObj_tmp
  seuratObj@project.name <- sce@metadata$project.name
}

seuratObj@tools$diagnostics <- c(seuratObj@tools$diagnostics,sce@metadata$diagnostics)
seuratObj@misc <- sce@metadata$misc

##### Add to diagnostics #####
seuratObj@tools$diagnostics[['dimAfterSeuratObj']]<-paste0(nrow(seuratObj@assays$RNA@counts),
                                                                   " genes - ",
                                                                   ncol(seuratObj@assays$RNA@counts),
                                                                   " cells")
# Outputs writing into rds file
if(is.na(opt$output)){
  filename = paste0(seuratObj@project.name,"_SEURAT__SCE_TO_SEURAT.rds")
  saveRDS(seuratObj,file=filename, compress = T)
} else {
  saveRDS(seuratObj,file=as.character(opt$output), compress = T)
}
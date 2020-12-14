#!/usr/bin/env Rscript

# Libraries loading
suppressPackageStartupMessages(library(reticulate))
use_python("/usr/bin/python3")
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(optparse))

# opt/arg parser definition
option_list = list(
  make_option(c("-i", "--inputSeuratRds"),
              default=NA,
              type='character',
              help="Path to a RDS file containing a Seurat object"),
  make_option(c("-o", "--output"),
              default=NA,
              type='character',
              help="Output name"),
  make_option(c("-d", "--ndims"),
              default=10,
              type='numeric',
              help="number of dimensions used to calculte the UMAP")

)

opt = parse_args(OptionParser(option_list=option_list))

# Input loading
sobj <- readRDS(opt$inputSeuratRds)

# Module code + report code
if(is.null(sobj@reductions$pca)){
  print(paste0("Cannot find PCA in the seurat object, a PCA with ", opt$ndims, " PCs is performed"))
  sobj <- RunPCA(sobj, npcs=opt$dims, verbose = F)
}
sobj <- RunUMAP(sobj,dims = 1:opt$ndims)

#pl2 <- DimPlot(sobj, reduction = "umap")

# Outputs writing into rds file
if(is.na(opt$output)){
  filename = paste0("SEURAT__DIM_REDUCTION_UMAP.rds")
  saveRDS(seuratobj.counts,file=filename, compress = T)
} else {
  saveRDS(seuratobj.counts,file=as.character(opt$output), compress = T)
}

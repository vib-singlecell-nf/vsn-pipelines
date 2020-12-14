#!/usr/bin/env Rscript

# Libraries loading
suppressPackageStartupMessages(library(reticulate))
use_python("/usr/bin/python3")
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(optparse))

# opt/arg parser definition
#####
option_list = list(
  make_option("--inputSeuratRds",
              default=NA,
              type='character',
              help="Path to a RDS file containing a Seurat object"),
  make_option("--output",
              default=NA,
              type='character',
              help="Output name"),
  make_option("--assay",
              default=NULL,
              type='character',
              help="Which assay to use for the tSNE."),
  make_option("--reduction",
              default="pca",
              type='character',
              help="Which dimensional reduction to use for the tSNE (e.g : pca, ica, dist)."),
  make_option("--dims",
              default=10,
              type='numeric',
              help="number of dimensions used to calculte the Tsne"),
  make_option("--seedUse",
              default=1,
              type='numeric',
              help="seed"),
  make_option("--tsneMethod",
              default="Rtsne",
              type='character',
              help="seed"),
  make_option("--addIter",
              default=0,
              type='numeric',
              help="If an existing tSNE has already been computed, uses the current tSNE to seed the algorithm and then adds additional iterations on top of this"),
  make_option("--dimEmbed",
              default=2,
              type='numeric',
              help="The dimensional space of the resulting tSNE"),
  make_option("--perplexity",
              default=30,
              type='numeric',
              help="The dimensional space of the resulting tSNE")
  )

#####
opt = parse_args(OptionParser(option_list=option_list))

# Input loading
seuratObj <- readRDS(opt$inputSeuratRds)

# Module code + report code
print(opt)
if(opt$reduction == "dist"){
  assay.data <- GetAssayData(seuratObj, assay = opt$assay, slot = "data")
  assay.dist <- dist(t(assay.data))
  
  seuratObj <- RunTSNE(seuratObj,
                  #cells = NULL,
                  #dims = 1:opt$dims,
                  #features = NULL,
                  perplexity = opt$perplexity,
                  seed.use = 1,
                  tsne.method = "Rtsne",
                  add.iter = 0,
                  dim.embed = 2,
                  distance.matrix = assay.dist,
                  reduction.name = paste0(opt$assay,"_tsne"),
                  reduction.key = paste0(tolower(opt$assay),"TSNE_"))
}else{
  if( opt$reduction == "pca" && !(paste0("pca_",opt$assay) %in% names(sobj@reductions))){
    print(paste0("Cannot find PCA in the seurat object, a PCA with ", opt$dims, " PCs is performed"))
    seuratObj <- RunPCA(seuratObj, npcs=opt$dims, verbose = F)
  }
  seuratObj <- RunTSNE(seuratObj,
                  assay = opt$assay,
                  reduction = opt$reduction,
                  #cells = NULL,
                  dims = 1:opt$dims,
                  #features = NULL,
                  perplexity = opt$perplexity,
                  seed.use = 1,
                  tsne.method = "Rtsne",
                  add.iter = 0,
                  dim.embed = 2,
                  reduction.name = paste0(opt$assay,"_tsne"),
                  reduction.key = paste0(tolower(opt$assay),"TSNE_"))
} 

# Outputs writing into rds file
if(is.na(opt$output)){
  filename = paste0("SEURAT__DIM_REDUCTION_TSNE.rds")
  saveRDS(seuratObj,file=filename, compress = T)
} else {
  saveRDS(seuratObj,file=as.character(opt$output), compress = T)
}

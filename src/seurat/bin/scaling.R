#!/usr/bin/env Rscript

# Libraries loading
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(optparse))

#opt/arg parser definition
#####
option_list = list(
  make_option("--inputSeuratRds",
              default=NA,
              type='character',
              help="Path to a RDS file containing a Seurat object"),
  make_option("--output",
              default=NULL,
              type='character',
              help="Output name"),
  make_option("--assay",
              default="RNA",
              type='character',
              help="Seurat assay to be processed"),
  make_option("--features",
              default=NULL,
              type='character',
              help="Vector of features names to scale/center. Default is variable features."),
  make_option("--varsToRegress",
              default=NULL,
              type='character',
              help="Variables to regress out (previously latent.vars in RegressOut). For example, nUMI, or percent.mito."),
  make_option("--splitBy",
              default=NULL,
              type='character',
              help="Name of variable in object metadata or a vector or factor defining grouping of cells."),
  make_option("--modelUse",
              default="linear",
              type='character',
              help="Used model for the regression. Options are 'linear' (default), 'poisson', and 'negbinom'"),
  make_option("--useUmi",
              default= FALSE,
              action="store_false",
              help="Regress on UMI count data, set to TRUE if modelUse is 'negbinom' or 'poisson'"),
  make_option("--doScale",
              default= TRUE,
              action="store_true",
              help="Whether to scale the data."),
  make_option("--doCenter",
              default= TRUE,
              action="store_true",
              help="Whether to center the data."),
  make_option("--scaleMax",
              default= 10,
              type='numeric',
              help="Max value to return for scaled data, The default is 10."),
  make_option("--blockSize",
              default= 1000,
              type='numeric',
              help="Default size for number of features to scale at in a single computation. Increasing blockSize may speed up calculations/increases memory cost."),
  make_option("--minCellsToBlock",
              default= 3000,
              type='numeric',
              help="If object contains fewer than this number of cells, don't block for scaling calculations.")
)
#####

opt = parse_args(OptionParser(option_list=option_list))

# Input loading<- 
sobj <- readRDS(opt$inputSeuratRds)

# Module code + report code
sobj <- ScaleData(sobj,
                  assay=opt$assay,
                  features = opt$features,
                  vars.to.regress = opt$varsToRegress,
                  split.by = opt$splitBy,
                  model.use = opt$modelUse,
                  use.umi = opt$useUmi,
                  do.scale = opt$doScale,
                  do.center = opt$doCenter,
                  scale.max = opt$scaleMax,
                  block.size = opt$blockSize,
                  min.cells.to.block = opt$minCellsToBlock,
                  verbose = FALSE)

# Outputs writing into rds file
if(is.null(opt$output)){
  filename = paste0(sobj@project.name,"_SEURAT__SCALING_",opt$assay,".rds")
  saveRDS(sobj,file=filename, compress = T)
} else {
  saveRDS(sobj,file=as.character(opt$output), compress = T)
}

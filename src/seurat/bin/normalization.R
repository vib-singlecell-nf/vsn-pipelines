#!/usr/bin/env Rscript

# Libraries loading
suppressPackageStartupMessages(library(reticulate))
use_python("/usr/bin/python3")
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

# opt/arg parser definition
option_list = list(
  make_option(c("-i", "--inputSeuratRds"),
              default=NA,
              type='character',
              help="Path to a RDS file containing a Seurat object"),
  make_option(c("-a", "--assay"),
              default="RNA",
              type='character',
              help="Seurat Assay to be processed"),
  make_option(c("-o", "--output"),
              default=NA,
              type='character',
              help="Output name"),
  make_option(c("-f", "--scalefactor"),
              default=10000,
              type='numeric',
              help="Scale factor value for the normalization"),
  make_option(c("-m", "--normalizationMethod"),
              default="LogNormalize",
              type='character',
              help="Normalization method"),
  make_option(c("-r", "--margin"),
              default=1,
              type='numeric',
              help="If performing CLR normalization, normalize across features (1) or cells (2)")
)

opt = parse_args(OptionParser(option_list=option_list))

# Input loading
sobj <- readRDS(opt$inputSeuratRds)

# Module code + report code

sobj <- NormalizeData(sobj,
                      assay = opt$assay,
                      normalization.method = opt$normalizationMethod,
                      scale.factor = opt$scalefactor,
                      margin = opt$margin)

# Outputs writing into rds file
if(is.na(opt$output)){
  filename = paste0("SEURAT__NORMALIZATION.rds")
  saveRDS(sobj,file=filename, compress = T)
} else {
  saveRDS(sobj,file=as.character(opt$output), compress = T)
}

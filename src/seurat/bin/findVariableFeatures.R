#!/usr/bin/env Rscript

# Libraries loading
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(patchwork))

# opt/arg parser definition
option_list = list(
  make_option("--seuratObj",
              default=NA,
              type='character',
              help="Path to a RDS file containing a Seurat object"),
  make_option("--output",
              default=NA,
              type='character',
              help="Output name"),
  make_option("--assay",
              default="RNA",
              type='character',
              help="Seurat assay to be processed"),
  make_option("--nfeatures",
              default=2000,
              type='numeric',
              help="Number of features to select as top variable features; only used when selectionMethod is set to 'dispersion' or 'vst'"),
  make_option("--selectionMethod",
              default="vst",
              type='character',
              help="Default is vst, other options are mean.var.plot and dispersion."),
  make_option("--loesSpan",
              default= 0.3,
              type='numeric',
              help="(vst method) Loess span parameter used when fitting the variance-mean relationship"),
  make_option("--clipMax",
              default="auto",
              type='character',
              help="(vst method) After standardization values larger than clip.max will be set to clip.max; default is 'auto' which sets this value to the square root of the number of cells"),
  make_option("--numBin",
              default= 20,
              type='numeric',
              help="Total number of bins to use in the scaled analysis (default is 20)"),
  make_option("--binningMethod",
              default="equal_width",
              type='character',
              help="Specifies how the bins should be computed. Default equal_width, other option : equal_frequency")
)

opt = parse_args(OptionParser(option_list=option_list))

# Input loading
sobj <- readRDS(opt$seuratObj)

if(nrow(sobj@assays[[opt$assay]]@counts) < opt$nfeatures){
  opt$nfeatures <- nrow(sobj@assays[[opt$assay]]@counts)
}
# Module code + report code
if(opt$assay != "SCT"){
  sobj <- FindVariableFeatures(sobj,
                               assay = opt$assay,
                               selection.method = opt$selectionMethod,
                               nfeatures = opt$nfeatures,
                               loess.span = opt$loesSpan,
                               clip.max = opt$clipMax,
                               num.bin = opt$numBin,
                               binning.method = opt$binningMethod)
}

# Identify the 10 most highly variable genes
sobj@misc[[paste0(opt$assay,"_HVG")]] <- VariableFeatures(sobj,assay= opt$assay)
top10 <- head(VariableFeatures(sobj,assay = opt$assay), 10)
plot1 <- VariableFeaturePlot(sobj, assay = opt$assay)
plot2 <- LabelPoints(plot = plot1, points = top10,label = gsub("(-{2,}|_{1,})","-",top10), repel = TRUE)

dir.create(paste0("Plots/",opt$assay), recursive = T)
png(file=paste0("Plots/",opt$assay,"/05_hvg.png"), width = 850, height = 642, type = 'cairo')
plot1+ plot2
dev.off()

sobj@tools$diagnostics[[paste0(opt$assay,'_varGenes')]] <- length(VariableFeatures(sobj))

# Outputs writing into rds file
if(is.na(opt$output)){
  filename = paste0(sobj@project.name,"SEURAT__FIND_VARIABLE_FEATURES_",opt$assay,".rds")
  saveRDS(sobj,file=filename, compress = T)
} else {
  saveRDS(sobj,file=as.character(opt$output), compress = T)
}

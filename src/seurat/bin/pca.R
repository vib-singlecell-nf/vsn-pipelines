#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))

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
              help="assay name"),
  make_option("--nPcs",
              default=50,
              type='numeric',
              help="number of calculated PCs"),
  make_option("--nPlotedPcs",
              default=2,
              type='numeric',
              help="number of ploted PCs"),
  make_option("--scaleData",
              default = FALSE,
              action = "store_true",
              help= "Scale data before the PCA "),
  make_option("--removeScaledData",
              default = FALSE,
              action = "store_true",
              help= "Remove Scaled data (use for RNA data). If TRUE, will speed up the pipeline and save RAM "),
  make_option("--excludePatternHVG",
              default = NA,
              type = "character",
              help = "Remove features from HVG matching this pattern"),
  make_option("--diagnosticPlots",
              default = FALSE,
              action = "store_true",
              help = "If true, generates DimPlots, ElbowPlot, HeatMap, ... ")
)

opt = parse_args(OptionParser(option_list=option_list))

# Input loading
seuratObj <- readRDS(opt$seuratObj)

# Module code + report code
if(opt$scaleData){
  seuratObj <- ScaleData(seuratObj, assay = opt$assay, verbose=F)
  seuratObj <- FindVariableFeatures(seuratObj, assay = opt$assay)
}

if( !is.na(opt$excludePatternHVG) ){
  hvg <- VariableFeatures(seuratObj, assay = opt$assay)[ !(grepl(opt$excludePatternHVG, VariableFeatures(seuratObj, assay = opt$assay)) )]
} else {
  hvg <- VariableFeatures(seuratObj, assay = opt$assay)
}
reduc.name <- paste0(opt$assay,"_pca")

seuratObj <- RunPCA(seuratObj,
               npcs=opt$nPcs,
               assay = opt$assay,
               reduction.name = reduc.name,
               reduction.key = paste0(tolower(opt$assay),"PCA_"),
               features = hvg,
               verbose = F)


if(opt$diagnosticPlots){
  dir.create(paste0("Plots/",opt$assay),recursive = T)
  dim.comb <- combn(1:opt$nPlotedPcs,2)
  pdf(file = paste0("Plots/",opt$assay,"/07_PCA_",opt$assay,".pdf"), width = 10)
  for(i in 1:ncol(dim.comb)){
    print(DimPlot(object = seuratObj, reduction = reduc.name, dims = dim.comb[,i], group.by= "subsample"))
  }
  dev.off()

  pdf(file=paste0("Plots/",opt$assay,"/07b_PCA_diag_",opt$assay,".pdf"), width = 10)
  print(ElbowPlot(seuratObj ,reduction = reduc.name))
  print(VizDimLoadings(seuratObj, dims = 1:opt$nPlotedPcs, reduction = reduc.name))
  print(DimHeatmap(seuratObj,assays = opt$assay ,dims = 1:opt$nPlotedPcs, reduction = reduc.name, fast=FALSE))
  dev.off()

  red.nb <- which(names(seuratObj@reductions) == reduc.name)
  assay.nb <- which(names(seuratObj@assays) == opt$assay)
  loadings <- rownames_to_column(as.data.frame(seuratObj@reductions[[red.nb]]@feature.loadings),var="feature")
  var.load <- seuratObj@assays[[assay.nb]]@meta.features %>%
    rownames_to_column(var="feature") %>%
    dplyr::filter(feature %in% hvg) %>%
    left_join(loadings,by="feature")

  print("Make PC plots for inspection")

  pc_list <- colnames(seuratObj@reductions[[red.nb]]@feature.loadings)
  pc_plots <- list()
  for (i in 1:50) {
    labels <-ifelse(abs(var.load[,2+i])>0.1,as.character(var.load$feature),'')
    xplot <- ifelse(opt$assay == "SCT", "sct.residual_variance", grep(".variance.standardized",colnames(var.load),value = T))
    p <- ggplot(var.load,aes_string(xplot,pc_list[[i]][1],color=pc_list[[i]][1])) +
      geom_point() +
      geom_abline(aes(slope = 0, intercept = 0)) +
      geom_label_repel(aes(label=ifelse(abs(.data[[pc_list[[i]]]])>0.1,as.character(feature),'')),color="black") +
      theme_classic() +
      scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red" )
    pc_plots[[i]] <- p
    print(i)
  }

  pdf(file = paste0("Plots/",opt$assay,"/08a_PCloadings_",opt$assay,".pdf"))
  for (i in 1:50) {
    print(pc_plots[[i]])
  }
  dev.off()

}

if(opt$removeScaledData){
  seuratObj@assays[[which(names(seuratObj@assays) == opt$assay)]]@scale.data <- matrix()
}

# Outputs writing into rds file
if(is.na(opt$output)){
  filename = paste0(seuratObj@project.name,"SEURAT__DIM_REDUCTION_PCA_",opt$assay,".rds")
  saveRDS(seuratObj,file=filename, compress = T)
} else {
  saveRDS(seuratObj,file=as.character(opt$output), compress = T)
}

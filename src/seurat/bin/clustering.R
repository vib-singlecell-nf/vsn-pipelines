suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(intrinsicDimension))
suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(ggrepel))

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
    help="Output name"),
  make_option(
    "--assay",
    default=NA,
    type='character',
    help="assay name"),
  make_option(
    "--dimsToUse",
    default=NA,
    type = "numeric",
    help= "Number of PCs on which to run find neighbors, TSNE and UMAP"
  ),
  make_option(
    "--resToUse",
    default=NA,
    type = "numeric",
    help= "Clustering's resolution"
  ),
  make_option(
    "--perplexity",
    default=NA,
    type = "numeric",
    help= "Perplexity"
  ),
  make_option(
    "--assayForCrossModalityGraphs",
    default=NA,
    type = "character",
    help= "Assay name for cross modality UMAP and TSNE"
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

 # opt <- list(seuratObj = "datasets/SAR13/SAR13.SEURAT__DIMENSIONALITY_REDUCTION_PCA_SCT.rds",
 #             assay = "SCT",
 #             resToUse=NA,
 #             dimsToUse=NA,
 #             perplexity=NA)

seuratObj <- readRDS(opt$seuratObj)

dir.create(paste0("Plots/",opt$assay), recursive = T)

reduc.name <- paste0(opt$assay,"_pca")
red.nb <- which(names(seuratObj@reductions) == reduc.name)
#assay.nb <- which(names(seuratObj@assays) == opt$assay)

####################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
####################################################

# based on Pipecomp paer
int.dim <- maxLikGlobalDimEst(seuratObj@reductions[[red.nb]]@cell.embeddings,
                              k=20,unbiased = TRUE, neighborhood.aggregation = 'robust')

est.PC <-round(int.dim[[1]])
seuratObj@tools$diagnostics[[paste0('est.dimsPC.maxlik.',opt$assay)]] <- est.PC

### Create PC Elbowplot
png(file = paste0("Plots/",opt$assay,"/08b_selectPC_",opt$assay,".png"), width = 850, height = 642, type = 'cairo')
ElbowPlot(object = seuratObj, ndims = 50, reduction = reduc.name) + geom_vline(aes(xintercept=est.PC))
dev.off()

if(is.na(opt$dimsToUse)){
  maxdimpca <- ifelse(dim(seuratObj@reductions[[red.nb]])[2]<40,dim(seuratObj@reductions[[red.nb]])[2],40)
  dimsToUse<-c(est.PC,seq(10,maxdimpca,by=10))
  seuratObj@tools$diagnostics[[paste0('dimsPC.',opt$assay)]] <- seuratObj@tools$diagnostics[[paste0('est.dimsPC.maxlik.',opt$assay)]]
}else{
  seuratObj@tools$diagnostics[[paste0('dimsPC.',opt$assay)]] <- opt$dimsToUse
  dimsToUse <- opt$dimsToUse
}

for(maxPCs in dimsToUse){
  
  print(paste0("Working on 1:",maxPCs))

  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, assay = opt$assay, reduction = reduc.name, dims = 1:maxPCs, verbose = F)
  if(is.na(opt$resToUse)){
    seuratObj <- FindClusters(object = seuratObj,  graph.name = paste0(opt$assay,"_snn"), resolution = 0.8, verbose = F)
  }else{
    seuratObj <- FindClusters(object = seuratObj,  graph.name = paste0(opt$assay,"_snn"), resolution = opt$resToUse, verbose = F)
  }
  
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj,
                               dims = 1:maxPCs,
                               assay = opt$assay,
                               reduction = reduc.name,
                               reduction.name = paste0(opt$assay,"_tsne"),
                               reduction.key = paste0(tolower(opt$assay),"TSNE_"),
                               verbose = F)
  
  tsnePlot <- DimPlot(seuratObj, reduction = paste0(opt$assay,"_tsne"), label=T, label.size = 8)
  tsnePlotSplit <- DimPlot(seuratObj, reduction = paste0(opt$assay,"_tsne"), label=F, group.by="ident", pt.size = 2)
  
  ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
         file=paste0("Plots/",opt$assay,"/10a_tSNE_1-",maxPCs,"_",opt$assay,".png"),
         width = 20, height=10,
         type = "cairo")
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj,
                       dims = 1:maxPCs,
                       assay = opt$assay,
                       reduction = reduc.name,
                       reduction.name = paste0(opt$assay,"_umap"),
                       reduction.key = paste0(tolower(opt$assay),"UMAP_"),
                       verbose = F)
  
  umapPlot <- DimPlot(seuratObj, reduction = paste0(opt$assay,"_umap"), label = T, label.size = 8)
  umapPlotSplit <- DimPlot(seuratObj, reduction = paste0(opt$assay,"_umap"), label = F, group.by="ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0("Plots/",opt$assay,"/10b_UMAP_1-",maxPCs,"_",opt$assay,".png"),
         width = 20, height=10,
         type = "cairo")
  
}

### Clustering: trying out clusTree

Perplexity <- ifelse(is.na(opt$perplexity),
                     ifelse(is.na(opt$dimsToUse),est.PC,opt$dimsToUse),
                     opt$perplexity)

seuratObj <- FindNeighbors(object = seuratObj, reduction = reduc.name, dims = 1:Perplexity)

if(is.na(opt$resToUse)){
  # Resolution estimation
  resolutions <- seq(0,1,by=0.1)
  seuratObj <- FindClusters(object = seuratObj,  graph.name = paste0(opt$assay,"_snn"), resolution = resolutions, verbose = F)
  
  pdf(file = paste0("Plots/",opt$assay,"/10c_Clustree_",opt$assay,".pdf"))
  print(clustree(seuratObj, prefix = paste0(opt$assay,"_snn_res.")))
  dev.off()
  
  # Use default value
  seuratObj@tools$diagnostics[[paste0('res.',opt$assay)]] <- "default"
  
}else{
  seuratObj@tools$diagnostics[[paste0('res.',opt$assay)]] <- opt$resToUse
  seuratObj <- FindClusters(object = seuratObj,  graph.name = paste0(opt$assay,"_snn"), resolution = opt$resToUse, verbose = F)
}

seuratObj <- RunTSNE(object = seuratObj,
                     dims = 1:Perplexity,
                     assay = opt$assay,
                     reduction = reduc.name,
                     reduction.name = paste0(opt$assay,"_tsne"),
                     reduction.key = paste0(tolower(opt$assay),"TSNE_"),
                     verbose = F)

seuratObj <- RunUMAP(seuratObj,
                     dims = 1:Perplexity,
                     assay = opt$assay,
                     reduction = reduc.name,
                     reduction.name = paste0(opt$assay,"_umap"),
                     reduction.key = paste0(tolower(opt$assay),"UMAP_"),
                     verbose = F)

if(is.na(opt$resToUse)){
  umapPlot<-DimPlot(seuratObj, reduction = paste0(opt$assay,"_umap"), label = T, group.by= paste0(opt$assay,"_snn_res.0.8"), label.size = 6)
  tsnePlot<-TSNEPlot(seuratObj, reduction = paste0(opt$assay,"_tsne"), label = T, group.by= paste0(opt$assay,"_snn_res.0.8"), label.size = 6)
  #seuratObj@active.assay
  seuratObj@meta.data[paste0(opt$assay,"_clusters")] <- seuratObj@meta.data[,which(colnames(seuratObj@meta.data) == paste0(opt$assay,"_snn_res.0.8"))]
  
  pdf(file = paste0("Plots/",opt$assay,"/11_tSNE_UMAP_default_resolution_",opt$assay,".pdf"), width = 17*0.45, height = 12.4*0.45)
  print(umapPlot)
  print(tsnePlot)
  dev.off()
} else{
  umapPlot <- DimPlot(seuratObj, reduction = paste0(opt$assay,"_umap"), label = T, group.by= paste0(opt$assay,"_snn_res.",opt$resToUse), label.size = 6)
  tsnePlot <- TSNEPlot(seuratObj, reduction = paste0(opt$assay,"_tsne"), label = T, group.by= paste0(opt$assay,"_snn_res.",opt$resToUse), label.size = 6)
  #seuratObj@active.assay
  seuratObj@meta.data[paste0(opt$assay,"_clusters")] <- seuratObj@meta.data[,colnames(seuratObj@meta.data) == paste0(opt$assay,"_snn_res.",opt$resToUse)]
  pdf(file=paste0("Plots/",opt$assay,"/11_tSNE_UMAP_",opt$assay,".pdf"), width = 17*0.45, height = 12.4*0.45)
  print(umapPlot)
  print(tsnePlot)
  dev.off()
}

if(!is.na(opt$assayForCrossModalityGraphs)){
  ###################################
  #### Cross-Modality UMAPS and tSNE
  ##################################
  umapAssay <- paste0(opt$assay,"_umap")
  tsneAssay <- paste0(opt$assay,"_tsne")
  clusterAssay <- paste0(opt$assay,"_clusters")
  umapCross <- paste0(opt$assayForCrossModalityGraphs,"_umap")
  tsneCross <- paste0(opt$assayForCrossModalityGraphs,"_tsne")
  clusterCross <- paste0(opt$assayForCrossModalityGraphs,"_clusters")
  
  umap_1_1 <- DimPlot(seuratObj, reduction = umapAssay, group.by=clusterAssay, label=T)  + ggtitle(paste0(opt$assay," clusters"))
  umap_1_2 <- DimPlot(seuratObj, reduction = umapAssay, group.by=clusterCross, label=T)  + ggtitle(paste0(opt$assayForCrossModalityGraphs," clusters"))
  umap_2_1 <- DimPlot(seuratObj, reduction = umapCross, group.by=clusterAssay, label=T)  + ggtitle(paste0(opt$assay," clusters"))
  umap_2_2 <- DimPlot(seuratObj, reduction = umapCross, group.by=clusterCross, label=T)  + ggtitle(paste0(opt$assayForCrossModalityGraphs," clusters"))
  
  tsne_1_1 <- DimPlot(seuratObj, reduction = tsneAssay, group.by=clusterAssay, label=T)  + ggtitle(paste0(opt$assay," clusters"))
  tsne_1_2 <- DimPlot(seuratObj, reduction = tsneAssay, group.by=clusterCross, label=T)  + ggtitle(paste0(opt$assayForCrossModalityGraphs," clusters"))
  tsne_2_1 <- DimPlot(seuratObj, reduction = tsneCross, group.by=clusterAssay, label=T)  + ggtitle(paste0(opt$assay," clusters"))
  tsne_2_2 <- DimPlot(seuratObj, reduction = tsneCross, group.by=clusterCross, label=T)  + ggtitle(paste0(opt$assayForCrossModalityGraphs," clusters"))
  
  pdf(file=paste0("Plots/",opt$assay,"/06_overview_UMAP_tSNE.pdf"), width = 17, height = 12.4)
  (umap_1_1 + umap_1_2)/(umap_2_1 + umap_2_2)
  (tsne_1_1 + tsne_1_2)/(tsne_2_1 + tsne_2_2)
  dev.off()
}


# Outputs writing into rds file
if(is.na(opt$output)){
  filename = paste0(seuratObj@project.name,"SEURAT__CLUSTERING_",opt$assay,".rds")
  saveRDS(seuratObj,file=filename, compress = T)
} else {
  saveRDS(seuratObj,file=as.character(opt$output), compress = T)
}
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library('DoubletFinder'))
suppressPackageStartupMessages(library('fields'))


option_list = list(
  make_option(
    "--seuratObj",
    default=NA,
    type='character',
    help="File path to the rds file containing a Seurat object."
  ),
  make_option(
    "--minPCT",
    default=1,
    type = "numeric",
    help= "Min number of PCs. Uses to calculate the pANN threshold used to make final doublet/singlet predictions."
  ),
  make_option(
    "--maxPCT",
    default=15,
    type = "numeric",
    help= "Max number of PCs. Uses to calculate the pANN threshold used to make final doublet/singlet predictions."
  ),
  make_option(
    "--pN",
    default=0.25,
    type = "numeric",
    help= "Defines the number of generated artificial doublets, expressed as a proportion of the merged real-artificial data"
  ),
  make_option(
    "--dimsToUse",
    default=NA,
    type = "numeric",
    help= "Number of PCs to run doubletFinder (optional if RNA_SCT_CLUSTERING_TSNE_UMAP was run before)"
  ),
  make_option(
    "--cores",
    default=1,
    type = "numeric",
    help= "Number of cores to run paramSweep_v3"
  ),
  make_option(
    "--output",
    default=NA,
    type='character',
    help="output file name"
  )
)

opt <- parse_args(OptionParser(option_list=option_list))


runDoubletFinder <- function(object, PCs, minPCT = 1, maxPCT = 10, pN = 0.25){

  DFPredictions <- c()

  # we need to add that we used sct and do not have a ground truth
  print("running paramSweep")
  findpK <- paramSweep_v3(object, PCs = 1:PCs, sct=T, num.cores=opt$cores) %>%
  summarizeSweep(GT=F) %>%
  find.pK()


  maxScore <- findpK %>% pull('BCmetric') %>% which.max()
  pKValue <- findpK[maxScore, 'pK'] %>% as.character() %>% as.numeric()

  # correct doublet ratio on class proportions
  homotypic.prop <- modelHomotypic(object@meta.data$SCT_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi.min <- round((minPCT/100)*length(colnames(object)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.min.adj <- round(nExp_poi.min*(1-homotypic.prop))
  nExp_poi.max <- round((maxPCT/100)*length(colnames(object)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.max.adj <- round(nExp_poi.max*(1-homotypic.prop))

  # running the 4 runs
  print("1st run")
  object <- doubletFinder_v3(seu = object,
                             pN = pN,
                             pK = pKValue,
                             nExp =  nExp_poi.min,
                             reuse.pANN = F,
                             PCs = 1:PCs,
                             sct = T)
  print("2nd run")
  object <- doubletFinder_v3(seu = object,
                             pN = pN,
                             pK = pKValue,
                             nExp =  nExp_poi.min.adj,
                             reuse.pANN = paste0("pANN_", pN, "_", pKValue, "_", nExp_poi.min),
                             PCs = 1:PCs,
                             sct = T)
  print("3rd run")
  object <- doubletFinder_v3(seu = object,
                             pN = pN,
                             pK = pKValue,
                             nExp =  nExp_poi.max,
                             reuse.pANN = paste0("pANN_", pN, "_", pKValue, "_",  nExp_poi.min),
                             PCs = 1:PCs,
                             sct = T)
  print("4th run")
  object <- doubletFinder_v3(seu = object,
                             pN = pN,
                             pK = pKValue,
                             nExp = nExp_poi.max.adj,
                             reuse.pANN = paste0("pANN_", pN, "_", pKValue, "_",  nExp_poi.min),
                             PCs = 1:PCs,
                             sct = T)

  object@meta.data$DFPrediction <- object@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_", nExp_poi.max)]
  object@meta.data$DFPrediction[ object@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_",  nExp_poi.max.adj)] == 'Doublet'] <- "Low Confidence adjusted"
  object@meta.data$DFPrediction[ object@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_",  nExp_poi.min)] == 'Doublet'] <- "High Confidence"
  object@meta.data$DFPrediction[ object@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_",  nExp_poi.min.adj)] == 'Doublet'] <- "High Confidence adjusted"
  object@meta.data$DFPrediction <- gsub('Doublet', "Low Confidence", object@meta.data$DFPrediction)

  return(object)
}

#################################################
### Doublet detection: doubletfinder
#################################################
## code was updated to accommodate SCT
## and Homotypic Doublet Proportion Estimate and adjustment
seuratObj <- readRDS(opt$seuratObj)
# test
#seuratObj <- readRDS("../work/df/b2fd54f7707e7584f27b231a7b8f82/seuratObj_diagPlots_SCT.rds")

dir.create("Plots/RNA", recursive = T)

if(is.na(opt$dimsToUse)){
  if(!is.null(seuratObj@tools$diagnostics[["dimsPC.SCT"]])){
    PCs <- seuratObj@tools$diagnostics[["dimsPC.SCT"]]
  }else if(!is.null(seuratObj@tools$diagnostics[["dimsPC.RNA"]])){
    PCs <- seuratObj@tools$diagnostics[["dimsPC.RNA"]]
  }else{
    stop("A fully pre-processed seurat object is needed to run DoubletFinder")
  }
}else{
  PCs <- opt$dimsToUse
}

seuratObj <- runDoubletFinder(seuratObj,
                              PCs = PCs,
                              minPCT = opt$minPCT,
                              maxPCT = opt$maxPCT,
                              pN = opt$pN)

# DimPlot(seuratObj, reduction = 'SCT_umap', group.by = 'DFPrediction')
# DimPlot(seuratObj, reduction = 'SCT_tsne', group.by = 'DFPrediction')

# additional visualization would be stacked bar plot
#dblf.plot <-seuratObj@meta.data

#### THIS IS NOT WORKING !!!! ####
# doublet.fill <- dplyr::select(dblf.plot,DFPrediction,SCT_clusters) %>%
#   group_by(DFPrediction) %>%
#   count() %>%
#   ggplot(aes(fill = DFPrediction, y = frequency(), x = n)) +
#   geom_bar(position="fill", stat="identity") +
#   theme_classic() +
#   labs(title="DoubletFinder prediction", x="Clusters", y="%")
#
# doublet.stack <- dplyr::select(dblf.plot,DFPrediction,SCT_clusters) %>%
#   group_by(DFPrediction) %>%
#   count() %>%
#   ggplot(aes(fill = DFPrediction,y=frequency() , x = n)) +
#   geom_bar(position="stack", stat="identity") +
#   theme_classic() +
#   labs(title="DoubletFinder prediction", x="Clusters", y="%")

png(filename = "Plots/RNA/16_DoubletFinder.png", width = 500, height=500, type= 'cairo')
  print(DimPlot(seuratObj, reduction = 'SCT_umap', group.by = 'DFPrediction'))
  #print(doublet.stack)
  #print(doublet.fill)
dev.off()

# Outputs writing into rds file
if(is.na(opt$output)){
  filename = paste0(seuratObj@project.name,".DOUBLETFINDER.rds")
  saveRDS(seuratObj,file=filename, compress = T)
} else {
  saveRDS(seuratObj,file=as.character(opt$output), compress = T)
}

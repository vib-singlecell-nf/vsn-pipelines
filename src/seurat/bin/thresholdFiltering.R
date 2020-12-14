#!/usr/bin/env Rscript

# Libraries loading
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))

# opt/arg parser definition
option_list = list(
  make_option("--seuratObj",
              default=NA,
              type='character',
              help="Path to a RDS file containing a Seurat object"),
  make_option("--output",
              default=NA,
              type='character',
              help="Output file name"),
  make_option(
    "--nmad_low_feature",
    default=5,
    type='numeric',
    help="lower bound of the minimum number of MADs away from median for feature filtering"
  ),
  make_option(
    "--nmad_high_feature",
    default=5,
    type='numeric',
    help="upper bound of the minimum number of MADs away from median for feature filtering"
  ),
  make_option(
    "--nmad_low_UMI",
    default=5,
    type='numeric',
    help="lower bound of the minimum number of MADs away from median for UMI filtering"
  ),make_option(
    "--nmad_high_UMI",
    default=5,
    type='numeric',
    help="upper bound of the minimum number of MADs away from median for UMI filtering"
  ),make_option(
    "--nmad_high_mito",
    default=NA,
    type='numeric',
    help="upper bound of the minimum number of MADs away from median for Mitochondrial genes filtering"
  ),
  make_option(
    "--scriptFunctions",
    default=NA,
    type='character',
    help="path of the script_funtions_covid.R file"
  )
)

opt = parse_args(OptionParser(option_list=option_list))

source(opt$scriptFunctions)
# Input loading
seuratObj <- readRDS(opt$seuratObj)

metaData <- seuratObj@meta.data
diagnostics <- seuratObj@tools$diagnostics

if("batch" %in% colnames(seuratObj@meta.data)){
  batch <- seuratObj$batch
}else{
  if("subsample" %in% colnames(seuratObj@meta.data)){
    seuratObj<- subset(seuratObj, cells = colnames(seuratObj)[grepl("#",seuratObj$subsample)])
    batch <- seuratObj$subsample
  }else{
    batch <- rep(seuratObj@project.name,ncol(seuratObj))
  }
}

feature.drop.low <- is.Outlier(metaData$nGene, nmads=opt$nmad_low_feature, type="lower", log=TRUE,batch = batch)
feature.drop.high <- is.Outlier(metaData$nGene, nmads=opt$nmad_high_feature, type="higher", log=TRUE,batch = batch)
metaData$nGene.drop <- as.logical(feature.drop.low + feature.drop.high,batch = batch)
##same as UMI in Seurat pipeline
libsize.drop.low <- is.Outlier(metaData$nUMI, nmads=opt$nmad_low_UMI, type="lower", log=TRUE,batch = batch)
libsize.drop.high <- is.Outlier(metaData$nUMI, nmads=opt$nmad_high_UMI, type="higher", log=TRUE,batch = batch)
metaData$nUMI.drop <- as.logical(libsize.drop.low+libsize.drop.high)

##### add to metaData matrix #####

metaData$final.drop = metaData$nGene.drop | metaData$nUMI.drop

##### Add to diagnostics #####
diagnostics[['nmad.low.feature']] <- opt$nmad_low_feature
diagnostics[['nmad.high.feature']] <- opt$nmad_high_feature
diagnostics[['nmad.low.libsize']] <- opt$nmad_low_UMI
diagnostics[['nmad.high.libsize']] <- opt$nmad_high_UMI

diagnostics[['feature.drop.low']]<-sum(feature.drop.low)
diagnostics[['feature.drop.high']]<-sum(feature.drop.high)
diagnostics[['feature.drop']] <- sum(metaData$nGene.drop)
diagnostics[['libsize.drop.low']] <- sum(libsize.drop.low)
diagnostics[['libsize.drop.high']] <- sum(libsize.drop.high)
diagnostics[['libsize.drop']] <- sum(metaData$nUMI.drop)
varCols <- c('nGene.drop','nUMI.drop')
varPlots <- c("nGene","nUMI")

if(!is.na(opt$nmad_high_mito)){
  ## mitochondrial genes
  mito.drop.high <- is.Outlier(metaData$percent.mito, nmads=opt$nmad_high_mito, type="higher",batch = batch)
  metaData$mito.drop <- as.logical(mito.drop.high)
  metaData$final.drop <- metaData$final.drop | metaData$mito.drop
  diagnostics[['nmad.high.mito']] <- opt$nmad_high_mito
  diagnostics[['mito.drop.high']]<-sum(mito.drop.high)
  diagnostics[['mito.drop']]<-sum(metaData$mito.drop)
  varCols <- c(varCols, "mito.drop")
  varPlots <- c(varPlots, "percent.mito")
}

########################################
########## Create violinPlots
########################################
dir.create("Plots/RNA",recursive = T)

drawVlnPlot_out(metaData,
                fileName = "Plots/RNA/02a_Outliers_nGene.png",
                colsToColor = list(nGene = 'nGene.drop',
                                   nUMI = 'nGene.drop',
                                   percent.mito = 'nGene.drop',
                                   percent.rbc = 'nGene.drop',
                                   percent.COVID = 'nGene.drop'),
                png.device.type = NULL)

drawVlnPlot_out(metaData,
                fileName = "Plots/RNA/02a_Outliers_nUMI.png",
                colsToColor = list(nGene = 'nUMI.drop',
                                   nUMI = 'nUMI.drop',
                                   percent.mito = 'nUMI.drop',
                                   percent.rbc = 'nUMI.drop',
                                   percent.COVID = 'nUMI.drop'),
                png.device.type = "cairo")
drawVlnPlot_out(metaData,
                fileName = "Plots/RNA/02a_Outliers_percMito.png",
                colsToColor = list(nGene = 'mito.drop',
                                   nUMI = 'mito.drop',
                                   percent.mito = 'mito.drop',
                                   percent.rbc = 'mito.drop',
                                   percent.COVID = 'mito.drop'),
                png.device.type = "cairo")

#filename can be added to write the output to png file
drawVlnPlot_color(metaData,"Plots/RNA/02a_Outliers_colors_base.png",
                  png.device.type = "cairo")
if("mito.drop" %in% colnames(metaData)){
  drawVlnPlot_color_nGene(metaData,"Plots/RNA/02a_Outliers_colors_nGene.png",
                          png.device.type = "cairo")
  drawVlnPlot_color_nGene_mito(metaData,"Plots/RNA/02a_Outliers_colors_nGene&mito.png",
                               png.device.type = "cairo")
}
if("percent.COVID" %in% colnames(metaData) && "percent.rbc" %in% colnames(metaData)){
  drawVlnPlot_color_mito_extra(metaData,"Plots/RNA/02a_Outliers_colors_Extra_mito.png",
                               png.device.type = "cairo")
  drawVlnPlot_color_nGene_extra(metaData,"Plots/RNA/02a_Outliers_colors_Extra_nGene.png",
                                png.device.type = "cairo")
  drawVlnPlot_color_nUMI_extra(metaData,"Plots/RNA/02a_Outliers_colors_Extra_nUMI.png",
                               png.device.type = "cairo")

  varCols <- c(varCols, "final.drop","final.drop")
  varPlots <- c(varPlots, "percent.rbc","percent.COVID")
}else if("percent.rbc" %in% colnames(metaData)){
  varCols <- c(varCols, "final.drop")
  varPlots <- c(varPlots, "percent.rbc")
}


### Before filtering
drawVlnPlot(metaData, fileName = "Plots/RNA/02b_beforeFiltering.png",
            varsToPlot = varPlots,
            colsToColor = varCols,
            png.device.type = "cairo")

### After filtering
#toPlot <- metaData[! metaData$final.drop,]
drawVlnPlot(metaData[! metaData$final.drop,], fileName = "Plots/RNA/02c_afterFiltering.png",
            varsToPlot = varPlots,
            colsToColor = varCols,
            png.device.type = "cairo")

####################################
#### Filtering
###################################

seuratObj <- seuratObj[,!(metaData$nUMI.drop | metaData$nGene.drop | metaData$mito.drop)]

### Number of cells removed
diagnostics[['firstRemove']]<-nrow(metaData)-ncol(seuratObj)
diagnostics[['dimFirstRemove']]<-paste0(nrow(seuratObj)," genes - ",ncol(seuratObj)," cells")

seuratObj@tools$diagnostics <-c(seuratObj@tools$diagnostics, diagnostics)
seuratObj@meta.data <- cbind(seuratObj@meta.data,metaData[!(metaData$nUMI.drop | metaData$nGene.drop | metaData$mito.drop),
                                                          which(!colnames(metaData) %in% colnames(seuratObj@meta.data))])



# Outputs writing into rds file
if(is.na(opt$output)){
  filename = paste0(seuratObj@project.name,".SEURAT__THRESHOLDING_FILTERING.rds")
  saveRDS(seuratObj,file=filename, compress = T)
} else {
  saveRDS(seuratObj,file=as.character(opt$output), compress = T)
}

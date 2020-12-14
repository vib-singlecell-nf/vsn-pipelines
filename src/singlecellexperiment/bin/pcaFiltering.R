suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(gridExtra))

option_list = list(
  make_option(
    "--sceObj",
    default=NA,
    type='character',
    help="File path to the rds file containing a SCE object."
  ),
  make_option(
    "--removeOutliers",
    default=FALSE,
    action = "store_true",
    help="Remove cells detected as outliers by the PCA"
  ),
  make_option(
	"--scriptFunctions",
	default=NA,
	type='character',
	help="path of the script_funtions_BioIT.R file"
  ),
  make_option(
    "--output",
    default=NA,
    type='character',
    help="output file name"
  )

)
opt <- parse_args(OptionParser(option_list=option_list))

source(opt$scriptFunctions)
sce <- readRDS(opt$sceObj)
metaData <- sce@colData
dir.create("Plots/RNA" ,recursive = T)
########################################
########## Create PCA
########################################

##(check via raw code of the function runPCA)
varsToUse <- c("nUMI","nGene", "percent.mito")
#setdiff(colnames(colData(sce)),varsToUse)
exprs_to_plot <- scale(colData(sce)[,varsToUse], scale = T)
x.mad <- apply(exprs_to_plot, 2, mad)
varsToUse<-setdiff(varsToUse, names(x.mad[x.mad==0]))
sceNew<-runColDataPCA(sce, outliers=T, variables = varsToUse)

#### Detect bad cells #####
#sceNew<-runPCA(sce,use_coldata=T, detect_outliers=T)
# table(sceNew$outlier)

outs <-colnames(sceNew)[sceNew$outlier]
### Add to metaData
metaData$pca.drop <- metaData$final.drop
metaData[outs,which(colnames(metaData)=="pca.drop")]<-TRUE
sce@colData$pca.drop <- metaData$pca.drop 
##### Color bad cells on PCA plot #####
colorDF<-as.data.frame(cbind(colnames(sceNew),"1"), stringsAsFactors=F)
rownames(colorDF)<-colorDF[,1]
colorDF[outs,2]<-"2"
colorDF[,2]<-as.factor(colorDF[,2])
tmp<-colorDF[,2,drop=F]

print("Plots")
png(file= "Plots/RNA/03a_PCA.png",  width = 850, height = 642, type = 'cairo')
plotReducedDim(sceNew, dimred = "PCA_coldata", colour_by='outlier',shape_by='outlier') + labs(title="PCA with outliers colored")
dev.off()

#### Add to metaData table ####
pca.drop<-metaData[colnames(sce),"pca.drop"]
# sum(pca.drop)

##### Create violinplots ####
##Before
metatemp <- as.data.frame(metaData[!metaData$final.drop,])
drawVlnPlot_out(metatemp, 
                fileName = "Plots/RNA/03b_beforePcaFiltering.png", 
                colsToColor = list(nGene = 'pca.drop',
                                   nUMI = 'pca.drop',
                                   percent.mito = 'pca.drop',
                                   percent.rbc = 'pca.drop',
                                   percent.COVID = 'pca.drop'),
                png.device.type = "cairo")

##After
metatemp <- as.data.frame(metaData[! metaData$pca.drop,])
drawVlnPlot_out(metatemp, 
                fileName = "Plots/RNA/03c_afterPcaFiltering.png", 
                colsToColor = list(nGene = 'pca.drop',
                                   nUMI = 'pca.drop',
                                   percent.mito = 'pca.drop',
                                   percent.rbc = 'pca.drop',
                                   percent.COVID = 'pca.drop'),
                png.device.type = "cairo")

if(opt$removeOutliers){
  sce@metadata$tools$diagnostics[['pcaRemove']]<-sum(pca.drop)
  sce@metadata$tools$diagnostics[['totalRemove']]<-nrow(metaData)-ncol(sce)
  ##### Remove outlier cells #####
  sce <- sce[,!(pca.drop)]
  sce@metadata$tools$diagnostics[['dimAfterPCA']]<-paste0(nrow(sce)," genes - ",ncol(sce)," cells")
}

if(is.na(opt$output)){
  filename = paste0(sce@metadata$project.name,"_SCE__PCA_FILTERED.rds")
  saveRDS(sce,file=filename, compress = T)
} else {
  saveRDS(sce,file=as.character(opt$output), compress = T)
}

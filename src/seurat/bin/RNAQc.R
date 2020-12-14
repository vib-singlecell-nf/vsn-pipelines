suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))

option_list = list(
  make_option(
    "--seuratObj",
    default=NA,
    type='character',
    help="File path to the rds file containing the seurat object."
  ),
  make_option(
    "--output",
    default=NA,
    type='character',
    help="output file name"
  ),
  make_option(
    "--mitoGenes",
    default=FALSE,
    type='logical',
    action="store_true",
    help="check for mitochondrial genes"
  ),
  make_option(
    "--covidGenes",
    default=FALSE,
    type='logical',
    action="store_true",
    help="check for covid genes"
  ),
  make_option(
    "--rbcGenes",
    default=FALSE,
    type='logical',
    action="store_true",
    help="check for red blood cell genes"
  ),
  make_option(
    "--genomeName1",
    default=NA,
    type='character',
    help="For multi genomes experiment, name of the human/mouse genome"
  ),
  make_option(
    "--genomeName2",
    default=NA,
    type='character',
    help="For multi genomes experiment, name of the other genome (if genome 1 is human, genome 2 can be mouse)"
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

seuratObj<- readRDS(file = opt$seuratObj)

if("subsample" %in% colnames(seuratObj@meta.data)){
  seuratObj<- subset(seuratObj, cells = colnames(seuratObj)[grepl("#",seuratObj$subsample)])
  batch <- seuratObj$subsample
}else{
  batch <- rep(seuratObj@project.name,ncol(seuratObj))
  seuratObj$subsample <- rep(seuratObj@project.name,ncol(seuratObj))
}

########################
########## PREP DATA
########################

#diagnostics
diagnostics<-list()
diagnostics[['dimRawDataRNA']]<-paste0(nrow(seuratObj@assays$RNA@counts)," genes - ",ncol(seuratObj@assays$RNA@counts)," cells")
diagnostics[['nrGenes']]<-nrow(seuratObj@assays$RNA@counts)
diagnostics[['nrCells']]<-ncol(seuratObj@assays$RNA@counts)
diagnostics[['UsedCells']]<-length(ncol(seuratObj@assays$RNA@counts))

nrZeros<-sum(seuratObj@assays$RNA@counts==0)/(nrow(seuratObj@assays$RNA@counts)*ncol(seuratObj@assays$RNA@counts))*100
##### In each cell: how many genes are expressed (count > 0) #####
cellCounts<-apply(seuratObj@assays$RNA@counts,2,function (x){sum(x>0)})
##### For each gene: in how many cells is it expressed? (count > 0) #####
geneCounts<-apply(seuratObj@assays$RNA@counts,1,function (x){sum(x>0)})
##### Add to diagnostics #####
diagnostics[['dimRawData']]<-paste0(nrow(seuratObj@assays$RNA@counts)," genes - ",ncol(seuratObj@assays$RNA@counts)," cells")
diagnostics[['nrGenes']]<-nrow(seuratObj@assays$RNA@counts)
diagnostics[['nrCells']]<-ncol(seuratObj@assays$RNA@counts)
diagnostics[['zeroInflation']]<-nrZeros
diagnostics[['minGenesPerCell']]<-min(cellCounts)
diagnostics[['maxGenesPerCell']]<-max(cellCounts)
diagnostics[['meanGenesPerCell']]<-mean(cellCounts)
diagnostics[['medianGenesPerCell']]<-median(cellCounts)
diagnostics[['cellsLess200genes']]<-length(cellCounts[cellCounts<200])
diagnostics[['genesNotExpressed']]<-length(geneCounts[geneCounts==0])
diagnostics[['genesLess3cells']]<-length(geneCounts[geneCounts<3])

###########################
########## QC: CELLS
###########################
##### Calculate QC metrics #####

##### Create metaData matrix (used for downstream analysis) #####
metaData<-data.frame("staticNr" = colnames(seuratObj@assays$RNA@counts),
                     "nGene" = cellCounts,
                     "nUMI" = Matrix::colSums(seuratObj),
                     "batch" = batch,
                     stringsAsFactors = F)
rownames(metaData)<-metaData$staticNr
metaData$staticNr <- 1


genome1Pattern <- ifelse(is.na(opt$genomeName1),"",paste0(opt$genomeName1,"(_|-){1,}"))
genome2Pattern <- ifelse(is.na(opt$genomeName2),"",paste0(opt$genomeName2,"(_|-){1,}"))

if(opt$mitoGenes){
  ##### Get mitochondrial genes #####
  is.mito <- grepl(paste0("(^|",genome2Pattern,"|",genome1Pattern,")MT-"), rownames(seuratObj), ignore.case = TRUE)
  metaData$percent.mito <- (Matrix::colSums(seuratObj@assays$RNA@counts[is.mito, ])/Matrix::colSums(seuratObj@assays$RNA@counts))*100
}

if(opt$covidGenes){
  #### Get COVID genes
  seuratObj@misc$covid.genes <- c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6",
                                  "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
  covidPatterns <- paste0("(^|",genome2Pattern,")",seuratObj@misc$covid.genes,"$")
  is.covidList <- lapply(covidPatterns,function(x) grepl(x, rownames(seuratObj), ignore.case = TRUE))
  is.covid <- as.logical(rowSums(simplify2array(is.covidList)))
  metaData$percent.COVID <- (Matrix::colSums(seuratObj@assays$RNA@counts[is.covid, ])/Matrix::colSums(seuratObj@assays$RNA@counts))*100
}

if(opt$rbcGenes){
  seuratObj@misc$rbc.genes <- c("ADIPOR1", "ALAS2", "ATP5E", "BAG1", "BCL2L1", "BNIP3L",
                                  "BPGM", "BTF3", "CA1", "DCAF12", "EPB42", "FBXO7", "FKBP8",
                                  "FTL", "GSPT1", "GUK1", "GYPC", "HBA1","HBB", "HBD", "HEMGN",
                                  "MPP1", "MYL6", "NCOA4", "OAZ1", "PFDN5", "RNF10", "RPL12",
                                  "RPL21", "RPL27A", "RPL30", "RPL31", "RPL32", "RPL38", "RPL41",
                                  "RPL7", "RPLP1", "RPLP2", "RPS11", "RPS12", "RPS13", "RPS14",
                                  "RPS24", "SELENBP1", "SERF2", "SLC25A37", "SLC25A39", "SNCA",
                                  "TPT1", "UBA52", "UBB", "YBX3")

  rbcPatterns <- paste0("(^|",genome1Pattern,")",seuratObj@misc$rbc.genes,"$")
  is.rbcList <- lapply(rbcPatterns,function(x) grepl(x, rownames(seuratObj), ignore.case = TRUE))
  is.rbc <- as.logical(rowSums(simplify2array(is.rbcList)))
  metaData$percent.rbc <- (Matrix::colSums(seuratObj@assays$RNA@counts[is.rbc, ])/Matrix::colSums(seuratObj@assays$RNA@counts))*100
}

print("########################## Plotting Data ##########################")
dir.create("Plots/RNA",recursive = T)
palette(c("#00BFC4","#F8766D","#7CAE00","#C77CFF"))
###F8766D=red
###00BFC4=cyan
###7CAE00=green
###C77CFF=purple

# calculating some limits for the plots

lognUMI <- log2(metaData$nUMI)
lognGene <- log2(metaData$nGene)
ll.cd <- 2^(median(lognUMI)-5*mad(lognUMI,na.rm = TRUE))
ul.cd <- 2^(median(lognUMI)+5*mad(lognUMI,na.rm = TRUE))
ll.ng <- 2^(median(lognGene)-5*mad(lognGene,na.rm = TRUE))
ul.ng <- 2^(median(lognGene)+5*mad(lognGene,na.rm = TRUE))

toPlot<-metaData

##nGene
png(file="Plots/RNA/01a_nGene.png", width=850, type = 'cairo')
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nGene),]
hist(tmp$nGene, breaks=30)
theColors<-as.factor(tmp$nGene < ll.ng | tmp$nGene > ul.ng )
barplot(tmp$nGene, col=theColors, border=theColors)
dev.off()

##nUMI
png(file="Plots/RNA/01a_nUMI.png", width=850, type = 'cairo')
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nUMI),]
hist(tmp$nUMI, breaks=30)
theColors<-as.factor(tmp$nUMI < ll.cd | tmp$nUMI > ul.cd)
barplot(tmp$nUMI, col=theColors, border=theColors)
dev.off()

png(file="Plots/RNA/01b_histogram_countdepth.png", width=850, type = 'cairo')
ggplot(metaData, aes(x = nUMI)) +
  geom_histogram(binwidth=100) +
  xlab("Count depth (total UMI count)") +
  ylab("Frequency") +
  geom_vline(xintercept=ul.cd,
             color = "red", size=1) +
  geom_vline(xintercept=ll.cd,
             color = "red", size=1) +
  facet_grid(.~batch) +
  theme_bw() # histogram of count depth
dev.off()

png(file="Plots/RNA/01b_histogram_genes.png", width=850, type = 'cairo')
ggplot(metaData, aes(x = nGene)) +
  geom_histogram(binwidth=20) +
  xlab("Number of Genes") +
  ylab("Frequency") +
  geom_vline(xintercept=ul.ng,
             color = "red", size=1) +
  geom_vline(xintercept=ll.ng,
             color = "red", size=1) +
  facet_grid(.~batch) +
  theme_bw() # histogram of nr of genes
dev.off()

if(opt$mitoGenes){
  ul.mito <- median(metaData$percent.mito)+5*mad(metaData$percent.mito,na.rm = TRUE)
  ##percent.mito

  png(file="Plots/RNA/01a_percMito.png", width=850, type = 'cairo')
  par(mfrow=c(1,2))
  tmp<-toPlot[order(toPlot$percent.mito),]
  hist(tmp$percent.mito, breaks=30)
  theColors<-as.factor(tmp$percent.mito > ul.mito)
  barplot(tmp$percent.mito, col=theColors, border=theColors)
  dev.off()

  png(file= "Plots/RNA/01b_histogram_percentageMito.png", width=850, type = 'cairo')
  g1b <- ggplot(metaData, aes(x = percent.mito)) +
    geom_histogram(binwidth=0.1) +
    xlab("% Mitochondrial counts") +
    ylab("Frequency") +
    geom_vline(xintercept=ul.mito,
               color = "red", size=1) +
    facet_grid(.~batch) +
    theme_bw()
  print(g1b)
  dev.off()

  png(file="Plots/RNA/01c_scatterplot_filtering.png", width=850, type = 'cairo')
  scatterplot1c <- ggplot(metaData, aes(x = nUMI,y=nGene,colour=percent.mito)) +
    geom_point(size=0.5) +
    scale_color_gradient2(midpoint=ul.mito, low="black", mid="white",
                          high="red", space ="Lab" )+
    xlab("Count depth (total UMI count)") +
    ylab("Nr of Genes") +
    geom_vline(xintercept=ul.cd,
               color = "red", size=1) +
    geom_vline(xintercept=ll.cd,
               color = "red", size=1) +
    geom_hline(yintercept=ul.ng,
               color = "red", size=1) +
    geom_hline(yintercept=ll.ng,
               color = "red", size=1) +
    geom_rug(col=rgb(0,0,0.5,alpha=.1)) +
    facet_grid(.~batch) +
    theme_bw()
  print(scatterplot1c)
  dev.off()

}

seuratObj@meta.data <- cbind(seuratObj@meta.data,metaData)
seuratObj@tools$diagnostics <- c(seuratObj@tools$diagnostics, diagnostics)

if(is.na(opt$output)){
  filename = paste0(seuratObj@project.name,".SEURAT__RNA_QC.rds")
  saveRDS(seuratObj,file=filename, compress = T)
} else {
  saveRDS(seuratObj,file=as.character(opt$output))
}

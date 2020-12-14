suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("ggplot2"))
option_list = list(
  make_option(
    c("-r", "--rdsFile"), 
    dest="rdsFile",
    default=NA, 
    type='character',
    help="Path to a .rds file containing a SCE object "),
  make_option(
    c("-o", "--output"),
    dest="output",
    type='character',
    help="Output file path to save the SingleCellExperiment object as an Rds file."),
  make_option(
    c("-m", "--nmads"),
    default=3,
    type='numeric',
    help="A numeric scalar, specifying the minimum number of MADs away from median required for a value to be called an outlier"),
  make_option(
    c("-i", "--mito"),
    default=FALSE,
    type='logical',
    action="store_true",
    help="Mitochondrial genes used to filter cells")
)
opt = parse_args(OptionParser(option_list=option_list))

sce <- readRDS(opt$rdsFile)

# Determine the outliers
if(is.null(sce$detected) || is.null(sce$sum)){
  if(opt$mito){
    is.mito <- grepl("^MT-", rownames(sce), ignore.case = TRUE)
    sce <- addPerCellQC(sce, subsets=list(Mito=is.mito))
  }else{
    sce <- addPerCellQC(sce)
  }
}

if(opt$mito){
  out.mito <- isOutlier(sce$subsets_Mito_percent, nmads=opt$nmads, type="higher")
}
out.ngene <- isOutlier(sce$detected, nmads=opt$nmads, type="both")
out.numi <- isOutlier(sce$sum, nmads=opt$nmads, type="both")

# Plots 

if(opt$mito){
  metaData<-data.frame("Experiment"=colData(sce)$Barcode,
                       "nGene"=sce$detected,
                       "nUMI"=sce$sum,
                       "percent.mito"=sce$subsets_Mito_percent,
                       "out.mito"=out.mito,
                       "out.ngene"=out.ngene,
                       "out.numi"=out.numi,
                       stringsAsFactors = F)
  ul.mito <- median(metaData$percent.mito)+3*mad(metaData$percent.mito,na.rm = TRUE)
}else{
  metaData<-data.frame("Experiment"=colData(sce)$Barcode,
                       "nGene"=sce$detected,
                       "nUMI"=sce$sum,
                       "out.ngene"=out.ngene,
                       "out.numi"=out.numi,
                       stringsAsFactors = F)
}

ll.cd <- median(metaData$nUMI)-3*mad(metaData$nUMI,na.rm = TRUE)
ul.cd <- median(metaData$nUMI)+3*mad(metaData$nUMI,na.rm = TRUE)
ll.ng <- median(metaData$nGene)-3*mad(metaData$nGene,na.rm = TRUE)
ul.ng <- median(metaData$nGene)+3*mad(metaData$nGene,na.rm = TRUE)

if(opt$mito){
  sce@metadata$plots$outliersDetection1 <- ggplot(metaData, aes(x = nUMI,y=nGene,colour=percent.mito)) + 
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
    theme_bw()
  
  sce@metadata$plots$outliersDetection2 <- ggplot(metaData, aes(x = nUMI,y=nGene,colour=out.mito)) + 
    geom_point(size=0.5) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
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
    theme_bw()
}else{
  sce@metadata$plots$outliersDetection1 <- ggplot(metaData, aes(x = nUMI,y=nGene)) + 
    geom_point(size=0.5) +
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
    theme_bw()
  
  sce@metadata$plots$outliersDetection2 <- ggplot(metaData, aes(x = nUMI,y=nGene)) + 
    geom_point(size=0.5) +
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
    theme_bw()
}

# Gather metadata for qc report
if(opt$mito){
  sce@metadata$outliersDetection$out.mito <- out.mito
}
sce@metadata$outliersDetection$out.ngene <- out.ngene
sce@metadata$outliersDetection$out.numi <- out.numi

# Exclude the outliers
if(opt$mito){
  keep <- !(out.ngene | out.numi | out.mito)
}else{
  keep <- !(out.numi | out.ngene)
}
sce$PassQC <- keep
sce <- sce[,keep]

saveRDS(sce, opt$output, compress = T)
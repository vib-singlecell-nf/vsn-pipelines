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
    c("-i", "--mito"),
    default=FALSE,
    type='logical',
    action="store_true",
    help="Check the mitochoindrial gene")
)
opt = parse_args(OptionParser(option_list=option_list))

sce <- readRDS(opt$rdsFile)

if(opt$mito){
  is.mito <- grepl("(^|_|-)MT-", rownames(sce), ignore.case = TRUE)
  sce <- addPerCellQC(sce,subsets=list(Mito=is.mito))
}else{
  sce <- addPerCellQC(sce)
}

# QC plots
if(opt$mito){
  sce@metadata$scater_qcMetrics$percent.mito <- sce$pct_counts_Mt
  metaData<-data.frame("Experiment"=colData(sce)$Barcode,
                       "nGene"=sce$detected,
                       "nUMI"=sce$sum,
                       "percent.mito"=sce$subsets_Mito_percent,
                       stringsAsFactors = F)
}else{
  metaData<-data.frame("Experiment"=colData(sce)$Barcode,
                       "nGene"=sce$detected,
                       "nUMI"=sce$sum,
                       stringsAsFactors = F)
}
metaData$Experiment <- 1 #needed for violin plots
sce@metadata$plots$qcMetrics1 <- ggplot(metaData, aes(x = Experiment, y = nUMI)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3,size=0.5) +
  ylab("Count depth (total UMI count)") +
  theme_bw() # violinplot of count depth

sce@metadata$plots$qcMetrics2 <- ggplot(metaData, aes(x = nUMI)) + 
  geom_histogram(binwidth=100) +
  xlab("Count depth (total UMI count)") +
  ylab("Frequency") +
  theme_bw() # histogram of count depth

sce@metadata$plots$qcMetrics3 <- ggplot(metaData, aes(x = nGene)) + 
  geom_histogram(binwidth=20) +
  xlab("Number of Genes") +
  ylab("Frequency") +
  theme_bw() # histogram of nr of genes

sce@metadata$plots$qcMetrics4 <- ggplot(metaData, aes(x = Experiment, y = nGene)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3,size=0.5) +
  ylab("Nr of Genes") +
  theme_bw() # violinplot of nr of genes

if(opt$mito){
  sce@metadata$plots$qcMetrics5 <- ggplot(metaData, aes(x = percent.mito)) + 
    geom_histogram(binwidth=0.1) +
    xlab("% Mitochondrial counts") +
    ylab("Frequency") +
    theme_bw() # histogram of % mito
  
  sce@metadata$plots$qcMetrics6 <- ggplot(metaData, aes(x = Experiment, y = percent.mito)) + 
    geom_violin(fill="gray80") + 
    geom_jitter(height = 0, width = 0.3,size=0.5) +
    ylab("% Mitochondrial counts") +
    theme_bw() # violinplot of % mito
}

saveRDS(sce, opt$output, compress = T)
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("DropletUtils"))
suppressPackageStartupMessages(library("ggplot2"))
option_list = list(
  make_option(
    c("-r", "--rdsFile"),
    dest="rds_file_path",
    default=NA,
    type='character',
    help="File path to the Rds file containing a SingleCellExperiment object."
  ),
  make_option(
    c("-o", "--output"),
    dest="output",
    type='character',
    help="Output file path to save the SingleCellExperiment object as an Rds file."),
  make_option(
    c("-f", "--fdr-threshold"),
    dest="fdr_threshold",
    default=0.01,
    type='numeric',
    help="False discovery rate to filter the empty droplets."
  ),
  make_option(
    c("-l", "--lower"),
    dest="lower",
    default=100,
    type='numeric',
    help="A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets."
  )
)
opt <- parse_args(OptionParser(option_list=option_list))

cleanEmptyDrops <- function(sce.object, lower.umi, fdr.threshold) {
  e_out <- DropletUtils::emptyDrops(m=counts(sce.object),lower=as.numeric(lower.umi))
  e_out$is_cell <- e_out$FDR <= as.numeric(fdr.threshold)
  sce.object <- sce.object[,which(e_out$is_cell == T)]
  counts_gene <- Matrix::rowSums(counts(sce.object), na.rm = FALSE, dims = 1)
  sce.object <- sce.object[which(counts_gene != 0),]
  # add metadata for qc plot
  sce.object@metadata$plots$emptydrops <- ggplot(na.omit(as.data.frame(e_out)),aes(x=Total,y=-LogProb,col=is_cell)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    xlab("Total UMI count") +
    ylab("-Log Probability") +
    ggtitle("Empty Droplets detection") +
    theme_bw()
  return (sce.object)
}

sce <- readRDS(opt$rds_file_path)
sce <- cleanEmptyDrops(sce.object=sce, lower.umi=opt$lower,fdr.threshold=opt$fdr_threshold)

saveRDS(sce, opt$output, compress=TRUE)
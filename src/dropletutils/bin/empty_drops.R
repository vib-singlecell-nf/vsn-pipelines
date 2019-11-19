#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("DropletUtils"))

option_list = list(
  make_option(
    c("-r", "--rds-file-path"),
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

args = parse_args(OptionParser(option_list=option_list))

cleanEmptyDrops <- function(sce, lower, fdr.threshold) {
  e_out <- DropletUtils::emptyDrops(
    m=counts(sce),
    lower=lower
  )
  is_cell <- e_out$FDR <= fdr.threshold
  sce <- sce[,which(is_cell == T)]
  counts_gene <- Matrix::rowSums(counts(sce), na.rm = FALSE, dims = 1)
  sce <- sce[which(counts_gene != 0),]
  return (sce)
}

sce <- readRDS(args$rds_file_path)
bcrank <- barcodeRanks(m=counts(sce))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
bcrank.plot <- as.data.frame(bcrank[uniq,]@listData)
sce@metadata$knee <- bcrank@metadata$knee
sce@metadata$inflection <- bcrank@metadata$inflection
sce@metadata$rank <- bcrank.plot$rank
sce@metadata$total <- bcrank.plot$total

sce <- cleanEmptyDrops(
  sce=sce, 
  lower=args$lower,
  fdr.threshold=args$fdr_threshold
)
saveRDS(sce, args$output, compress=TRUE)

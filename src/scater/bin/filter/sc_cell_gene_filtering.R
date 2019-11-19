#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("scater"))

option_list = list(
  make_option(c("-r", "--rdsFile"), default=NA, type='character',
              help="Path to a .rds file containing a SCE object "),
  make_option(c("-o", "--output"), default=NA, type='character',
              help="Name of output file"),
  make_option(c("-m", "--nmads"), default=3, type='numeric',
              help="A numeric scalar, specifying the minimum number of MADs away from median required for a value to be called an outlier")
)
opt = parse_args(OptionParser(option_list=option_list))

sce <- readRDS(opt$rdsFile)

# TODO: Make work for multiple species
is.mito <- grepl("^MT-", rownames(sce), ignore.case = TRUE)
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito))

out.mito <- isOutlier(sce$pct_counts_Mt, nmads=opt$nmads, type="higher")
out.ngene <- isOutlier(sce$total_features_by_counts, nmads=opt$nmads, type="both")
out.numi <- isOutlier(sce$total_counts, nmads=opt$nmads, type="both")

# Test
#out.mito <- isOutlier(sce$pct_counts_Mt, nmads=3, type="higher")
#out.ngene <- isOutlier(sce$total_features_by_counts, nmads=3, type="both")
#out.numi <- isOutlier(sce$total_counts, nmads=3, type="both")

sce@metadata$out.mito <- out.mito
sce@metadata$out.ngene <- out.ngene
sce@metadata$out.numi <- out.numi
sce@metadata$Experiment <- colData(sce)$Barcode
sce@metadata$nGene <- sce$total_features_by_counts
sce@metadata$nUMI <- sce$total_counts
sce@metadata$percent.mito <- sce$pct_counts_Mt

keep <- !(out.ngene | out.numi | out.mito)
sce$PassQC <- keep
sce <- sce[,keep]

saveRDS(sce, opt$output, compress=TRUE)
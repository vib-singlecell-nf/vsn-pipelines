#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("scater"))

option_list = list(
  make_option(c("-i", "--input-format"), default=NA, type='character',
              help="Input file format"),
  make_option(c("-o", "--output-format"), default=NA, type='character',
              help="Output file format")
)
arguments = parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)
opt = arguments$options
files = arguments$args

if (opt$i != '10x_mtx') {
    stop("Input file format not yet implemented")
}

if (opt$o != 'sce.rds') {
    stop("Output file format not yet implemented")
}

if(!is.na(files[1])){
  sce <- DropletUtils::read10xCounts(files[1],col.names = T)
} else {
  stop("Please provide a h5 file path or a path to a mtx file folder")
}

if(!is.null(rowData(sce)$Type)){
  is.RNA <- rowData(sce)$Type=="Gene Expression" 
  sce <- sce[is.RNA,] # only Gene Expression data
}
# getting proper names for genes

saveRDS(sce, files[2], compress=TRUE)
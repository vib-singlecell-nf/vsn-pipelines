#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

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

if (opt$o == 'sce.rds') {
  suppressPackageStartupMessages(library("scater"))

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

} else if (opt$o == 'seurat.rds') {
  suppressPackageStartupMessages(library("Seurat"))

  if (file.exists(file.path(files[1], 'genes.tsv'))) {
    data <- Read10X(data.dir = files[1])
    seurat_object = CreateSeuratObject(counts = data)

  } else if (file.exists(file.path(files[1], 'features.tsv.gz'))) {
    data <- Read10X(data.dir = files[1])
    seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)

  } else {
    
    stop("Could not find either 'genes.tsv' or 'features.tsv.gx' to determine 10X output version")
  }
  saveRDS(seurat_object, files[2], compress=TRUE)

} else {
  stop("Output file format not yet implemented")
}
#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option("--inputMatrix", 
              default=NA, 
              type='character',
              help="Path to a .h5 count matrix or a folder containing the files : matrix.mtx(.gz), features.tsv(.gz or genes.tsv) and barcodes.tsv(.gz)"),
  make_option("--output", 
              default=NA, 
              type='character',
              help="Output name"),
  make_option("--sampleName", 
              default=NA, 
              type='character',
              help="sample name")
)

opt = parse_args(OptionParser(option_list=option_list))

readSTARsoloCounts <- function(sample) {
  dir.data <- normalizePath(sample, mustWork = TRUE)
  m <- methods::as(Matrix::readMM(paste0(dir.data, "/","matrix.mtx")), "dgTMatrix")
  gene_info <- data.table::fread(paste0(dir.data, "/", "features.tsv"), header = FALSE,col.names=c("ID","Symbol"),stringsAsFactors = F)  
  cell_info <- data.table::fread(paste0(dir.data, "/", "barcodes.tsv"), header = FALSE,stringsAsFactors = F,col.names = c("Barcode"))
  rownames(gene_info) <- gene_info$ID
  colnames(m) <- cell_info$Barcode
  sce.obj <- SingleCellExperiment(list(counts = m), rowData = gene_info, colData = cell_info,metadata=list(Sample_path=sample))
  return(sce.obj)
}

if(grepl("\\.h5",as.character(opt$inputMatrix))){
  sce.counts <- read10xCounts(as.character(opt$inputMatrix),col.names = T)
} else{
  file.names <- list.files(as.character(opt$inputMatrix))
  if("features.tsv.gz" %in% file.names){
    sce.counts <- read10xCounts(as.character(opt$inputMatrix),col.names = T)
  } else if ("genes.tsv" %in% file.names){
    sce.counts <- read10xCounts(as.character(opt$inputMatrix),col.names = T)
  } else if ("features.tsv" %in% file.names){
    sce.counts <- readSTARsoloCounts(as.character(opt$inputMatrix))
  } else {
    stop("unrecognized input format")
  }
}

if(!is.null(rowData(sce.counts)$Type)){
  is.RNA <- rowData(sce.counts)$Type=="Gene Expression" 
  sce.counts <- sce.counts[is.RNA,] # only Gene Expression data
}
# getting proper names for genes
rownames(sce.counts) <- uniquifyFeatureNames(rowData(sce.counts)$ID, rowData(sce.counts)$Symbol)
sce.counts@metadata$project.name <- opt$sampleName
if(is.na(opt$output)){
  filename = "SINGLE_CELL_EXPERIMENT__SCE_BUILDER_output.rds"
  saveRDS(sce.counts,file=filename, compress = T)
} else {
  saveRDS(sce.counts,file=as.character(opt$output), compress = T)
}

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
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
    help="Output csv file path to save the Cell Barcode names")
)
opt <- parse_args(OptionParser(option_list=option_list))

sce <- readRDS(opt$rds_file_path)

cbarcodes <- as.data.frame(colnames(sce))

write.table(cbarcodes, opt$output,col.names = F,row.names = F)

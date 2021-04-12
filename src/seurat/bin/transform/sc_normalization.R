#!/usr/bin/env Rscript

library("argparse")
suppressMessages(library("Seurat", quietly = TRUE))

parser <- ArgumentParser(description = "Normalize and scale data using SCT")

parser$add_argument(
    "--input",
    type = "character",
    dest = "input",
    action = "store",
    help = "A Rds file containing a Seurat object."
)
parser$add_argument(
    "--output",
    type = "character",
    dest = "output",
    action = "store",
    help = "Output filename."
)
parser$add_argument(
    "--method",
    type = "character",
    dest = "method",
    action = "store",
    default = "LogNormalize",
    help = "Method for normalization. Choose one of: LogNormalize, CLR, RC"
)
parser$add_argument(
    "--scale-factor",
    type = "integer",
    dest = "scale_factor",
    action = "store",
    default = 10000,
    help = "Sets the scale factor for cell-level normalization."
)

args <- parser$parse_args()
print(args)

seuratObj <- tryCatch({
    readRDS(file = args$input)
}, warning = function(w) {
}, error = function(e) {
    stop("VSN ERROR: Cannot read the given Rds file.")
}, finally = {})

if (class(x = seuratObj) != "Seurat") {
    stop("VSN ERROR: The object contained in the Rds file is not a Seurat object.")
}

seuratObj <- NormalizeData(
    object = seuratObj,
    normalization.method = args$method,
    scale.factor = args$scale_factor
)

saveRDS(
    object = seuratObj,
    file = args$output
)

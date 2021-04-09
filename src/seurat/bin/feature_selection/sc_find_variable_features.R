#!/usr/bin/env Rscript

library("argparse")
suppressMessages(library("Seurat", quietly = TRUE))

parser <- ArgumentParser(description = "Find high variable features")

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
    help = "Method to choose top variable features. Choose on off: vst, mvp (mean.var.plot), disp (dispersion)"
)
parser$add_argument(
    "--n-features",
    type = "integer",
    dest = "n_features",
    action = "store",
    default = 2000,
    help = "Number of features to select as top variable features; only used when method is one of: vst, disp"
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

seuratObj <- FindVariableFeatures(
    object = seuratObj,
    method = args$method,
    nfeatures = args$n_features,
)

saveRDS(
    object = seuratObj,
    file = args$output
)

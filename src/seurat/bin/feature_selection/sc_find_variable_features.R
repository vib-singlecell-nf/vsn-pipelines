#!/usr/bin/env Rscript

library("optparse")
suppressMessages(library("Seurat", quietly = TRUE))

option_list <- list(
    make_option(
        "--input",
        type = "character",
        dest = "input",
        help = "A Rds file containing a Seurat object."
    ),
    make_option(
        "--output",
        type = "character",
        dest = "output",
        help = "Output filename."
    ),
    make_option(
        "--method",
        type = "character",
        dest = "method",
        help = "Method to choose top variable features. Choose on off: vst, mvp (mean.var.plot), disp (dispersion)"
    ),
    make_option(
        "--n-features",
        type = "integer",
        dest = "n_features",
        default = 2000,
        help = "Number of features to select as top variable features; only used when method is one of: vst, disp"
    )
)

args <- parse_args(OptionParser(option_list = option_list))
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
    nfeatures = args$n_features
)

saveRDS(
    object = seuratObj,
    file = args$output
)

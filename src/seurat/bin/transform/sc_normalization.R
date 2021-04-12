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
        default = "LogNormalize",
        help = "Method for normalization. Choose one of: LogNormalize, CLR, RC"
    ),
    make_option(
        "--scale-factor",
        type = "integer",
        dest = "scale_factor",
        default = 10000,
        help = "Sets the scale factor for cell-level normalization."
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

seuratObj <- NormalizeData(
    object = seuratObj,
    normalization.method = args$method,
    scale.factor = args$scale_factor
)

saveRDS(
    object = seuratObj,
    file = args$output
)

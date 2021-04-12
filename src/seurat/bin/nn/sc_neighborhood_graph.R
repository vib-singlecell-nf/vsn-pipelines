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
        "--n-pcs",
        type = "integer",
        dest = "n_pcs",
        help = "Use this many PCs."
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

seuratObj <- FindNeighbors(
    object = seuratObj,
    dims = 1:args$n_pcs
)

saveRDS(
    object = seuratObj,
    file = args$output
)

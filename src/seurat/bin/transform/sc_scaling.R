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
        "--only-var",
        type = "logical",
        dest = "only_var",
        default = TRUE,
        help = "Only scale variable features."
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

IsTrue <- function(x) {
    x <- tolower(x = as.character(x))
    return (
        x == 't' | x == 'true' | x == 'y' | x == 'yes'
    )
}

if (IsTrue(args$only_var)) {
    features <- rownames(seuratObj)
} else {
    features <- Seurat::VariableFeatures(seuratObj)
}

seuratObj <- ScaleData(
    object = seuratObj,
    features = features
)

saveRDS(
    object = seuratObj,
    file = args$output
)

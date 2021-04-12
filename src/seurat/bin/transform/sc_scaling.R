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
    "--only-var",
    type = "logical",
    dest = "only_var",
    action = "store",
    default = TRUE,
    help = "Only scale variable features."
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

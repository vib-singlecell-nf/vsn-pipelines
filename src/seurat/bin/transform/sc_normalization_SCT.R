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
    "--new-assay",
    type = "character",
    dest = "new_assay",
    action = "store",
    default = "SCT",
    help = "Name for the assay created by SCT."
)
parser$add_argument(
    "--n-cells",
    type = "integer",
    dest = "n_cells",
    action = "store",
    default = 5000,
    help = "Number of subsampling cells used to build NB regression."
)
parser$add_argument(
    "--only-var",
    type = "logical",
    dest = "only_var",
    action = "store",
    default = TRUE,
    help = "Only output variable features in scaled data. This can drastically reduce object size."
)
parser$add_argument(
    "--seed",
    type = "integer",
    dest = "seed",
    action = "store",
    default = 0,
    help = "Use this integer seed for reproducibility."
)
parser$add_argument(
    "--n-variable-features",
    type = "character",
    dest = "n_variable_features",
    action = "store",
    default = 3000,
    help = "Use this many features as variable features after ranking by residual variance"
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

seuratObj <- SCTransform(
    object = seuratObj,
    new.assay.name = args$new_assay,
    ncells = args$n_cells,
    return.only.var.genes = IsTrue(args$only_var),
    variable.features.n = args$n_variable_features,
    seed.use = args$seed,
)

saveRDS(
    object = seuratObj,
    file = args$output
)

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
        "--new-assay",
        type = "character",
        dest = "new_assay",
        default = "SCT",
        help = "Name for the assay created by SCT."
    ),
    make_option(
        "--n-cells",
        type = "integer",
        dest = "n_cells",
        default = 5000,
        help = "Number of subsampling cells used to build NB regression."
    ),
    make_option(
        "--only-var",
        type = "logical",
        dest = "only_var",
        default = TRUE,
        help = "Only output variable features in scaled data. This can drastically reduce object size."
    ),
    make_option(
        "--seed",
        type = "integer",
        dest = "seed",
        default = 0,
        help = "Use this integer seed for reproducibility."
    ),
    make_option(
        "--n-variable-features",
        type = "character",
        dest = "n_variable_features",
        default = 3000,
        help = "Use this many features as variable features after ranking by residual variance"
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

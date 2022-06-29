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
        "--algorithm",
        type = "character",
        dest = "algorithm",
        help = "Cluster cells. Choose one of: Louvain, Louvain_MLR, SLM, Leiden"
    ),
    make_option(
        "--resolution",
        type = "double",
        dest = "resolution",
        default = 0.8
    ),
    make_option(
        "--seed",
        type = "integer",
        dest = "seed",
        default = 0,
        help = "Use this numeric seed for reproducibility."
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

algStringToNumber <- function (algorithm) {
    algorithm <- tolower(x = algorithm)

    if (algorithm == 'louvain') {
        return(1)
    } else if (algorithm == 'louvain_mlr') {
        return(2)
    } else if (algorithm == 'slm') {
        return(3)
    } else if (algorithm == 'leiden') {
        return(4)
    } else {
        stop(paste0("VSN ERROR: Invalid clustering algorithm ", algorithm))
    }
}

seuratObj <- Seurat::FindClusters(
    object = seuratObj,
    resolution = args$resolution,
    algorith = algStringToNumber(args$algorithm),
    random.seed = args$seed
)

saveRDS(
    object = seuratObj,
    file = args$output
)

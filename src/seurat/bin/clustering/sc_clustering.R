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
    "--algorithm",
    type = "character",
    dest = "algorithm",
    action = "store",
    help = "Cluster cells. Choose one of Louvain, Louvain_MLR, SLM, Leiden"
)
parser$add_argument(
    "--resolution",
    type = "double",
    dest = "resolution",
    action = "store",
    default = 0.8
)
parser$add_argument(
    "--seed",
    type = "integer",
    dest = "seed",
    action = "store",
    default = 0,
    help = "Use this integer seed for reproducibility."
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

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
    help = "Reduce the dimensionality of the data. Choose one of : PCA, UMAP, tSNE"
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
    "--n-comps",
    type = "integer",
    dest = "n_comps",
    action = "store",
    default = 50,
    help = "[PCA] Number of principal components to compute."
)
parser$add_argument(
    "--n-pcs",
    type = "integer",
    dest = "n_pcs",
    action = "store",
    default = 30,
    help = "[UMAP, tSNE] Use this many PCs."
)
parser$add_argument(
    "--n-neighbors",
    type = "integer",
    dest = "n_neighbors",
    action = "store",
    default = 30,
    help = "[UMAP] Use this many neighbors."
)
parser$add_argument(
    "--algorithm",
    type = "character",
    dest = "algorithm",
    action = "store",
    default = NULL,
    help = "[UMAP, tSNE] Algorithm to use"
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

if (tolower(args$method) == "pca") {
    seuratObj <- Seurat::RunPCA(
        object = seuratObj,
        npcs = args$n_comps,
        verbose = FALSE
    )
} else if (tolower(args$method == "umap")) {

    algorithm <- if (is.null(args$algorithm)) "uwot" else args$algorithm

    seuratObj <- Seurat::RunUMAP(
        object = seuratObj,
        umap.method = algorithm,
        dims = 1:args$n_pcs,
        n.neighbors = args$n_neighbors,
        seed.use = args$seed,
        verbose = FALSE
    )
} else if (tolower(args$method == "tsne")) {
    algorithm <- if (is.null(args$algorithm)) "Rtsne" else args$algorithm

    seuratObj <- Seurat::RunTSNE(
        object = seuratObj,
        tsne.method = algorithm,
        dims = 1:args$n_pcs,
        seed.use = args$seed,
        verbose = FALSE
    )
} else {
    stop(paste0("VSN ERROR: The dimensionality reduction method ", args$method, " does not exits."))
}

saveRDS(
    object = seuratObj,
    file = args$output
)

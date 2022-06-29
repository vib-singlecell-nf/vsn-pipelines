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
        help = "Reduce the dimensionality of the data. Choose one of : PCA, UMAP, tSNE"
    ),
    make_option(
        "--seed",
        type = "integer",
        dest = "seed",
        default = 0,
        help = "Use this integer seed for reproducibility."
    ),
    make_option(
        "--n-comps",
        type = "integer",
        dest = "n_comps",
        default = 50,
        help = "[PCA] Number of principal components to compute."
    ),
    make_option(
        "--n-pcs",
        type = "integer",
        dest = "n_pcs",
        default = 30,
        help = "[UMAP, tSNE] Use this many PCs."
    ),
    make_option(
        "--n-neighbors",
        type = "integer",
        dest = "n_neighbors",
        default = 30,
        help = "[UMAP] Use this many neighbors."
    ),
    make_option(
        "--algorithm",
        type = "character",
        dest = "algorithm",
        default = NULL,
        help = "[UMAP, tSNE] Algorithm to use"
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

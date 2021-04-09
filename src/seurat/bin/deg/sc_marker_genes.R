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
    help = "Which test to use. Choose one of : wilcox, bimod, roc, t, negbinom, poisson, LR, MAST, DESeq2, "
)
parser$add_argument(
    "--logfc-threshold",
    type = "double",
    dest = "logfc_threshold",
    action = "store",
    default = 0.25,
    help = "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells."
)
parser$add_argument(
    "--min-pct",
    type = "double",
    dest = "min_pct",
    action = "store",
    default = 0.1,
    help = "Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations."
)
parser$add_argument(
    "--only-pos",
    type = "logical",
    dest = "only_pos",
    action = "store",
    default = FALSE,
    help = "Only return positive markers."
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

IsTrue <- function(x) {
    x <- tolower(x = as.character(x))
    return (
        x == 't' | x == 'true' | x == 'y' | x == 'yes'
    )
}

markers <- Seurat::FindAllMarkers(
    object = seuratObj,
    test.use = args$method,
    logfc.threshold = args$logfc_threshold,
    min.pct = args$min_pct,
    only.pos = IsTrue(args$only_pos)
)

saveRDS(
    object = markers,
    file = args$output
)

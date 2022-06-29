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
        help = "Which test to use. Choose one of : wilcox, bimod, roc, t, negbinom, poisson, LR, MAST, DESeq2, "
    ),
    make_option(
        "--logfc-threshold",
        type = "double",
        dest = "logfc_threshold",
        default = 0.25,
        help = "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells."
    ),
    make_option(
        "--min-pct",
        type = "double",
        dest = "min_pct",
        default = 0.1,
        help = "Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations."
    ),
    make_option(
        "--only-pos",
        type = "logical",
        dest = "only_pos",
        default = FALSE,
        help = "Only return positive markers."
    ),
    make_option(
        "--seed",
        type = "integer",
        dest = "seed",
        default = 0,
        help = "Use this integer seed for reproducibility."
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

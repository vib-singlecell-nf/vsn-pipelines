#!/usr/bin/env Rscript

library("argparse")
suppressMessages(library("Seurat", quietly = TRUE))
suppressMessages(library("SCopeLoomR", quietly = TRUE))

parser <- ArgumentParser(description = "Convert Seurat RDS file to SCope compatible loom file")

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
    "--sampleId",
    type = "character",
    dest = "sample_id",
    action = "store",
    help = "Sample name."
)
parser$add_argument(
    "--assay",
    type = "character",
    dest = "assay",
    action = "store",
    help = "Assay to use"
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

if (! args$assay %in% names(seuratObj@assays)) {
    stop(paste0("VSN ERROR: Assay ", args$assay, " not available in Seurat object"))
}

default.reduction <- Seurat:::DefaultDimReduc(seuratObj)

# Create loom file with GEX and default embedding
build_loom(
    file.name = args$output,
    dgem = seuratObj@assays$[[args$assay]]@counts,
    title = args$sample_id,
    default.embedding = seuratObj@reductions[[default.reduction]]@cell.embeddings,
    default.embedding.name = default.reduction
)

loom <- open_lom(args$output, mode = "r+")

# Add all remaining embeddings
foreach (names(seuratObj@reductions) as reduction) {
    if (reduction == default.reduction) next

    add_embedding(
        loom = loom,
        embedding = seuratObj@reductions[[reduction]]@cell.embeddings,
        name = reduction
    )
}

add_seurat_clustering(
    loom = loom,
    seurat = seuratObj,
    seurat.assay = args$assay
)

close_loom(loom)

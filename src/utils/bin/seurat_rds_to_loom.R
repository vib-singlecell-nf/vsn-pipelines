#!/usr/bin/env Rscript

library("optparse")
suppressMessages(library("Seurat", quietly = TRUE))
suppressMessages(library("SCopeLoomR", quietly = TRUE))

option_list <- list(
    make_option(
        "--input",
        type = "character",
        dest = "input",
        action = "store",
        help = "A Rds file containing a Seurat object."
    ),
    make_option(
        "--output",
        type = "character",
        dest = "output",
        action = "store",
        help = "Output filename."
    ),
    make_option(
        "--sampleId",
        type = "character",
        dest = "sample_id",
        action = "store",
        help = "Sample name."
    ),
    make_option(
        "--assay",
        type = "character",
        dest = "assay",
        action = "store",
        default = "RNA",
        help = "Assay to use"
    ),
    make_option(
        "--clustering-prefix",
        dest = "clustering_prefix",
        action = "store",
        default = "RNA_snn_res.",
        help = "Prefix for metadata columns with clustering results."
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

if (! args$assay %in% names(seuratObj@assays)) {
    stop(paste0("VSN ERROR: Assay ", args$assay, " not available in Seurat object"))
}

default.reduction <- Seurat:::DefaultDimReduc(seuratObj)

# Create loom file with GEX and default embedding
build_loom(
    file.name = args$output,
    dgem = seuratObj@assays[[args$assay]]@counts,
    title = args$sample_id,
    default.embedding = seuratObj@reductions[[default.reduction]]@cell.embeddings,
    default.embedding.name = default.reduction
)

loom <- open_loom(args$output, mode = "r+")

# Add all remaining embeddings
for (reduction in names(seuratObj@reductions)) {
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
    seurat.assay = args$assay,
    seurat.clustering.prefix = args$clustering_prefix
)

close_loom(loom)

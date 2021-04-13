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
        "--type",
        type = "character",
        dest = "type",
        help = "What to filter, Choose one of: cell, feature"
    ),
    make_option(
        "--min-number-cells",
        type = "integer",
        dest = "min_number_cells",
        help = "[FEATURE] Minimal number of cells the feature should be expressed in (counts > 0)."
    ),
    make_option(
        "--min-n-counts",
        type = "integer",
        dest = "min_n_counts",
        default = -1,
        help = "[CELL] Minimal number of counts for a cell to be kept."
    ),
    make_option(
        "--max-n-counts",
        type = "integer",
        dest = "max_n_counts",
        default = -1,
        help = "[CELL] Maximal number of counts for a cell to be kept."
    ),
    make_option(
        "--min-n-features",
        type = "integer",
        dest = "min_n_features",
        default = -1,
        help = "[CELL] Minimal number of features for a cell to be kept."
    ),
    make_option(
        "--max-n-features",
        type = "integer",
        dest = "max_n_features",
        default = -1,
        help = "[CELL] Maximal number of features for a cell to be kept."
    ),
    make_option(
        "--max-percent-mito",
        type = "double",
        dest = "max_percent_mito",
        default = -1,
        help = "[CELL] Maximal percent mitochondrial genes a cell can have."
    ),
    make_option(
        "--mito-prefix",
        type = "character",
        dest = "mito_prefix",
        default = "MT-",
        help = "[CELL] Prefix to deterimine of a feature is mitochondrial or not (case sensitive!)"
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

filterCells <- function(
    object,
    min_n_counts,
    max_n_counts,
    min_n_features,
    max_n_features,
    max_percent_mito,
    mito_prefix
) {
    if (min_n_counts < 0) {
        min_n_counts <- -Inf
    }
    if (max_n_counts < 0) {
        max_n_counts <- Inf
    }
    if (min_n_features < 0) {
        min_n_features <- -Inf
    }
    if (max_n_features < 0) {
        max_n_features <- Inf
    }
    if (max_percent_mito < 0) {
        max_percent_mito <- Inf
    } else if (max_percent_mito < 1) {
        max_percent_mito <- max_percent_mito * 100
    }

    mito.pattern <- paste0("^", mito_prefix)

    mito.genes <- grep(mito.pattern, rownames(object))
    if (length(mito.genes) == 0 && max_percent_mito != Inf) {
        stop(paste0("VSN ERROR: Could not find any mitochondrial genes with the prefix: '" , mito_prefix, "'"))
    }

    object[["pct.mito"]] <- Seurat::PercentageFeatureSet(object, pattern = mito.pattern)
    object <- subset(
        object,
        subset = nFeature_RNA > min_n_features & nFeature_RNA < max_n_features & nCount_RNA > min_n_counts & nCount_RNA < max_n_counts & pct.mito < max_percent_mito
    )

    return(object)
}

# FIXME: this is a very dirty approach for feature filtering
# AFAIK it is not possible to do this in a existing object, since the features might be referenced internally (not enough to remove some rows from the matrices)
# The current approach assumes there is no information in the object of clustering/SNN/dimred... other than some QC metadata.
# If there woulde be other data, this will simply be lost since we only copy the metadata over to the new filtered object
filterFeatures <- function(
    object,
    min_number_cells
) {
    if (min_number_cells < 0) {
        return(object)
    }
    
    new_object <- Seurat::CreateSeuratObject(
        object@assays$RNA@data,
        min.cells = min_number_cells
    )
    new_object@meta.data <- object@meta.data

    return(new_object)
}

if (args$type == 'cell') {
    seuratObj <- filterCells(
        seuratObj,
        args$min_n_counts,
        args$max_n_counts,
        args$min_n_features,
        args$max_n_features,
        args$max_percent_mito,
        args$mito_prefix
    )
} else if (args$type == 'feature') {
    seuratObj <- filterFeatures(
        seuratObj,
        args$min_number_cells
    )
} else {
    stop(paste0("VSN ERROR: invalid filter type ", args$type))
}

saveRDS(
    object = seuratObj,
    file = args$output
)

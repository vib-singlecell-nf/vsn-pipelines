#!/usr/bin/env Rscript

in_formats = c("seurat_rds", "10x_cellranger_mex")
out_formats = c("h5ad", "sce_rds", "seurat_rds")

library("argparse")

parser <- ArgumentParser(description='Conversion script between single-cell data formats.')

parser$add_argument(
    'input',
    metavar='INPUT',
    action="store",
    type="character",
    help='A Rds file containing a SingleCellExperiment object.'
)
parser$add_argument(
    'output',
    metavar='OUTPUT',
    action="store",
    type="character",
    default=NULL,
    help='Output filename.'
)
parser$add_argument(
    '--input-format',
    type="character",
    dest='input_format',
    action = "store",
    default = NULL,
    help = paste0("Input format of the file to be converted. Choose one of: ", paste(in_formats, collapse = ", "), ".")
)
parser$add_argument(
    '--output-format',
    type="character",
    dest='output_format',
    action = "store",
    default = NULL,
    help = paste0("Output format which the file should be converted to. Choose one of: ", paste(out_formats, collapse = ", "), ".")
)
parser$add_argument(
    '--sample-id',
    type="character",
    dest='sample_id',
    action = "store",
    default = NULL,
    help = "Sample ID of the given input file."
)
parser$add_argument(
    '--group-name',
    type="character",
    dest='group_name',
    action = "store",
    default = NULL,
    help = "Name of the group which the given input files are from. A new column named by this group_name will be added in anndata.obs"
)
parser$add_argument(
    '--group-value',
    type="character",
    dest='group_value',
    action = "store",
    default = NULL,
    help = "Value of the group which the given input files are from. The group_name column in anndata.obs will be populated with this group_value."
)
parser$add_argument(
    '--tag-cell-with-sample-id',
    type="character",
    dest='tag_cell_with_sample_id',
    action = "store",
    default = TRUE,
    help = "Append the sample ID to the cell barcodes."
)
parser$add_argument(
    '--remove-10x-gem-well',
    type="character",
    dest='remove_10x_gem_well',
    action = "store",
    default = FALSE,
    help = "If tag_cell_with_sample_id is passed, remove the GEM well number from the barcode."
)
parser$add_argument(
    '--seurat-assay',
    type="character",
    dest='seurat_assay',
    action = "store",
    default = "RNA",
    help = "[Seurat] The assay name which the --seurat-layer will get it from."
)
parser$add_argument(
    '--seurat-main-layer',
    type="character",
    dest='seurat_main_layer',
    action = "store",
    default = "counts",
    help = "[Seurat] The layer name of array to put as main matrix."
)
parser$add_argument(
    '--seurat-reset',
    dest='seurat_reset',
    type="character",
    action = "store",
    default = TRUE,
    help = "[Seurat] If true, the following slots will be removed: @graphs, @neighbors, @reductions, $[[args$`seurat-assay`]]@data, $[[args$`seurat-assay`]]@scale.data."
)
parser$add_argument(
    '--sce-main-layer',
    type="character",
    dest='sce_main_layer',
    action = "store",
    default = "counts",
    help = "[SingleCellExperiment] The layer name of array to put as main matrix. Used only when --input-format sce_rds."
)

args <- parser$parse_args()
print(args)

isTrue <- function(x) {
	x <- tolower(x = x)
	return (
		x == FALSE | x == F | x == 't' | x == 'true' | x == 'y' | x == 'yes'
	)
}

UpdateSeuratCellMetadata <- function(seurat, args) {
    # Add sample ID as obs entry
	seurat <- AddMetaData(
		object = seurat,
		metadata = as.character(x = args$`sample_id`),
		col.name = 'sample_id'
    )

    # Tag cell with sample ID
    if(isTrue(x = args$`tag_cell_with_sample_id`)) {
        if(isTrue(x = args$`remove_10x_gem_well`)) {
            new.names <- gsub(
                pattern = "-([0-9]+)$",
                replace = paste0("-", args$`sample_id`),
                x = colnames(x = seurat)
            )
        } else {
            new.names <- paste0(colnames(x = seurat), "___", args$`sample_id`)
        }
		seurat <- Seurat::RenameCells(
			object = seurat,
			new.names = new.names
		)
    }
    
    # Add group_value as obs entry with group_name as column name
    if(!is.null(args$`group_name`) & !is.null(args$`group_value`)) {
        seurat <- AddMetaData(
            object = seurat,
            metadata = as.character(x = args$`group_value`),
            col.name = as.character(x = args$`group_name`)
        )
    }
    return (seurat)
}

UpdateSCECellMetadata <- function(sce, args) {
    # Add sample ID as colData entry
	col_data <- SummarizedExperiment::colData(x = sce)
    col_data$sample_id <- args$`sample_id`
    
    # Tag cell with sample ID
    if(isTrue(x = args$`tag_cell_with_sample_id`)) {
        if(isTrue(x = args$`remove_10x_gem_well`)) {
            new.names <- gsub(
                pattern = "-([0-9]+)$",
                replace = paste0("-", args$`sample_id`),
                x = colnames(x = sce)
            )
        } else {
            new.names <- paste0(colnames(x = sce), "___", args$`sample_id`)
        }
        colnames(x = sce) <- new.names
    }
    # Add group_value as obs entry with group_name as column name
    if(!is.null(args$`group_name`) & !is.null(args$`group_value`)) {
        col_data[[args$`group_name`]] <-args$`group_value`
    }
    SummarizedExperiment::colData(x = sce) <- col_data
    return (sce)
}

# Define the arguments properly
FILE_PATH_IN <- args$`input`
FILE_NAME_SPLITTED <- strsplit(
	x = basename(path = args$`output`),
	split = "\\."
)[[1]]
FILE_PATH_OUT_BASENAME <- paste0(FILE_NAME_SPLITTED[1], ".", FILE_NAME_SPLITTED[2])
INPUT_FORMAT <- args$`input_format`
OUTPUT_FORMAT <- args$`output_format`

if(INPUT_FORMAT == 'seurat_rds' & OUTPUT_FORMAT == 'h5ad') {
    seurat <- tryCatch({
      	readRDS(file = FILE_PATH_IN)
    }, warning = function(w) {
    }, error = function(e) {
      	stop("VSN ERROR: Cannot read the given Rds file.")
    }, finally = {})

    if(class(x = seurat) != "Seurat") {
      	stop("VSN ERROR: The object contained in the Rds file is not a Seurat object.")
    }
    # Reset Seurat slots
    if(isTrue(x = args$`seurat_reset`)) {
		print("Resetting @graphs in seurat object...")
		seurat@graphs <- list()
		print("Resetting @reductions in seurat object...")
		seurat@reductions <- list()
		print("Resetting @neighbors in seurat object...")
		seurat@neighbors <- list()
        # Set @assays`<assay>`@data to @assays`<assay>`@counts (as in CreateSeuratObject
        # @assays`<assay>`@data cannot be set to empty matrix otherwise row.names(object) will be NULL
		print(paste0("Setting @assays$",args$`seurat_assay`,"@data to @assays$",args$`seurat_assay`,"@counts (as in CreateSeuratObject) in seurat object..."))
		seurat@assays[[args$`seurat_assay`]]@data <- seurat@assays[[args$`seurat_assay`]]@counts
		print(paste0("Resetting @assays$",args$`seurat_assay`,"@scale.data in seurat object..."))
		seurat@assays[[args$`seurat_assay`]]@scale.data <- matrix(ncol = 0, nrow = 0)
    }

    # Update Seurat object cell metadata
    seurat <- UpdateSeuratCellMetadata(seurat, args)
    # Check if all meta.data columns are 1-dimensional (i.e.: dim(x = ...) should not return NULL)
    are_metadata_cols_dim_not_null <- do.call(
		what="c", 
		args=lapply(
			X=colnames(x = seurat@meta.data),
			FUN=function(x) { 
			    !is.null(x=dim(x = seurat@meta.data[[x]]))
			}
		)
    )
    if(any(are_metadata_cols_dim_not_null)) {
        metadata_cols_dim_not_null_colnames <- colnames(x = seurat@meta.data[, which(x = are_metadata_cols_dim_not_null)])
        stop(paste0("VSN ERROR: Some columns (", paste(metadata_cols_dim_not_null_colnames, collapse=" and "),") from the given Seurat object in the meta.data slot are not 1-dimensional."))
    }

    # Sort genes
    seurat <- seurat[sort(x = row.names(x = seurat)),]
    
    sceasy::convertFormat(
		seurat,
		from="seurat",
		to="anndata",
		outFile=paste0(FILE_PATH_OUT_BASENAME, ".h5ad"),
		main_layer = args$`seurat_main_layer`,
		assay = args$`seurat_assay`,
		drop_single_values = FALSE
    )
} else if(INPUT_FORMAT == 'sce_rds' & OUTPUT_FORMAT == 'h5ad') {
    sce <- tryCatch({
      	readRDS(file = FILE_PATH_IN)
    }, warning = function(w) {
    }, error = function(e) {
      	stop("VSN ERROR: Cannot read the given Rds file.")
    }, finally = {})

    if(class(x = sce) != "SingleCellExperiment") {
      	stop("VSN ERROR: The object contained in the Rds file is not a SingleCellExperiment object.")
    }

    # Set/update row.names with gene symbols
    row_data <- SummarizedExperiment::rowData(x = sce)
    if("Symbol" %in% colnames(x = row_data)) {
	    row.names(x = sce) <- row_data$Symbol
    }

    # Update SingleCellExperiment object cell metadata
    sce <- UpdateSCECellMetadata(sce, args)

    # Sort genes
    sce <- sce[sort(x = row.names(x = sce)),]
    sceasy::convertFormat(
		sce,
		from="sce",
		to="anndata",
		outFile=paste0(FILE_PATH_OUT_BASENAME, ".h5ad"),
		main_layer = args$`sce_main_layer`,
		drop_single_values = FALSE
    )
} else if(INPUT_FORMAT == '10x_cellranger_mex' & OUTPUT_FORMAT == 'sce_rds') {
    sce <- DropletUtils::read10xCounts(
      	samples = FILE_PATH_IN
    )

    # Set/update row.names with gene symbols
    row_data <- SummarizedExperiment::rowData(x = sce)
    if("Symbol" %in% colnames(x = row_data)) {
	    row.names(x = sce) <- row_data$Symbol
    }

    # Set col.names with barcode ID
    colnames(x = sce) <- SummarizedExperiment::colData(x = sce)$Barcode

    # Update SingleCellExperiment object cell metadata
    sce <- UpdateSCECellMetadata(sce, args)
    # Sort genes
    sce <- sce[sort(x = row.names(x = sce)),]
    saveRDS(
		object = sce,
		file = paste0(FILE_PATH_OUT_BASENAME, ".SCE.Rds"),
		compress = TRUE
    )
} else if(INPUT_FORMAT == '10x_cellranger_mex' & OUTPUT_FORMAT == 'seurat_rds') {
    if (grepl("\\.h5", FILE_PATH_IN)) {
        counts_matrix <- Seurat::Read10X_h5(FILE_PATH_IN)
    } else {
        counts_matrix <- Seurat::Read10X(FILE_PATH_IN)
    }
    # Create seurat object without any filtering
    # TODO: handle multi modal data. Make sure we only get the RNA
    seurat <- Seurat::CreateSeuratObject(
        counts_matrix,
        min.cells = 0,
        min.features = 0,
        project = args$`sample_id`
    )
    saveRDS(
        object = seurat,
        file = paste0(FILE_PATH_OUT_BASENAME, ".seurat.Rds"),
        compress = TRUE
    )
} else {
    stop(paste0("VSN ERROR: File format conversion ", INPUT_FORMAT," --> ", OUTPUT_FORMAT," hasn't been implemented yet.)"))
}

print("Done!")

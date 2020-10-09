#!/usr/bin/env Rscript

in_formats = c("seurat_rds")
out_formats = c("h5ad")

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
    '--tag-cell-with-sample-id',
    type="character",
    dest='tag_cell_with_sample_id',
    action = "store",
    default = TRUE,
    help = "Sample ID of the given input file."
)
parser$add_argument(
    '--seurat-assay',
    type="character",
    dest='seurat_assay',
    action = "store",
    default = "RNA",
    help = "The assay name which the --seurat-layer will get it from."
)
parser$add_argument(
    '--seurat-main-layer',
    type="character",
    dest='seurat_main_layer',
    action = "store",
    default = "counts",
    help = "The layer name of array to put as main matrix."
)
parser$add_argument(
    '--seurat-reset',
    dest='seurat_reset',
    type="character",
    action = "store",
    default = TRUE,
    help = "If true, the following slots will be removed: @graphs, @neighbors, @reductions, $[[args$`seurat-assay`]]@data, $[[args$`seurat-assay`]]@scale.data."
)

args <- parser$parse_args()

isTrue <- function(x) {
	x <- tolower(x = x)
	return (
		x == FALSE | x == F | x == 't' | x == 'true' | x == 'y' | x == 'yes'
	)
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
      	stop("Cannot read the given Rds file.")
    }, finally = {})

    if(class(x = seurat) != "Seurat") {
      	stop("The object contained in the Rds file is not a Seurat object.")
    }
    # Tag cell with sample ID
    if(isTrue(x = args$`tag_cell_with_sample_id`)) {
		new.names <- gsub(
			pattern = "-([0-9]+)$",
			replace = paste0("-", args$`sample_id`),
			x = colnames(x = seurat)
		)
		seurat <- Seurat::RenameCells(
			object = seurat,
			new.names = new.names
		)
    }
    # Reset Seurat slots
    if(isTrue(x = args$`seurat_reset`)) {
		print("Resetting @graphs in seurat object...")
		seurat@graphs<-list()
		print("Resetting @reductions in seurat object...")
		seurat@reductions<-list()
		print("Resetting @neighbors in seurat object...")
		seurat@neighbors<-list()
		print(paste0("Resetting @assays$",args$`seurat_assay`,"@data in seurat object..."))
		seurat@assays[[args$`seurat_assay`]]@data<-matrix(ncol = 0, nrow = 0)
		print(paste0("Resetting @assays$",args$`seurat_assay`,"@scale.data in seurat object..."))
		seurat@assays[[args$`seurat_assay`]]@scale.data<-matrix(ncol = 0, nrow = 0)
    }
    # Add sample ID as obs entry
	seurat <- AddMetaData(
		object = seurat,
		metadata = as.character(x = args$`sample_id`),
		col.name = 'sample_id'
    )
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
    if(any(!are_metadata_cols_dim_not_null)) {
        metadata_cols_dim_not_null_colnames <- colnames(x = seurat@meta.data[, which(x = are_metadata_cols_dim_not_null)])
        stop(paste0("Some columns (", paste(metadata_cols_dim_not_null_colnames, collapse=" and "),") from the given Seurat object in the meta.data slot are not 1-dimensional."))
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
} else if(INPUT_FORMAT == '10x_cellranger_mex' & OUTPUT_FORMAT == 'sce_rds') {
} else {
  stop(paste0("File format conversion ", INPUT_FORMAT," --> ", OUTPUT_FORMAT," hasn't been implemented yet.)"))
}

print("Done!")
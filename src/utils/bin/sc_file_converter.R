#!/usr/bin/env Rscript

in_formats = c("seurat_rds")
out_formats = c("h5ad")

library("optparse")
parser <- OptionParser(
  prog = "sc_file_converter.R",
  description = "Conversion script between single-cell data formats."
)
parser <- add_option(
  parser,
  c("-v", "--input-file"),
  action = "store",
  default = NULL,
  help = "Input file [default]"
)
parser <- add_option(
  parser,
  c("-w", "--output-file"),
  action = "store",
  default = NULL,
  help = "Output file [default]"
)
parser <- add_option(
  parser,
  c("-i", "--input-format"),
  action = "store",
  default = NULL,
  help = paste0("Input format of the file to be converted. Choose one of: ", paste(in_formats, collapse = ", "), ".")
)
parser <- add_option(
  parser,
  c("-o", "--output-format"),
  action = "store",
  default = NULL,
  help = paste0("Output format which the file should be converted to. Choose one of: ", paste(out_formats, collapse = ", "), ".")
)

parser <- add_option(
  parser,
  c("-s", "--sample-id"),
  action = "store",
  default = NULL,
)
parser <- add_option(
  parser,
  c("-t", "--tag-cell-with-sample-id"),
  action = "store",
  default = TRUE,
  help = "Tag each cell with the given sample_id."
)
parser <- add_option(
  parser,
  c("-a", "--seurat-assay"),
  action = "store",
  default = "RNA",
  help = "The assay name which the --seurat-layer will get it from."
)
parser <- add_option(
  parser,
  c("-l", "--seurat-main-layer"),
  action = "store",
  default = "counts",
  help = "The layer name of array to put as main matrix."
)
parser <- add_option(
  parser,
  c("-r", "--seurat-reset"),
  action = "store",
  default = TRUE,
  help = "If true, the following slots will be removed: @graphs, @neighbors, @reductions, $[[args$`seurat-assay`]]@data, $[[args$`seurat-assay`]]@scale.data."
)

args <- parse_args(parser)

isTrue <- function(x) {
  x <- tolower(x = x)
  return (
    x == FALSE | x == F | x == 't' | x == 'true' | x == 'y' | x == 'yes'
  )
}

# Define the arguments properly
FILE_PATH_IN <- args$`input-file`
FILE_NAME_SPLITTED <- strsplit(x = basename(path = args$`output-file`), split = "\\.")[[1]]
FILE_PATH_OUT_BASENAME <- paste0(FILE_NAME_SPLITTED[1], ".", FILE_NAME_SPLITTED[2])
INPUT_FORMAT <- args$`input-format`
OUTPUT_FORMAT <- args$`output-format`

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
  if(isTrue(x = args$`tag-cell-with-sample-id`)) {
    new.names <- gsub(
      pattern = "-([0-9]+)$",
      replace = paste0("-", args$`sample-id`),
      x = colnames(x = seurat)
    )
    seurat <- Seurat::RenameCells(
      object = seurat,
      new.names = new.names
    )
  }
  if(isTrue(x = args$`seurat-reset`)) {
    print("Resetting @graphs in seurat object...")
    seurat@graphs<-list()
    print("Resetting @reductions in seurat object...")
    seurat@reductions<-list()
    print("Resetting @neighbors in seurat object...")
    seurat@neighbors<-list()
    print(paste0("Resetting @assays$",args$`seurat-assay`,"@data in seurat object..."))
    seurat@assays[[args$`seurat-assay`]]@data<-matrix(ncol = 0, nrow = 0)
    print(paste0("Resetting @assays$",args$`seurat-assay`,"@scale.data in seurat object..."))
    seurat@assays[[args$`seurat-assay`]]@scale.data<-matrix(ncol = 0, nrow = 0)
  }
  # Add sample ID as obs entry
  seurat <- AddMetaData(
    object = seurat,
    metadata = as.character(x = args$`sample-id`),
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
  
  sceasy::convertFormat(
    seurat,
    from="seurat",
    to="anndata",
    outFile=paste0(FILE_PATH_OUT_BASENAME, ".h5ad"),
    main_layer = args$`seurat-main-layer`,
    assay = args$`seurat-assay`,
    drop_single_values = FALSE
  )
} else {
  stop(paste0("File format conversion ", INPUT_FORMAT," --> ", OUTPUT_FORMAT," hasn't been implemented yet.)"))
}

print("Done!")
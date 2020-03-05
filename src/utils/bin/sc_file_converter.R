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
  c("-l", "--seurat-layer"),
  action = "store",
  default = "counts",
  help = "The layer name of array to put as main matrix."
)

args <- parse_args(parser)

isTrue <- function(x) {
  x <- tolower(x = x)
  return (
    x | x == 't' | x == 'true' | x == 'y' | x == 'yes'
  )
}

# Define the arguments properly
FILE_PATH_IN = args$input
FILE_PATH_OUT_BASENAME = strsplit(x = args$output, split = "\\.")[[1]][1]
INPUT_FORMAT = args$`input-format`
OUTPUT_FORMAT = args$`output-format`

if(INPUT_FORMAT == 'seurat_rds' & OUTPUT_FORMAT == 'h5ad') {
  seurat <- tryCatch({
    readRDS(file = FILE_PATH_IN)
  }, warning = function(w) {
  }, error = function(e) {
    stop("Cannot read the given Rds file.")
  }, finally = {}

  if(class(seurat) != "Seurat") {
    stop("The object contained in the Rds file is not a Seurat object.")
  }
  if(isTrue(x = args$`tag-cell-with-sample-id)) {
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
  sceasy::convertFormat(
    seurat,
    from="seurat",
    to="anndata",
    assay = args$`seurat-assay`,
    main_layer = args$`seurat-layer`,
    outFile=paste0(FILE_PATH_OUT_BASENAME, ".h5ad")
  )
} else {
  stop(paste0("File format conversion ", INPUT_FORMAT," --> ", OUTPUT_FORMAT," hasn't been implemented yet.)"))
}

print("Done!")
#!/usr/bin/env Rscript

print("##################################################")
print("# Harmony: Algorithm for single cell integration #")
print("##################################################")

# Loading dependencies scripts

library("optparse")
parser <- OptionParser(
  prog = "run_harmony.R",
  description = "Scalable integration of single cell RNAseq data for batch correction and meta analysis"
)
parser <- add_option(
  parser,
  c("-i", "--input-file"),
  action = "store",
  default = NULL,
  help = "Input file [default]"
)
parser <- add_option(
  parser,
  c("-a", "--vars-use"),
  action = "store",
  default = NULL,
  help = "If meta_data is dataframe, this defined which variable(s) to remove (character vector)."
)
parser <- add_option(
  parser,
  c("-p", "--do-pca"),
  action = "store",
  default = FALSE,
  help = "Whether to perform PCA on input matrix."
)
parser <- add_option(
  parser,
  c("-o", "--output-prefix"),
  action = "store",
  default = "foo",
  help="Prefix path to save output files. [default %default]"
)

args <- parse_args(parser)

cat("Parameters: \n")
print(args)

if(is.null(args$`vars-use`)) {
	stop("The parameter --vars-use has to be set.")
}

input_ext <- tools::file_ext(args$`input-file`)

if(input_ext == "h5ad") {
  # Current fix until https://github.com/satijalab/seurat/issues/2485 is fixed
  file <- hdf5r::h5file(filename = args$`input-file`, mode = 'r')
  if(!("X_pca" %in% names(x = file[["obsm"]]))) {
    stop("HI")
  }
  obs <- file[['obs']][]
  pca_embeddings <- t(x = file[["obsm"]][["X_pca"]][,])
  row.names(x = pca_embeddings) <- obs$index
  colnames(x = pca_embeddings) <- paste0("PCA_", seq(from = 1, to = ncol(x = pca_embeddings)))
  metadata <- obs
  # seurat <- Seurat::ReadH5AD(file = args$`input-file`)
  # if(!("pca" %in% names(seurat@reductions)) || is.null(x = seurat@reductions$pca))
  #   stop("Expects a PCA embeddings data matrix but it does not exist.")
  # data <- seurat@reductions$pca
  # pca_embeddings <- data@cell.embeddings
  # metadata <- seurat@meta.data
} else {
  stop(paste0("Unrecognized input file format: ", input_ext, "."))
}

print(paste0("PCA embeddings matrix has ", dim(x = data)[1], " rows, ", dim(x = data)[2], " columns."))

if(sum(args$`vars-use` %in% colnames(x = metadata)) != length(x = args$`vars-use`)) {
	stop("Some argument value from the parameter(s) --vars-use are not found in the metadata.")
}

# Run Harmony
# Expects PCA matrix (Cells as rows and PCs as columns.)
harmony_embeddings <- harmony::HarmonyMatrix(data_mat = pca_embeddings
                                             , meta_data = metadata
                                             , vars_use = args$`vars-use`
                                             , do_pca = args$`do-pca`
                                             , verbose = FALSE
)

# Save the results

## PCA corrected embeddings

write.table(
	x = harmony_embeddings,
	file = paste0(args$`output-prefix`, ".tsv"),
	quote = FALSE,
	sep = "\t",
	row.names = TRUE,
	col.names = NA
)

#!/usr/bin/env Rscript

print("#############################################################")
print("# Harmony: Algorithm for single cell integration            #")
print('# GitHub: https://github.com/immunogenomics/harmony         #')
print('# Paper: https://www.nature.com/articles/s41592-019-0619-0  #')
print("#############################################################")

# Loading dependencies scripts
library("argparse")

parser <- ArgumentParser(description='Scalable integration of single cell RNAseq data for batch correction and meta analysis')
parser$add_argument(
    'input',
    metavar='INPUT',
    type="character",
    help='Input file [default]'
)
parser$add_argument(
    '--output-prefix',
    type="character",
    dest='output_prefix',
    default = "foo",
  	help="Prefix path to save output files. [default %default]"
)
parser$add_argument(
    '--seed',
    type="character",
    dest='seed',
    default=617,
    help='Seed. [default %default]'
)
parser$add_argument(
    "--vars-use",
    type="character",
    dest='vars_use',
    action="append",
	  default=NULL,
    help='If meta_data is dataframe, this defined which variable(s) to remove (character vector).'
)
parser$add_argument(
	'--do-pca',
	type="logical",
	dest='do_pca',
	action="store",
	default=FALSE,
	help='Whether to perform PCA on input matrix.'
)
parser$add_argument(
    '--theta',
    type="double",
    dest='theta',
    default=NULL,
    help='Diversity clustering penalty parameter. Specify for each variable in vars_use Default theta=2. theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters. [default %default]'
)
parser$add_argument(
    '--lambda',
    type="double",
    dest='lambda',
    default=NULL,
    help='Ridge regression penalty parameter. Specify for each variable in vars_use. Default lambda=1. Lambda must be strictly positive. Smaller values result in more aggressive correction. [default %default]'
)
parser$add_argument(
    '--epsilon-harmony',
    type="double",
    dest='epsilon_harmony',
    default=1e-04,
    help='Convergence tolerance for Harmony. Set to -Inf to never stop early. [default %default]'
)


args <- parser$parse_args()

if(args$epsilon_harmony < 0) {
  args$epsilon_harmony <- -Inf
  print("Setting epsilon.harmony argument to -Inf...")
}

cat("Parameters: \n")
print(args)

if(is.null(args$vars_use)) {
	stop("The parameter --vars-use has to be set.")
}

# Required by irlba::irlba (which harmony depends on) for reproducibility
if(!is.null(args$seed)) {
  set.seed(args$seed)
} else {
  warnings("No seed is set, this will likely give none reproducible results.")
}

# Required for reproducibility in case numeric parameters are passed (e.g.: theta, lambda)
args <- lapply(X = args, FUN = function(arg) {
  if(is.numeric(x = arg)) {
    if(arg %% 1 == 0) {
      return (as.integer(x = arg))
    } else {
      return (arg)
    }
  }
  return (arg)
})


input_ext <- tools::file_ext(args$input)

if(input_ext == "h5ad") {
  # Current fix until https://github.com/satijalab/seurat/issues/2485 is fixed
  file <- hdf5r::h5file(filename = args$input, mode = 'r')
  if(!("X_pca" %in% names(x = file[["obsm"]]))) {
    stop("X_pca slot is not found in the AnnData (h5ad).")
  }
  obs <- file[['obs']][]
  pca_embeddings <- t(x = file[["obsm"]][["X_pca"]][,])
  row.names(x = pca_embeddings) <- obs$index
  colnames(x = pca_embeddings) <- paste0("PCA_", seq(from = 1, to = ncol(x = pca_embeddings)))
  metadata <- obs
  # seurat <- Seurat::ReadH5AD(file = args$input)
  # if(!("pca" %in% names(seurat@reductions)) || is.null(x = seurat@reductions$pca))
  #   stop("Expects a PCA embeddings data matrix but it does not exist.")
  # data <- seurat@reductions$pca
  # pca_embeddings <- data@cell.embeddings
  # metadata <- seurat@meta.data
} else {
  stop(paste0("Unrecognized input file format: ", input_ext, "."))
}

print(paste0("PCA embeddings matrix has ", dim(x = pca_embeddings)[1], " rows, ", dim(x = pca_embeddings)[2], " columns."))

if(sum(args$vars_use %in% colnames(x = metadata)) != length(x = args$vars_use)) {
	stop("Some argument value from the parameter(s) --vars-use are not found in the metadata.")
}

print(paste0("Batch variables used for integration: ", paste0(args$vars_use, collapse=", ")))

# Run Harmony
# Expects PCA matrix (Cells as rows and PCs as columns.)
harmony_embeddings <- harmony::HarmonyMatrix(
  data_mat = pca_embeddings,
  meta_data = metadata,
  vars_use = args$vars_use,
  do_pca = args$do_pca,
  theta = args$theta,
  lambda = args$lambda,
  epsilon.harmony = args$epsilon_harmony,
  verbose = FALSE
)

# Save the results

## PCA corrected embeddings

write.table(
	x = harmony_embeddings,
	file = paste0(args$output_prefix, ".tsv"),
	quote = FALSE,
	sep = "\t",
	row.names = TRUE,
	col.names = NA
)

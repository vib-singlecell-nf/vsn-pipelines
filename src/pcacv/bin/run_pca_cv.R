#!/usr/bin/env Rscript

print("##################################################")
print("# Calculating Optimal Number of PCs using PCA CV #")
print("##################################################")

# This a R implementation of the algorithm described at:
# - https://stats.stackexchange.com/questions/93845/how-to-perform-leave-one-out-cross-validation-for-pca-to-determine-the-number-of

library("optparse")
parser <- OptionParser(
  prog = "dim_reduction_pca_cv.R",
  description = "Perform PCA through cross-validation to find the optimal number of principal components."
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
  c("-a", "--accessor"),
  action = "store",
  default = NULL,
  help = "Character defining the path to access the data matrix in S4 object."
)
parser <- add_option(
  parser,
  c("-y", "--use-variable-features"),
  action = "store",
  default = TRUE,
  help = "Use the highly features only otherwise use all available features. [default %default]"
)
parser <- add_option(
  parser,
  c("-p", "--n-pc-fallback"),
  action = "store",
  default = 0,
  help="Minimum principal components to return. If no optimal number of PCs is found, this minimum number of PCs will be returned. [default %default]"
)
parser <- add_option(
  parser,
  c("-k", "--k-fold"),
  action = "store",
  default = 10,
  help="Set value k for the k-fold cross-validation. [default %default]"
)
parser <- add_option(
  parser,
  c("-f", "--from-n-pc"),
  action = "store",
  default = 2,
  help="Number of principal components to start with. [default %default]"
)
parser <- add_option(
  parser, c("-t", "--to-n-pc"),
  action = "store",
  default = 150,
  help = "Number of principal components to compute [default %default]"
)
parser <- add_option(
  parser, 
  c("-b", "--by-n-pc"),
  action = "store",
  default = 5,
  help="Number of principal components to increment by. [default %default]"
)
parser <- add_option(
  parser,
  c("-m", "--max-iters"), 
  action="store", 
  default = 500, 
  help="Maximum number of iteration. [default %default]")
parser <- add_option(
  parser, 
  c("-s", "--seed"), 
  action = "store", 
  default = NULL,
  help="Seed. [default %default]"
)
parser <- add_option(
  parser,
  c("-c", "--n-cores"),
  action = "store",
  default = 1,
  help="Number of cores to use. [default %default]"
)
parser <- add_option(
  parser,
  c("-d", "--default-svd"),
  action = "store",
  default = FALSE,
  help="Perform SVD instead of approximate PCA. [default %default]"
)
parser <- add_option(
  parser,
  c("-v", "--verbose"),
  action = "store",
  default = TRUE,
  help="Print some informational messages. [default %default]"
)
parser <- add_option(
  parser,
  c("-o", "--output-prefix"),
  action = "store",
  default = "foo",
  help="Prefix path to save output files (table of the PCA CV, PCA PRESS plot, optimal number of PCs). [default %default]"
)

args <- parse_args(parser)

# Dollar sign needs to be slashed ($ -> \$ when running in BATCH MODE)
GetDataFromS4ObjectByAccessPath <- function(object, path) {
  acc <- matrix(data = strsplit(x = path, "(?=[@$])", perl = TRUE)[[1]], nrow = 2)
  S4Parse <- function(object, accessor, i = 1) {
    if(i > ncol(x = accessor)) 
      return (object)
    o <- switch(EXPR = acc[1, i],
                "@" = slot(object = object, name = acc[2, i]), 
                "$" = object[[acc[2, i]]])
    return (S4Parse(object = o, accessor = accessor, i = i + 1))
  }
  return (S4Parse(object = object, acc))
}

ConvertH5GroupToDataFrame <- function(h5.group) {
  h5_external_names <- names(h5.group)[Vectorize(function(x) !startsWith(x = x, prefix = "_"))(names(x = h5.group))]
  df <- do.call(what = "cbind", args = lapply(X = h5_external_names, FUN = function(h5_name) {
    return (h5.group[[h5_name]][])
  }))
  colnames(x = df) <- h5_external_names
  return (as.data.frame(x = df))
}

ExtractDataFromH5Object <- function(h5.object) {
  if(h5.object$get_obj_type() == "H5I_DATASET")
    return (h5.object[])
  # Handle the case for H5AD build with python package anndata>=0.7.1
  if(h5.object$get_obj_type() == "H5I_GROUP")
    return (
      ConvertH5GroupToDataFrame(h5.group = h5.object)
    )
  stop("Unrecognized H5 type.")
}

#'
#' @title RunPCACV
#' @description       Return the PRESS errors of the PCA cross-validation 
#' @param data        M-by-N matrix where M = features, N = cells. Expects a feature-scaled data matrix.
#' @param k           Number of partitions in cross-validation.
#' @param from        Number of principal components to start with.
#' @param to          Number of principal components to compute.
#' @param by          Number of principal components to increment by.
#' @param maxit       Maximum number of iterations.
#' @param seed        Seed number for reproducibility.
#' @param n.cores     Number of cores to use.
#' @param verbose     If TRUE, prints status messages during the computation.
#' @param default.svd If TRUE, perform base::svd instead of approximate PCA.
#' @return data.frame containing of the PRESS error for each number of principal components tested.
RunPCACV <- function(
	data, 
	k = 10,
	from = 2,
	to = 150,
	by = 5,
	maxit = 200,
	seed = NULL,
	n.cores = 1,
	verbose = T,
	default.svd = F) {
		# Load the libraries required for parallelization
		library(doFuture)
		# Increase max size of global objects passed to an item in the future package
		options(future.globals.maxSize= 2000*1024^2)
		library(doRNG)

		if(is.null(seed)) {
			stop("No seed is set, this will likely give none reproducible results. Please set one.")
		}
		# Required by irlba::irlba for reproducibility
		set.seed(seed)

		# Setup the parallelization
		print(paste0("Number of workers used: ", n.cores))
		registerDoFuture()
		plan(multiprocess, workers = n.cores)
		registerDoRNG(as.numeric(seed))

		if(default.svd) {
			print("Default singular value decomposition (SVD) is used for estimating the number PCs.") 
		}

		pc <- seq(from = from, to = to, by = by)

		if(length(x = pc) > 50) {
			stop("Too many PCs to test please adjust: --k-fold, --from-n-pc, --to-n-pc, --by-n-pc.")
		}

		# Init the error matrices
		error <- matrix(0, nrow = length(c(1:k)), ncol = length(x = pc))

		# Partition the data into K-folds
		print(paste0(k,"-fold paritioning..."))
		data_kfold <- dismo::kfold(x = Matrix::t(x = data), k=k)

		# Perform PCA cross-validation
		for(i in c(1:k)) {
			print(paste0("Processing ", i, "th fold..."))
			data_train <- Matrix::t(x = data[, data_kfold!=i])
			data_test <- Matrix::t(x = data[, data_kfold==i])
			print("Running SVD...")
			if(!default.svd) {
				print("...fast truncated SVD")
				pca_results <- irlba::irlba(
					A = data_train
					, nv = to
					, maxit = maxit
					, verbose = verbose
				)
			} else {
				print("...regular SVD")
				pca_results <- base::svd(
					x = data_train
					, nv = to
				)
			}
			print("==========================================")
			gl <- pca_results$v
			res <- foreach(j=1:length(x = pc), .combine='rbind') %dopar% {
			# for(j in 1:length(x = pc)) {
				print(paste0("...with ", pc[j], " PCs."))
				P <- gl[,1:pc[j]]%*%Matrix::t(gl[,1:pc[j]])
				err <- data_test %*% (diag(x = dim(x = P)[1]) - P + diag(x = diag(x = P)))
				return (
					data.frame("i"=i, "j"=j, err=sum(err^2))
				)
			}
			apply(X = res, MARGIN = 1, FUN = function(r) {
				error[r[["i"]],r[["j"]]] <<- r[["err"]]
			})
			# Free memory
			# https://stat.ethz.ch/R-manual/R-patched/library/base/html/gc.html
			rm(data_train)
			rm(data_test)
			rm(pca_results)
			rm(gl)
			rm(res)
			gc()
		}
		errors <- colSums(x = error)
		return (data.frame("PC"=pc,"error"=log(x = errors)))
	}
	cat("Parameters: \n")
	print(args)

	input_ext <- tools::file_ext(args$`input-file`)

	hv_column_names <- hv_column_names <- c("highly.variable","highly_variable")
	if(input_ext == "loom") {
		print("Input type: Loom...")
		# Get the data matrix from loom 
		warnings("ASSUMPTION: Expecting a scaled data matrix in the main layer of the loom file!")
		loom  <- SCopeLoomR::open_loom(file.path = args$`input-file`)
		data <- SCopeLoomR::get_dgem(loom = loom)
	} else if(input_ext == "Rds") {
		print("Input type: Rds...")
		object <- readRDS(file = args$`input-file`)
		# Check if the accessor argument is NULL
		if(is.null(args$`accessor`))
			stop("The --accessor argument need to be set.")
		data <- GetDataFromS4ObjectByAccessPath(
			object = object
			, path = args$`accessor`
		)
		meta_features <- object@assays$RNA@meta.features
		# FIXME: what do we do when there are more than 1 variable feature column in the data 
		# e.g. after running SCT and FindVariableFeatures on the same dataset. 
		# This is highly unlikely, but we should still catch it here.
		variable_column <- grep(".variable", colnames(meta_features), value = T)[[1]]
		hv_column_names <- variable_column
	} else if(input_ext == "h5ad") {
		# Current fix until https://github.com/satijalab/seurat/issues/2485 is fixed
		file <- hdf5r::h5file(filename = args$`input-file`, mode = 'r')
		if (is(object = file[['X']], class2 = 'H5Group')) {
			data <- Seurat::as.sparse(x = file[['X']])
		} else {
			data <- file[['X']][, ]
		}
		# x will be an S3 matrix if X was scaled, otherwise will be a dgCMatrix
		# Pull cell- and feature-level metadata
		obs <- ExtractDataFromH5Object(file[['obs']])
		meta_features <- ExtractDataFromH5Object(file[['var']])
		rownames(x = data) <- rownames(x = meta_features) <- meta_features$index
		colnames(x = data) <- rownames(x = obs) <- obs$index
		file$close_all()
		# seurat <- Seurat::ReadH5AD(file = args$`input-file`)
		# if(!methods::.hasSlot(object = seurat@assays$RNA, name = "scale.data") || is.null(x = slot(object = seurat@assays$RNA, name = "scale.data")))
		#   stop("Expects a feature-scaled data matrix but it does not exist.")
		# if(is.null(args$`accessor`))
		#   stop("The --accessor argument need to be set.")
		# data <- GetDataFromS4ObjectByAccessPath(
		#   object = seurat
		#   , path = args$`accessor`
		# )
	} else {
		stop(paste0("Unrecognized input file format: ", input_ext, "."))
	}

	use_variable_features <- as.logical(args$`use-variable-features`)

	if(use_variable_features) {
		hv_column_mask <- hv_column_names %in% colnames(x = meta_features)
		if(sum(x = hv_column_mask) == 0) {
			stop("Cannot subset the matrix with the 'highly variable features since 'highly.variable' or 'highly_variable' is not present does not exist in meta.features ('seurat@assays$RNA@meta.features') data.frame.'")
		} else {
			hv_column_name <- hv_column_names[hv_column_mask]
			meta_features[[hv_column_name]][is.na(x = meta_features[[hv_column_name]])] <- FALSE
			hvf <- row.names(x = meta_features)[meta_features[[hv_column_name]]]
		}
		if(all(hvf == row.names(x = data))) {
			print("Input contains already data with highly variable features.")
			print("Skipping the subsetting.")
		} else {
			print(paste0("Subsetting matrix with highly variable features... "))
			data <- data[hvf,]
		}
	}

print(paste0("Data matrix has ", dim(x = data)[1], " rows (features), ", dim(x = data)[2], " columns (observations)."))

k <- args$`k-fold`
from_npcs <- args$`from-n-pc`
to_npcs <- args$`to-n-pc`
by_npcs <- args$`by-n-pc`

data_min_dim <- min(nrow(x = data), ncol(x = data) / k)
# Check if --to-n-pc is not violating requirement and adapt the parameters
if(args$`to-n-pc` > data_min_dim) {
	print("The --to-n-pc parameter is greater than the smallest dimension of the data given the k-fold setting.")
	if(k < 10) {
		stop("Please adapt the following parameters: --k-fold, --from-n-pc, --to-n-pc, --by-n-pc.")
	}
	k <- 5
	print(paste0("The k-fold parameter is decreased to ", k, "."))
	from_npcs <- 2
	to_npcs <- data_min_dim
	by_npcs <- 1
	repeat {
		pc <- seq(from = from_npcs, to = to_npcs, by = by_npcs)
		if(length(x = pc) < 50) {
			break
		} else {
			by_npcs <- by_npcs + 1
		}
	}
	print(paste0("Setting --from-n-pc parameter to ", from_npcs, " instead of ", args$`from-n-pc`, "."))
	print(paste0("Setting --to-n-pc parameter to ", to_npcs, " instead of ", args$`to-n-pc`, "."))
	print(paste0("Setting --by-n-pc parameter to ", by_npcs, " instead of ", args$`by-n-pc`, "."))
}

# Run
out <- RunPCACV(
	data = data,
	k = k,
	from = from_npcs,
	to = to_npcs,
	by = by_npcs,
	maxit = args$`max-iters`,
	seed = args$`seed`,
	n.cores = args$`n-cores`,
	verbose = args$`verbose`,
	default.svd = args$`default-svd`
)


pcs_pred <- NULL

# Fit to get optimal number of PCs
if(length(x = out$error) >= 4) {
	print("Find the number of optimal number of PCs by fitting spline...")
	out_fit <- smooth.spline(x = out$PC, out$error)
	pcs_pred <- predict(object = out_fit, x = 1:to_npcs)
	optimum_npcs <- pcs_pred$x[which.min(pcs_pred$y)]
} else {
	print("Find the number of optimal number of PCs without fitting a spline (too few values) ...")
	optimum_npcs <- out$PC[which.min(out$error)]
}

print(paste0("Optimal number of PCs: ", optimum_npcs))

out_npcs <- NULL

if(optimum_npcs == 1) {
	if(args$`n-pc-fallback` < 0) {
		stop(paste0("Invalid value for --n-pc-fallback parameter: value should be > 0."))
	}
	if(args$`n-pc-fallback` == 0) {
 		write.table(
			x = out,
			file = paste0(args$`output-prefix`, ".PRESS_ERRORS.tsv"),
			quote = F,
			sep = "\t",
			row.names = F,
			col.names = T
		)
		stop(paste0("Could not find an optimal number of PCs. You can either set --n-pc-fallback parameter (> 0) in order to return a minimum number of PCs or adapt the following parameters: --k-fold, --from-n-pc, --to-n-pc, --by-n-pc."))
	}
	msg <- paste0("No optimal number of PCs found. The number of PCs returned is defined by --n-pc-fallback: ", args$`n-pc-fallback`)
	warning(msg)
	file_conn <- file(paste0(args$`output-prefix`, ".WARNING.txt"))
	writeLines(c(msg), file_conn)
	close(file_conn)
	out_npcs <- args$`n-pc-fallback`
} else {
	out_npcs <- optimum_npcs
}

# Save the results

## Save parameters
yaml::write_yaml(x = yaml::as.yaml(x = args), file = ".PARAMETERS.yaml")

## PRESS errors
pcs_data <- NULL

if(!is.null(x = pcs_pred)) {
	pcs_data <- as.data.frame(x = pcs_pred)
} else {
	pcs_data <- data.frame(x = out$PC, y = out$error)
}
colnames(x = pcs_data) <- c("PC", "error")
write.table(
	x = pcs_data,
	file = paste0(args$`output-prefix`, ".PRESS_ERRORS.tsv"),
	quote = F,
	sep = "\t",
	row.names = F,
	col.names = T
)

## PRESS errors plot

pdf(file = paste0(args$`output-prefix`, ".PRESS_ERROR_PLOT.pdf"))
plot(pcs_data)
lines(pcs_data, col = "blue")
dev.off()

## Optimal number of principal components

write.table(
	out_npcs,
	paste0(args$`output-prefix`, ".OPTIMAL_NPCS.txt"),
	append = FALSE,
	row.names = FALSE,
	col.names = FALSE
)

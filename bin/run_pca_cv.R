#!/usr/bin/env Rscript

print("##################################################")
print("# Calculating Optimal Number of PCs using PCA CV #")
print("##################################################")

# Loading dependencies scripts

# matrix <- GetDataFromS4ObjectByAccessPath(object = combined, path = "@assays$integrated@scale.data")
# e.g.: Rscript /ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/scripts/src_dwmax/ngs_digest/bin/scrna_seq/dim_reduction_pca_cv.R -i 10x_REW_Seurat_v3_Integrate_MaxFeatures.Rds -a "@assays\$integrated@scale.data"

library("optparse")
parser <- OptionParser(prog = "dim_reduction_pca_cv.R", description = "Perform PCA through cross-validation to find the optimal number of principal components.")
parser <- add_option(parser, c("-i", "--input-file"), action="store", default = NULL, help="Input file [default]")
parser <- add_option(parser, c("-a", "--accessor"), action = "store", default = NULL, help="Character defining the path to access the data matrix in S4 object.")
parser <- add_option(parser, c("-k", "--k-fold"), action="store", default = 10, help="Set value k for the k-fold cross-validation. [default %default]")
parser <- add_option(parser, c("-f", "--from-n-pc"), action="store", default = 2, help="Number of principal components to start with. [default %default]")
parser <- add_option(parser, c("-t", "--to-n-pc"), action="store", default = 150, help="Number of principal components to compute [default %default]")
parser <- add_option(parser, c("-b", "--by-n-pc"), action="store", default = 5, help="Number of principal components to increment by. [default %default]")
parser <- add_option(parser, c("-m", "--max-iters"), action="store", default = 200, help="Maximum number of iteration. [default %default]")
parser <- add_option(parser, c("-s", "--seed"), action="store", default = 617, help="Seed. [default %default]")
parser <- add_option(parser, c("-c", "--n-cores"), action="store", default = 1, help="Number of cores to use. [default %default]")
parser <- add_option(parser, c("-d", "--default-svd"), action="store", default = FALSE, help="Perform SVD instead of approximate PCA. [default %default]")
parser <- add_option(parser, c("-v", "--verbose"), action="store", default = TRUE, help="Print some informational messages. [default %default]")
parser <- add_option(parser, c("-o", "--output-prefix"), action="store", default = "dim_reduction_pca_cv_out_table.tsv", help="Prefix path to save output files (table of the PCA CV, PCA PRESS plot, optimal number of PCs). [default %default]")

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
  seed = 617,
  n.cores = 1,
  verbose = T,
  default.svd = F) {
  set.seed(seed)
  library(foreach)
  library(doMC)
  registerDoMC(cores = n.cores) 
  library(irlba)
  if(default.svd) {
    print("Default SVD is using for estimating the number PCs.") 
  }
  
  pc <- seq(from = from, to = to, by = by)
  # Init the error matrices
  error <- matrix(0, nrow = length(c(1:k)), ncol = length(x = pc))

  print(paste0(k,"-fold paritioning..."))
  # K-Fold Partitioning
  data_kfold <- dismo::kfold(x = Matrix::t(x = data), k=k)
  
  # Reference: https://stats.stackexchange.com/questions/93845/how-to-perform-leave-one-out-cross-validation-for-pca-to-determine-the-number-of
  for(i in c(1:k)) {
    print(paste0("k:",i))
    data_train <- Matrix::t(x = data[, data_kfold!=i])
    data_test <- Matrix::t(x = data[, data_kfold==i])
    print("Running SVD...")
    if(!default.svd) {
      pca_results <- irlba::irlba(A = data_train, nv = to, maxit = maxit, verbose = verbose)
    } else {
      pca_results <- base::svd(x = data_train, nv = to)
    }
    gl <- pca_results$v
    res <- foreach(j=1:length(x = pc), .combine='rbind') %dopar% {
      print(paste0("Ndims:", pc[j]))
      P <- gl[,1:pc[j]]%*%Matrix::t(gl[,1:pc[j]])
      err <- data_test %*% (diag(x = dim(x = P)[1]) - P + diag(x = diag(x = P)))
      data.frame("i"=i, "j"=j, err=sum(err^2))
    }
    apply(X = res, MARGIN = 1, FUN = function(r) {
      error[r[["i"]],r[["j"]]] <<- r[["err"]]
    })
  }
  errors <- colSums(x = error)
  return (data.frame("PC"=pc,"error"=log(x = errors)))
}

cat("Parameters: \n")
print(args)

input_ext <- tools::file_ext(args$`input-file`)

if(input_ext == "loom") {
  print("Input type: Loom...")
  # Get the data matrix from loom 
  warnings("Assumption: Expecting a scaled data matrix in the main layer of the loom file!")
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
} else if(input_ext == "h5ad") {
  seurat <- Seurat::ReadH5AD(file = args$`input-file`)
  if(!methods::.hasSlot(object = seurat@assays$RNA, name = "scale.data") || is.null(x = slot(object = seurat@assays$RNA, name = "scale.data")))
    stop("Expects a feature-scaled data matrix but it does not exist.")
  if(is.null(args$`accessor`))
    stop("The --accessor argument need to be set.")
  data <- GetDataFromS4ObjectByAccessPath(
    object = seurat
    , path = args$`accessor`
  )
} else {
  stop(paste0("Unrecognized input file format: ", input_ext, "."))
}

print(paste0("Data matrix has ", dim(x = data)[1], " rows, ", dim(x = data)[2], " columns."))

# Run
out <- RunPCACV(
  data = data,
  k = args$`k-fold`,
  from = args$`from-n-pc`,
  to = args$`to-n-pc`,
  by = args$`by-n-pc`,
  maxit = args$`max-iters`,
  seed = args$`seed`,
  n.cores = args$`n-cores`,
  verbose = args$`verbose`,
  default.svd = args$`default-svd`
)

# Fit to get optimal number of PCs
out_fit <- smooth.spline(x = out$PC, out$error)
pcs_pred <- predict(object = out_fit, x = 1:args$`to-n-pc`)

optimum_npcs <- pcs_pred$x[which.min(pcs_pred$y)]
print(paste0("Optimal number of PCs: ", optimum_npcs))

# Save the results

## Save parameters
yaml::write_yaml(x = yaml::as.yaml(x = args), file = ".PARAMETERS.yaml")

## PRESS errors
pcs_pred <- as.data.frame(x = pcs_pred)
colnames(x = pcs_pred) <- c("PC", "error")
write.table(
  x = pcs_pred,
  file = paste0(args$`output-prefix`, ".PRESS_ERRORS.tsv"),
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)

## PRESS errors plot

pdf(file = paste0(args$`output-prefix`, ".PRESS_ERROR_PLOT.pdf"))
plot(pcs_pred)
lines(pcs_pred, col = "blue")
dev.off()

## Optimal number of principal components

write.table(
  optimum_npcs,
  paste0(args$`output-prefix`, ".OPTIMAL_NPCS.txt"),
  append = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

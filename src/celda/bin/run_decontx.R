#!/usr/bin/env Rscript

print("##################################################")
print("# Contamination estimation with decontX          #")
print("##################################################")

library("argparse")

parser <- ArgumentParser(description='Identifies contamination from factors such as ambient RNA in single cell genomic datasets.')
parser$add_argument(
    'input',
    metavar='INPUT',
    type="character",
    help='A Rds file containing a SingleCellExperiment object.'
)
parser$add_argument(
    '--output-prefix',
    type="character",
    dest='output_prefix',
    default = "foo",
  	help="Prefix path to save output files (plots and tables). [default %default]"
)
parser$add_argument(
    '--sample-id',
    type="character",
    dest='sample_id',
    help='Sample ID of the given input file.'
)
parser$add_argument(
    '--seed',
    type="character",
    dest='seed',
    default=12345,
    help='Passed to with_seed. For reproducibility, a default value of 12345 is used. If NULL, no calls to with_seed are made.'
)
parser$add_argument(
    '--num-mads-threshold',
    type="integer",
    dest='num_mads_thresholds',
    action="append",
	default=NULL,
    help='Number of median absolution deviation (MAD) to decide for outlier assignment.'
)
parser$add_argument(
    '--custom-threshold',
    type="double",
    dest='custom_thresholds',
    action="append",
    default=NULL,
    help='Threshold to decide for outlier assignment. Should be float between 0 and 1.'
)
parser$add_argument(
	'--round-to-int',
	type="logical",
	dest='round_to_int',
	action="store",
	default=FALSE,
	help='Round the corrected matrix by DecontX.'
)
parser$add_argument(
	'--filter-empty-cells',
	type="logical",
	dest='filter_empty_cells',
	action="store",
	default=FALSE,
	help='Filter empty cells after DecontX decntamination'
)

args <- parser$parse_args()
print(args)

isTrue <- function(x) {
	x <- tolower(x = x)
	return (
		x == FALSE | x == F | x == 't' | x == 'true' | x == 'y' | x == 'yes'
	)
}

# IO
print("Reading SingleCellExperiment Rds...")
sce <- tryCatch({
	readRDS(file = args$input)
}, warning = function(w) {
}, error = function(e) {
	stop(paste0("Cannot read the given Rds file: ", args$input))
}, finally = {})

# Load libraries
library(ggplot2)

# Run
print("Running DecontX...")
sce <- celda::decontX(
	x = sce,
	seed = args$seed
)

# Post-process
print("Generating the different outlier masks...")
col_data <- SummarizedExperiment::colData(x = sce)

thresholds <- list()

# Outler assignment based on DecontX contamination and nMAD deviation threshold
if(!is.null(x = args$num_mads_thresholds)) {
	for(num_nmad_threshold in args$num_mads_thresholds) {
		decontX_contamination <- scater::isOutlier(
			metric = col_data$decontX_contamination,
			log=FALSE,
			type="higher",
			nmads = num_nmad_threshold
		)
		decontX_contamination_threshold <- attr(
			x = decontX_contamination,
			which = "thresholds"
		)[2]
		outlier_threshold_name <- paste0("scater_isOutlier_", num_nmad_threshold, "MAD")
		thresholds[[length(thresholds)+1]] <- list("name" = outlier_threshold_name, "value" = decontX_contamination_threshold)
		col_data[[paste0("celda_decontx__", outlier_threshold_name, "_predicted_outliers")]] <- decontX_contamination
	}
}

# Outler assignment based on DecontX contamination and DoubleMAD threshold
# Source: https://aakinshin.net/posts/harrell-davis-double-mad-outlier-detector/
getOutliersDoubleMADThresholds <- function(x) {
	hdmedian <- function(u) as.numeric(Hmisc::hdquantile(u, 0.5))
	x <- x[!is.na(x)]
	m <- hdmedian(x)
	deviations <- abs(x - m)
	lowerMAD <- 1.4826 * hdmedian(deviations[x <= m])
	upperMAD <- 1.4826 * hdmedian(deviations[x >= m])
	lowerMAD_threshold <- m - 3 * lowerMAD
	upperMAD_threshold <- m + 3 * upperMAD
	return (c(lowerMAD_threshold, upperMAD_threshold))
}
outlier_doublemad_thresholds <- getOutliersDoubleMADThresholds(
    x = col_data$decontX_contamination
)
outlier_threshold_name <- "doublemad"
thresholds[[length(thresholds)+1]] <- list("name" = outlier_threshold_name, "value" = outlier_doublemad_thresholds[2])
col_data[[paste0("celda_decontx__", outlier_threshold_name, "_predicted_outliers")]] <- col_data$decontX_contamination < outlier_doublemad_thresholds[1] | col_data$decontX_contamination > outlier_doublemad_thresholds[2]

# Outler assignment based on DecontX contamination and custom threshold(s)
if(!is.null(x = args$custom_thresholds)) {
    for(ct in args$custom_thresholds) {
        outlier_threshold_name <- paste0("custom_gt_", ct)
		thresholds[[length(thresholds)+1]] <- list("name" = outlier_threshold_name, "value" = ct)
        col_data[[paste0("celda_decontx__", outlier_threshold_name, "_predicted_outliers")]] <- col_data$decontX_contamination > ct
    }
}

# Output

## Plot DecontX clusters on UMAP
print("Plotting DecontX clusters on UMAP...")
umap <- SingleCellExperiment::reducedDim(x = sce, type = "decontX_UMAP")
celda::plotDimReduceCluster(umap[,1], umap[,2], cluster = sce$decontX_clusters)
ggsave(
    filename=paste0(args$`output_prefix`, ".UMAP_Clusters.pdf")
)

## Plot DecontX contamination score on UMAP
print("Plotting DecontX contamination score on UMAP...")
celda::plotDecontXContamination(sce)
ggsave(
    filename=paste0(args$`output_prefix`, ".UMAP_Contamination_Score.pdf")
)

## Plot DecontX contamination score density with outliers
print("Plotting DecontX contamination score density with outliers...")
for(threshold_idx in seq_along(x = thresholds)) {
	threshold <- thresholds[[threshold_idx]]
	threshold_accessor_name <- paste0("celda_decontx__", threshold$name, "_predicted_outliers")
	num_outliers <- sum(col_data[[threshold_accessor_name]])
	p <- ggplot(as.data.frame(x = col_data, stringsAsFactors = FALSE), aes(x=decontX_contamination)) + 
		geom_density(color="darkblue", fill="lightblue") +
		ggtitle(
			label = paste(args$sample_id, "-", "DecontX Contamination Density Plot"),
			subtitle = paste0("Num. Outlier Cells: ", round(x = num_outliers/nrow(col_data)*100, digits = 2), "%", " (", num_outliers, " / ", nrow(col_data) ,")")
		)
	# Check if there is at least 1 outlier otherwise
	# Error: `data` must be uniquely named but has duplicate columns
	if(num_outliers > 0) {
		d <- ggplot_build(plot = p)$data[[1]]
		p <- p + 
			geom_area(data = subset(d, x > threshold$value), aes(x=x,y=y), fill="orange")
	}
	print(p)
	ggsave(
		filename=paste0(args$`output_prefix`, paste0(".Contamination_Score_Density_with_", threshold$name,".pdf"))
	)
}

## Outlier tresholds
print("Saving outlier thresholds...")
thresholds  <-  as.data.frame(
	matrix(
		unlist(x = thresholds),
		nrow=length(x = unlist(x = thresholds[1]))
	)
)
row.names(x = thresholds)<- c("name","value")
write.table(
	thresholds,
	paste0(args$`output_prefix`, ".Contamination_Outlier_Thresholds.tsv"),
	append = FALSE,
	row.names = TRUE,
	col.names = FALSE,
	quote = FALSE,
	sep = "\t"
)


## Outlier table
print("Saving outlier table...")
col_data$Sample <- NULL
# Renames decontx_ colnames to decontx_celda__ to harmony augmented data generated by this tool
colnames(x = col_data) <- gsub(
	pattern = "decontX_",
	replacement = "celda_decontx__",
	x = colnames(x = col_data)
)
write.table(
	cbind(
		data.frame("index"=row.names(x = col_data), stringsAsFactors=FALSE),
		col_data
	),
	paste0(args$`output_prefix`, ".Contamination_Outlier_Table.tsv"),
	append = FALSE,
	row.names = FALSE,
	col.names = TRUE,
	quote = FALSE,
	sep = "\t"
)

## SingleCellExperiment with DecontX results
print("Saving SingleCellExperiment with DecontX results as Rds...")
if(isTrue(x = args$round_to_int)) {
	print("Rounding the DecontX count matrix...")
	mat <- celda::decontXcounts(object = sce)
	mat <- round(
		x = mat,
		digits = 0
	)
	celda::decontXcounts(object = sce) <- mat
	rm(mat)
	gc()
}

# Check if cells has no counts
mat <- celda::decontXcounts(object = sce)
if(any(Matrix::colSums(x = mat) == 0)) {
	empty_cells <- Matrix::colSums(x = mat) == 0
	num_empty_cells <- sum(empty_cells)
	if(!args$filter_empty_cells) {
		stop(paste0("VSN ERROR: Some cells (", num_empty_cells,") have no counts after DecontX decontamination. You can overcome this by specifying --filter-empty-cells (filterEmptyCells = true) true."))
	}
	print(paste0("VSN MSG: Some cells (", num_empty_cells,") have no counts after DecontX decontamination. Filtering those empty cells out..."))
	print(paste0("VSN MSG: Number of cells before the empty cells filter: ", ncol(x = sce)))
	col_data <- col_data[!empty_cells, ]
	sce <- sce[, !empty_cells]
	print(paste0("VSN MSG: Number of cells after the empty cells filter: ", ncol(x = sce)))
}

SummarizedExperiment::colData(x = sce) <- col_data
saveRDS(
	object = sce,
	paste0(args$`output_prefix`, ".Rds"),
	compress = TRUE
)

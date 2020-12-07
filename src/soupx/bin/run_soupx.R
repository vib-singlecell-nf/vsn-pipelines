#!/usr/bin/env Rscript

print("################################################################################################################")
print("# Estimation and removal of cell free mRNA contamination in droplet based single cell RNA-seq data with SoupX  #")
print("################################################################################################################")

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
	'--round-to-int',
	type="logical",
	dest='round_to_int',
	action="store",
	default=FALSE,
	help='Round the corrected matrix by DecontX.'
)

args <- parser$parse_args()
print(args)

set.seed(seed = args$seed)

# IO
print("VSN MSG: Reading SingleCellExperiment Rds...")
sce <- tryCatch({
	readRDS(file = args$input)
}, warning = function(w) {
}, error = function(e) {
	stop(paste0("VSN ERROR: Cannot read the given Rds file: ", args$input))
}, finally = {})

# Run

Is10xOuts <- function(path) {
    return (basename(path = path) == "outs")
}

IsTrue <- function(x) {
	x <- tolower(x = x)
	return (
		x == FALSE | x == F | x == 't' | x == 'true' | x == 'y' | x == 'yes'
	)
}

data_dir <- args$input
if(!Is10xOuts(path = data_dir)) {
    parent <- dirname(path = data_dir)
    if(!Is10xOuts(path = parent)) {
        stop("VSN ERROR: Cannot infer top level cellranger output directory (the directory that contains the raw_gene_bc_matrices folder).")
    }
    data_dir <- parent
}

soupchannel = SoupX::load10X(dataDir = data_dir)
soupchannel = SoupX::autoEstCont(sc = soupchannel)
if(IsTrue(x = args$round_to_int)) {
	print("VSN MSG: The SoupX corrected count matrix will be rounded...")
}
soupx_counts = SoupX::adjustCounts(
    sc = soupchannel,
    roundToInt = args$round_to_int
)

## Contamination Fraction Plot
source(file.path(Sys.getenv("NXF_BIN_DIR"), "plot.R"))

pdf(file = paste0(args$output_prefix, ".Contamination_Fraction_Density_Plot.pdf"))
PlotContaminationFraction(object = soupchannel, title = paste0(args$sample_id, " - Contamination Fraction Density Plot"))
dev.off()

## SingleCellExperiment with SoupX results
print("VSN MSG: Saving SingleCellExperiment with SoupX results as Rds...")
# colData
col_data <- soupchannel$metaData
tsne <- as.matrix(x = col_data[, c("tSNE1", "tSNE2")])
col_data <- col_data[, !(colnames(x = col_data) %in% c("tSNE1", "tSNE2"))]
## Rename columns
colnames(x = col_data) <- paste0("soupx__", colnames(x = col_data))
# rowData
row_data <- soupchannel$soupProfile
# Rename columns
colnames(x = row_data) <- paste0("soupx__", colnames(x = row_data))
sce <- SingleCellExperiment::SingleCellExperiment(
    assays=list(soupXcounts=soupx_counts),
    reducedDims=S4Vectors::SimpleList(tSNE=tsne),
    colData=col_data,
    rowData=row_data
)
saveRDS(
	object = sce,
	paste0(args$`output_prefix`, ".Rds"),
	compress = TRUE
)

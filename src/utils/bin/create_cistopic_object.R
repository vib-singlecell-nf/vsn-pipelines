#!/usr/bin/env Rscript

print("##################################################")
print("# cisTopic #")
print("##################################################")

# Loading dependencies scripts

library("optparse")
parser <- OptionParser(
  prog = "create_cistopic_object.R",
  description = "Create cisTopic object from 10x Cell Ranger MEX output"
)
parser <- add_option(
  parser,
  c("-i", "--tenx_path"),
  action = "store",
  default = NULL,
  help = "Path to Cell Ranger 10x output containing filtered_peak_bc_matrix/ directory"
)
parser <- add_option(
  parser,
  c("-m", "--metrics_fname"),
  action = "store",
  default = "singlecell.csv",
  help = "Filename of Cell Ranger 10x output per barcode metrics"
)
parser <- add_option(
  parser,
  c("-s", "--sampleId"),
  action = "store",
  default = "",
  help = "sample ID"
)
parser <- add_option(
  parser,
  c("-o", "--output"),
  action = "store",
  default = NULL,
  help = "Output file, rds format"
)

args <- parse_args(parser)

cat("Parameters: \n")
print(args)

################################################################################

suppressWarnings(library(cisTopic))

data_folder = file.path(args$tenx_path, 'filtered_peak_bc_matrix')
metrics = file.path(args$tenx_path, args$metrics_fname)

cisTopicObject <- createcisTopicObjectFrom10Xmatrix(data_folder, metrics, project.name='VSN-ATAC')

saveRDS(cisTopicObject,file=args$output)


#!/usr/bin/env Rscript

print("##################################")
print("# ArchR: Fragments -> Arrow file #")
print("##################################")

################################################################################

library("optparse")
parser <- OptionParser(
  prog = "createArrow_unfiltered.R",
  description = "Load a fragments file, create an Arrow file without discarding non-cell barcodes."
)
parser <- add_option(
  parser,
  c("-f", "--fragments_file"),
  action = "store",
  default = NULL,
  help = "Input file [default]"
)
parser <- add_option(
  parser,
  c("--output_directory"),
  action = "store",
  default = NULL,
  help = "output directory prefix for plot file"
)
parser <- add_option(
  parser,
  c("--sample_name"),
  action = "store",
  default = "sample_id",
  help = "Sample ID to name the output file [%default]"
)
parser <- add_option(
  parser,
  c("--min_frags"),
  action = "store",
  default = 10,
  help = "Minimum number of fragments to keep per barcode [%default]"
)
parser <- add_option(
  parser,
  c("--filter_tss"),
  action = "store",
  default = 4,
  help = "Filter on TSS Enrichment [%default]"
)
parser <- add_option(
  parser,
  c("--filter_frags"),
  action = "store",
  default = 1000,
  help = "Filter on fragments/cell [%default]"
)
parser <- add_option(
  parser,
  c("-g", "--genome"),
  action = "store",
  default = "hg38",
  help = "Genome to use, choose from (hg19, hg38, mm9, and mm10) [%default]"
)

parser <- add_option(
  parser,
  c("-t", "--threads"),
  action = "store",
  default = 1,
  help="Number of threads to use. [default %default]"
)
parser <- add_option(
  parser,
  c("-s", "--seed"),
  action = "store",
  default = 617,
  help="Seed. [default %default]"
)

args <- parse_args(parser)

cat("Parameters: \n")
print(args)


################################################################################

suppressMessages(library(ArchR, quietly=TRUE))

set.seed(args$seed)
addArchRThreads(threads=args$threads)
addArchRGenome(args$genome)

ArrowFiles = createArrowFiles(
    inputFiles=args$fragments_file,
    sampleNames=args$sample_name,
    #outputNames=paste0(args$output_directory,"/",args$sample_name),
    removeFilteredCells=FALSE,
    minFrags=args$min_frags,
    maxFrags=1000000,
    #QCDir=paste0(args$output_directory,"/QualityControl"),
    nChunk=1,
    filterTSS=args$filter_tss,
    filterFrags=args$filter_frags, 
    force=TRUE,
    addTileMat=FALSE,
    addGeneScoreMat=FALSE
)



#!/usr/bin/env Rscript 

args <- base::commandArgs(trailingOnly=TRUE)

nomenclature <- read.table(
  file = args[1], 
  header = T, 
  sep = "\t", 
  stringsAsFactors = F
)

library(flybaseR)
nomenclature_converted <- id.converter(x = nomenclature[[args[2]]], symbols = TRUE, convert.into = "g")

if(length(x = nomenclature[[args[2]]]) != length(nomenclature_converted)) {
  stop("Problem occurred during FBgn to gene symbol conversion: dimensions are not the same.")
}

nomenclature$gene_symbol <- nomenclature_converted

write.table(
  x = nomenclature,
  file = args[3],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

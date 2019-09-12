#!/usr/bin/env python
import os
from optparse import OptionParser
import scanpy as sc
import anndata as ad

parser = OptionParser(usage="usage: %prog [options] h5ad_file_path",
                      version="%prog 1.0")
parser.add_option("-c", "--min-number-cells",
                  type=int,
                  action="store",
                  dest="min_number_cells",
                  default=-1,
                  help="Filter out genes that are detected in less than the minimum number of cells.")
(options, args) = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args[0]
FILE_PATH_OUT_BASENAME = os.path.splitext(args[1])[0]

# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(FILE_PATH_IN))

#
# Filter on min number of cells
#

if options.min_number_cells > -1:
    sc.pp.filter_genes(adata, min_cells=options.min_number_cells)

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

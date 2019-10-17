#!/usr/bin/env python
import os
from optparse import OptionParser

import numpy as np

import scanpy as sc

parser = OptionParser(
    usage="usage: %prog [options] h5ad_file_path",
    version="%prog 1.0"
)
parser.add_option(
    "-c", "--min-n-counts",
    type=int,
    action="store",
    dest="min_n_counts",
    default=-1,
    help="Filter out cells with less than the minimum n of counts."
)
parser.add_option(
    "-C", "--max-n-counts",
    type=int,
    action="store",
    dest="max_n_counts",
    default=-1,
    help="Filter out cells with more than the maximum n of counts."
)
parser.add_option(
    "-g", "--min-n-genes",
    type=int,
    action="store",
    dest="min_n_genes",
    default=-1,
    help="Filter out cells with less than the minimum n of genes expressed."
)
parser.add_option(
    "-G", "--max-n-genes",
    type=int,
    action="store",
    dest="max_n_genes",
    default=-1,
    help="Filter out cells with more than the maximum n of genes expressed."
)
parser.add_option(
    "-M", "--max-percent-mito",
    type=float,
    action="store",  # optional because action defaults to "store"
    dest="max_percent_mito",
    default=-1,
    help="Filter out cells with more than the maximum percentage of mitochondrial genes expressed."
)
parser.add_option(
    "--min-number-cells",
    type=int,
    action="store",
    dest="min_number_cells",
    default=-1,
    help="Filter out genes that are detected in less than the minimum number of cells."
)
#
parser.add_option(
    "-v", "--verbose",
    action="store_false",  # optional because action defaults to "store"
    dest="verbose",
    default=False,
    help="Show messages."
)

(options, args) = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args[0]
FILE_PATH_OUT_BASENAME = os.path.splitext(args[1])[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))


#
# Compute number of counts
#

def compute_n_counts(adata):
    adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))

compute_n_counts(adata=adata)

#
# Compute number of genes
#

def compute_n_genes(adata):
    # simply compute the number of genes per cell (computes 'n_genes' column)
    sc.pp.filter_cells(adata, min_genes=0)

compute_n_genes(adata=adata)

#
# Compute fraction of mitochondrial genes
#

def compute_percent_mito(adata):
    # mito and genes/counts cuts
    mito_genes = adata.var_names.str.contains('^MT-|^mt[-:]')
    # for each cell compute fraction of counts in mito genes vs. all genes
    adata.obs['percent_mito'] = np.ravel(np.sum(adata[:, mito_genes].X, axis=1)) / np.ravel(np.sum(adata.X, axis=1))
    return adata

adata = compute_percent_mito(adata=adata)

#
# Write filtering value into adata.uns
#
adata.uns['sc']={
        'scanpy': {
            'filter': {
                'cellFilterMinNGenes': options.min_n_genes,
                'cellFilterMaxNGenes': options.max_n_genes,
                'cellFilterMaxPercentMito': options.max_percent_mito,
                'geneFilterMinNCells': options.min_number_cells,
            }

        }
    }

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))


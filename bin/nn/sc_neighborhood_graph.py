#!/usr/bin/env python3

import argparse
import os
import warnings
import scanpy as sc

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input h5ad file.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output h5ad file.'
)

parser.add_argument(
    "-n", "--n-neighbors",
    type=int,
    action="store",
    dest="n_neighbors",
    default=15,
    help="[UMAP], The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation."
)

parser.add_argument(
    "-p", "--n-pcs",
    type=int,
    action="store",
    dest="n_pcs",
    default=None,
    help="Use this many PCs. If n_pcs==0 use .X if use_rep is None."
)

parser.add_argument(
    "-s", "--seed",
    type=int,
    action="store",
    dest="seed",
    default=0,
    help="Use this integer seed for reproducibility."
)


args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

#
# Calculate the neighboorhood graph
#

sc.pp.neighbors(
    adata=adata,
    n_neighbors=args.n_neighbors,
    n_pcs=args.n_pcs,
    random_state=args.seed
)

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

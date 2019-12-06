#!/usr/bin/env python3

import argparse
import os
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
    "-x", "--method",
    type=str,
    action="store",
    dest="method",
    default="Louvain",
    help="Cluster cells using the Louvain algorithm [Blondel et al. (2008)] in the implementation of [Traag (2017)]."
)

parser.add_argument(
    "-r", "--resolution",
    type=float,
    action="store",
    dest="resolution",
    default=1.0,
    help="[louvain], For the default flavor ('vtraag'), you can provide a resolution (higher resolution means finding"
         " more and smaller clusters) (Default: 1.0)."
)

parser.add_argument(
    "-n", "--n-neighbors",
    type=int,
    action="store",
    dest="n_neighbors",
    default=15,
    help="[Louvain], The size of local neighborhood (in terms of number of neighboring data points) used for manifold"
         " approximation."
)

parser.add_argument(
    "-p", "--n-pcs",
    type=int,
    action="store",
    dest="n_pcs",
    default=30,
    help="[Louvain], Use this many PCs."
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
    raise Exception("Can only handle .h5ad files.")

#
# Transform the distribution of the data
#

if args.method == "Louvain":
    # Run Louvain clustering
    if "neighbors" not in adata.uns.keys():
        raise Exception(
            "The neighborhood graph of observations has not been computed. Please do so before running {} clustering".format(
                args.method)
        )
        # sc.pp.neighbors(adata, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)
    sc.tl.louvain(adata, resolution=args.resolution)
else:
    raise Exception("The dimensionality reduction method {} does not exist.".format(args.method))

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

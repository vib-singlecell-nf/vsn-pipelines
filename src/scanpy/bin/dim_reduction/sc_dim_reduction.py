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
    "-x", "--method",
    type=str,
    action="store",
    dest="method",
    default="PCA",
    help="Reduce the dimensionality of the data. Choose one of : PCA, UMAP, t-SNE"
)

parser.add_argument(
    "-c", "--n-comps",
    type=int,
    action="store",
    dest="n_comps",
    default=50,
    help="[PCA], Number of principal components to compute."
)

parser.add_argument(
    "-s", "--svd-solver",
    type=str,
    action="store",
    dest="svd_solver",
    default="arpack",
    help="[PCA], SVD solver to use. Choose one of : arpack (Default), randomized, auto."
)

parser.add_argument(
    "-n", "--n-neighbors",
    type=int,
    action="store",
    dest="n_neighbors",
    default=15,
    help="[Louvain], The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation."
)

parser.add_argument(
    "-p", "--n-pcs",
    type=int,
    action="store",
    dest="n_pcs",
    default=30,
    help="[Louvain], Use this many PCs."
)

parser.add_argument(
    "-j", "--n-jobs",
    type=int,
    action="store",
    dest="n_jobs",
    default=1,
    help="The number of jobs. When set to None, automatically uses the number of cores."
)

parser.add_argument(
    "-f", "--use-fast-tsne",
    action="store_true",
    dest="use_fast_tsne",
    default=False,
    help="Use the MulticoreTSNE package by D. Ulyanov if it is installed."
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

if args.method == "PCA":
    # Run PCA
    sc.tl.pca(
        data=adata,
        n_comps=min(adata.shape[0], args.n_comps),
        svd_solver=args.svd_solver
    )
elif args.method == "UMAP":
    # Run UMAP
    # Notes:
    # - /!\ BBKNN is slotting into the sc.pp.neighbors() => sc.pp.neighbors() should not be run afterwards otherwise results will be overwritten
    if "neighbors" not in adata.uns.keys():
        warnings.warn("The neighborhood graph of observations has not been computed. Computing...")
        sc.pp.neighbors(
            adata=adata,
            n_neighbors=args.n_neighbors,
            n_pcs=args.n_pcs
        )
    sc.tl.umap(adata)
elif args.method == "t-SNE":
    # Run t-SNE
    sc.tl.tsne(adata=adata, n_jobs=args.n_jobs, use_fast_tsne=args.use_fast_tsne)
else:
    raise Exception("The dimensionality reduction method {} does not exist.".format(args.method))

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

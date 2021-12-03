#!/usr/bin/env python3

import argparse
import os
import warnings
import scanpy as sc


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


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
    default="pca",
    help="Reduce the dimensionality of the data. Choose one of : PCA, UMAP, t-SNE"
)

parser.add_argument(
    "-p", "--perplexity",
    type=int,
    action="store",
    dest="perplexity",
    default=30,
    help="[t-SNE], The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms."
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
    "-v", "--svd-solver",
    type=str,
    action="store",
    dest="svd_solver",
    default="arpack",
    help="[PCA], SVD solver to use. Choose one of : arpack (Default), randomized, auto."
)

parser.add_argument(
    "-n", "--n-pcs",
    type=int,
    action="store",
    dest="n_pcs",
    default=None,
    help="[UMAP, t-SNE], Use this many PCs. If n_pcs==0 use .X if use_rep is None."
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
    type=str2bool,
    action="store",
    dest="use_fast_tsne",
    default=True,
    help="Use the MulticoreTSNE package by D. Ulyanov if it is installed."
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
# Transform the distribution of the data
#

if args.method.lower() == "pca":
    # Run PCA
    sc.tl.pca(
        data=adata,
        n_comps=min(adata.shape[0], args.n_comps),
        svd_solver=args.svd_solver,
        random_state=args.seed
    )
elif args.method.lower() == "umap":
    # Run UMAP
    # Notes:
    # - /!\ BBKNN is slotting into the sc.pp.neighbors() => sc.pp.neighbors() should not be run afterwards otherwise results will be overwritten
    if "neighbors" not in adata.uns.keys():
        raise Exception("VSN ERROR: The neighborhood graph of observations has not been computed. Computing...")
    sc.tl.umap(
        adata=adata,
        random_state=args.seed
    )
elif args.method.lower() == "tsne":
    # Run t-SNE
    # Source: https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.tsne.html
    # If n_pcs is None and X_pca has been computed, X_pca with max computed pcs will be used
    # ||
    # n_pcs : int, None (default: None)
    #     Use this many PCs. If n_pcs==0 use .X if use_rep is None.
    # use_rep : str, None (default: None)
    #     Use the indicated representation. 'X' or any key for .obsm is valid. If None, the representation is chosen automatically: For .n_vars < 50, .X is used, otherwise ‘X_pca’ is used. If ‘X_pca’ is not present, it’s computed with default parameters.
    sc.tl.tsne(
        adata=adata,
        perplexity=args.perplexity,
        n_pcs=args.n_pcs,
        random_state=args.seed,
        use_fast_tsne=args.use_fast_tsne,
        n_jobs=args.n_jobs
    )
else:
    raise Exception("VSN ERROR: The dimensionality reduction method {} does not exist.".format(args.method))

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

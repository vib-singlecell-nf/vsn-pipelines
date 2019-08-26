#!/usr/bin/env python
import os
from optparse import OptionParser
import scanpy as sc
import anndata as ad

parser = OptionParser(usage="usage: %prog [options] h5ad_file_path",
                    version="%prog 1.0")
parser.add_option("-x", "--method",
                    type="string",
                    action="store",
                    dest="method",
                    default="PCA",
                    help="Reduce the dimensionality of the data. Choose one of : PCA, UMAP, t-SNE")
parser.add_option("-c", "--n-comps",
                    type="int",
                    action="store",
                    dest="n_comps",
                    default=50,
                    help="[PCA], Number of principal components to compute.")
parser.add_option("-s", "--svd-solver",
                    type="string",
                    action="store",
                    dest="svd_solver",
                    default="arpack",
                    help="[PCA], SVD solver to use. Choose one of : arpack (Default), randomized, auto.")
parser.add_option("-n", "--n-neighbors",
                    type="int",
                    action="store",
                    dest="n_neighbors",
                    default=15,
                    help="[UMAP], The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation.")
parser.add_option("-p", "--n-pcs",
                    type="int",
                    action="store",
                    dest="n_pcs",
                    default=30,
                    help="[UMAP], Use this many PCs.")
parser.add_option("-j", "--n-jobs",
                    type="int",
                    action="store",
                    dest="n_jobs",
                    default=1,
                    help="[t-SNE], The number of jobs. When set to None, automatically uses the number of cores.")           
(options, args) = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args[0]
FILE_PATH_OUT_BASENAME = os.path.splitext(args[1])[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except:
    raise Exception("Can only handle .h5ad files.")

#
# Transform the distribution of the data
#

if options.method == "PCA":
    # Run PCA
    sc.tl.pca(adata, n_comps=options.n_comps, svd_solver=options.svd_solver)
elif options.method == "UMAP":
    # Run UMAP
    sc.pp.neighbors(adata, n_neighbors=options.n_neighbors, n_pcs=options.n_pcs)
    sc.tl.umap(adata)
elif options.method == "t-SNE":
    # Run t-SNE
    from MulticoreTSNE import MulticoreTSNE as TSNE
    tsne = TSNE(n_jobs=options.n_jobs)
    adata.obsm['X_tsne'] = tsne.fit_transform(adata.X)
else:
    raise Exception("The dimensionality reduction method {} does not exist.".format(options.method))

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import scanpy as sc
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix
from collections import OrderedDict

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    nargs='+',
    type=argparse.FileType('r'),
    help='The path to the input h5ad files.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='The path to the merged h5ad output.'
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATHS_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

adatas = []

for FILE_PATH_IN in FILE_PATHS_IN:
    try:
        adata = sc.read_h5ad(FILE_PATH_IN.name)
        adatas.append(adata)
    except IOError:
        raise Exception("VSN ERROR: Wrong input format. Expects .h5ad files, got .{}".format(FILE_PATH_IN))

index_unique = None


def is_same_matrix(mat1, mat2):
    if type(mat1) == np.ndarray and type(mat2) == np.ndarray:
        return np.array_equal(mat1, mat2)
    return (mat1 != mat2).nnz == 0


if not all([is_same_matrix(mat1=adatas[0].X, mat2=_adata.X) for _adata in adatas]):
    raise Exception("VSN ERROR: adata.X slots are not the same between h5ad files.")
if not all([is_same_matrix(mat1=adatas[0].raw.X, mat2=_adata.raw.X) for _adata in adatas]):
    # Source: https://stackoverflow.com/questions/30685024/check-if-two-scipy-sparse-csr-matrix-are-equal
    raise Exception("VSN ERROR: adata.raw.X slots are not the same between h5ad files.")
if not all([np.array_equal(adatas[0].obs.index, _adata.obs.index) for _adata in adatas]):
    raise Exception("VSN ERROR: adata.obs.index slots are not the same between h5ad files.")
if not all([np.array_equal(adatas[0].var.index, _adata.var.index) for _adata in adatas]):
    raise Exception("VSN ERROR: adata.var.index slots are not the same between h5ad files.")

#
# Merge multiple h5ads with different clustering data
#
# Cannot use AnnData.concatenate since it's not joining multiple h5ad that have the same Cell IDs.
adata = ad.AnnData(
    X=adatas[0].X,
    obs=adatas[0].obs,
    var=adatas[0].var,
    uns=adatas[0].uns,
    raw=adatas[0].raw
)

adata.var.index = adata.var.index.astype(str)
adata = adata[:, np.sort(adata.var.index)]
print(f"Total number of cells: {adata.obs.shape[0]}, genes: {adata.var.shape[0]}.")


def get_clustering_algorithm(adata):
    if 'louvain' in adata.uns.keys():
        return 'louvain'
    if 'leiden' in adata.uns.keys():
        return 'leiden'
    raise Exception("No clustering found in the given AnnData.")

#######################
# Update the uns slot #
#######################


print("Update the uns slot in the merged h5ad...")
merged_uns = OrderedDict()

for _adata in adatas:
    for uns_key in _adata.uns.keys():
        clustering_algorithm = get_clustering_algorithm(adata=_adata)
        clustering_params = _adata.uns[clustering_algorithm]['params']
        new_key = f"{clustering_algorithm}_res{clustering_params['resolution']}__{uns_key}"
        if new_key in merged_uns:
            continue
        merged_uns[new_key] = _adata.uns[uns_key]

adata.uns = merged_uns

#######################
# Update the var slot #
#######################

#######################
# Update the obs slot #
#######################

# Since we merge multiple AnnData with different clustering results we rename them

print("Update the obs slot in the merged h5ad...")

if "louvain" in adata.obs.columns:
    adata.obs = adata.obs.drop(columns=["louvain"])
if "leiden" in adata.obs.columns:
    adata.obs = adata.obs.drop(columns=["leiden"])

for _adata in adatas:
    new_key = f"{clustering_algorithm}_res{clustering_params['resolution']}__{uns_key}"
    clustering_algorithm = get_clustering_algorithm(adata=_adata)
    clustering_params = _adata.uns[clustering_algorithm]['params']
    clustering = pd.DataFrame(
        {
            "clustering": _adata.obs[clustering_algorithm]
        },
        index=_adata.obs.index
    )
    clustering_name = f"{clustering_algorithm}_res{clustering_params['resolution']}"
    clustering.columns = [clustering_name]

    # Remove clustering if already exists
    if clustering_name in adata.obs:
        print(f"VSN MSG: Already existing {clustering_name} obs column in AnnData. Removing this column before joining clustering data...")
        del adata.obs[clustering_name]
    adata.obs = adata.obs.join(clustering)

####################
# Update obsm slot #
####################

# X_pca

if "X_pca" in adatas[0].obsm and all([np.array_equal(adatas[0].obsm["X_pca"], _adata.obsm["X_pca"]) for _adata in adatas]):
    print("Adding adata.obsm.X_pca to the merged h5ad...")
    adata.obsm["X_pca"] = adatas[0].obsm["X_pca"]
# X_tsne
if "X_tsne" in adatas[0].obsm and all([np.array_equal(adatas[0].obsm["X_tsne"], _adata.obsm["X_tsne"]) for _adata in adatas]):
    print("Adding adata.obsm.X_tsne to the merged h5ad...")
    adata.obsm["X_tsne"] = adatas[0].obsm["X_tsne"]
# X_umap
if "X_umap" in adatas[0].obsm and all([np.array_equal(adatas[0].obsm["X_umap"], _adata.obsm["X_umap"]) for _adata in adatas]):
    print("Adding adata.obsm.X_umap to the merged h5ad...")
    adata.obsm["X_umap"] = adatas[0].obsm["X_umap"]

####################
# Update varm slot #
####################

# PCs
if "PCs" in adatas[0].varm and all([np.array_equal(adatas[0].varm["PCs"], _adata.varm["PCs"]) for _adata in adatas]):
    print("Adding adata.varm.PCs to the merged h5ad...")
    adata.varm["PCs"] = adatas[0].varm["PCs"]

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

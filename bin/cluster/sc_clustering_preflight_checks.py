#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
import numpy as np
import itertools

parser = argparse.ArgumentParser(
    description='',
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input h5ad file.'
)

parser.add_argument(
    "-x", "--method",
    type=str,
    action="append",
    dest="methods",
    default=None,
    choices=["leiden", "louvain"],
    help="Clustering algorithm to check."
)

parser.add_argument(
    "-r", "--resolution",
    type=float,
    action="append",
    dest="resolutions",
    default=None,
    help="Resolution parameter to check."
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

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")


def check_neighborhood_graph_exists(adata):
    if "neighbors" not in adata.uns.keys():
        raise Exception(
            "The neighborhood graph of observations has not been computed.")


def has_single_cluster(adata, method):
    num_clusters = len(np.unique(adata.obs[method]))
    return num_clusters == 1


def has_any_singlet_cluster(adata, method):
    num_cell_by_clusters = adata.obs[method].value_counts()
    num_singlet_clusters = np.sum(num_cell_by_clusters == 1)
    return num_singlet_clusters > 0

#
# Check valid resolution when clustering the data
#


# [('louvain', 0.2), ('louvain', 0.4), ...]
args_grid = list(itertools.product(args.methods, args.resolutions))

args_valid = []
args_to_remove = []

print("VSN MSG: Checking valid methods and resolutions... ")

for method, resolution in args_grid:
    if method.lower() == "louvain":
        # Run Louvain clustering
        # Source: https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.louvain.html
        check_neighborhood_graph_exists(adata=adata)
        sc.tl.louvain(
            adata,
            resolution=resolution,
            random_state=args.seed
        )

    elif method.lower() == "leiden":
        # Run Leiden clustering
        # Source: https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.leiden.html
        check_neighborhood_graph_exists(adata=adata)
        sc.tl.leiden(
            adata,
            resolution=resolution,
            random_state=args.seed
        )
    else:
        raise Exception("VSN ERROR: The given clustering algorithm {} does not exist or is not implemeted.".format(args.method))

    if has_single_cluster(adata=adata, method=method):
        args_to_remove += [(method, resolution)]
        continue

    if has_any_singlet_cluster(adata=adata, method=method):
        args_to_remove += [(method, resolution)]
        continue

    args_valid += [(method, resolution)]

if len(args_valid) > 0:
    msg = "\n".join([f'- method: {m}, resolution: {r}'for m, r in args_valid])
    print(f"VSN MSG: The following arguments are valid: \n\033[92m{msg}\033[0m \033[91m.")

if len(args_to_remove) > 0:
    msg = "\n".join([f'- method: {m}, resolution: {r}'for m, r in args_to_remove])
    raise Exception(f"VSN ERROR: The following arguments should be removed from the grid:\n{msg}.")
else:
    print("VSN MSG: All preflight clustering checks have passed successfully ")

# I/O
# None

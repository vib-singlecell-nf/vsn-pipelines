#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='The path to the input h5ad file '
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='The path to the updated h5ad output'
)

parser.add_argument(
    '-p', "--x-pca",
    type=argparse.FileType('r'),
    dest="x_pca",
    required=False,
    help='The path the (compressed) TSV file containing the new PCA embeddings.'
)
parser.add_argument(
    '-r', "--empty-x",
    action="store_true",
    dest="remove_x",
    default=False,
    help="Empty X in the given h5ad file."
)

args = parser.parse_args()

FILE_PATH_IN = args.input.name
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

#
# Update the feature/observation-based metadata with all the columns present within the look-up table.
#

if args.remove_x:
    print("Emptying X slot of the given AnnData...")
    adata.X = csr_matrix(
        (adata.shape[0], adata.shape[1]),
        dtype=np.uint8
    )

if args.x_pca is not None:
    print("Updating X_pca slot of the given AnnData...")
    x_pca = pd.read_csv(
        filepath_or_buffer=args.x_pca,
        sep="\t",
        index_col=0
    ).values
    adata.obsm["X_pca"] = x_pca
    adata.uns["harmony"] = {
        "computed": True,
        "location": "AnnData.obsm['X_pca']"
    }


# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

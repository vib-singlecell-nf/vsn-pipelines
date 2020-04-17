#!/usr/bin/env python3

import argparse
import loompy as lp
import numpy as np
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

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

row_attrs = {
    "Gene": np.array(adata.var.index),
}
col_attrs = {
    "CellID": np.array(adata.obs.index),
    "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
    "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),
}

matrix = (adata.X).T

lp.create(
    filename=f"{FILE_PATH_OUT_BASENAME}.loom",
    layers=matrix if type(matrix) == np.ndarray else matrix.toarray(),
    row_attrs=row_attrs,
    col_attrs=col_attrs,
)

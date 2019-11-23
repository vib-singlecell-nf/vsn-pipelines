#!/usr/bin/env python
import os
from optparse import OptionParser

import loompy as lp
import numpy as np

import scanpy as sc

parser = OptionParser(
    usage="usage: %prog [options] h5ad_file_path",
    version="%prog 1.0"
)
(options, args) = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args[0]
FILE_PATH_OUT_BASENAME = os.path.splitext(args[1])[0]

try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

row_attrs = {
    "Gene": np.array(adata.var.index),
}
col_attrs = {
    "CellID": np.array(adata.obs.index),
    "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
    "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),
}

lp.create(
    filename=f"{FILE_PATH_OUT_BASENAME}.loom",
    layers=(adata.X).T.toarray(),
    row_attrs=row_attrs,
    col_attrs=col_attrs,
)

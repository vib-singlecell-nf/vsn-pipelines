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
    default="zscore_scale",
    help="Scale the data. Choose one of : zscore_scale"
)

parser.add_argument(
    "-M", "--max-sd",
    type=float,
    action="store",
    dest="max_sd",
    default=-1,
    help="Clip values greater than maximum number of standard deviation."
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]
SCALE__MAX_SD = args.max_sd if args.max_sd > 0 else None

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

#
# Transform the distribution of the data
#

if args.method == "zscore_scale":
    # scale each gene to unit variance, clip values exceeding SD max_sd.
    sc.pp.scale(
        adata,
        max_value=SCALE__MAX_SD
    )
else:
    raise Exception("VSN ERROR: Method does not exist.")

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

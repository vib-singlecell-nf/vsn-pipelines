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
    default="cpx",
    help="Normalize the data. Choose one of : cpx, regress_out"
)

parser.add_argument(
    "-f", "--counts-per-cell-after",
    type=int,
    action="store",
    dest="counts_per_cell_after",
    default=1e4,
    help="Multiplying factor used when running 'cpx' method."
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
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

#
# Normalize the data
#

adata.raw = adata

if args.method == 'cpx':
    # Total-count normalize (library-size correct) to '-r' reads/cell
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=args.counts_per_cell_after)
else:
    raise Exception("Method does not exist.")

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

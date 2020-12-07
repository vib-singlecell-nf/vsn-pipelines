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
    action="store",
    dest="method",
    default="linear_regression",
    help="Normalize the data. Choose one of : regress_out"
)

parser.add_argument(
    "-r", "--variable-to-regress-out",
    action="append",
    dest="vars_to_regress_out",
    default=None,
    help="Variable to regress out. To regress multiple variables, add that many -v arguments. Used when running 'regress_out"
)

parser.add_argument(
    "-j", "--n-jobs",
    type=int,
    action="store",
    dest="n_jobs",
    default=1,
    help="The number of jobs. When set to None, automatically uses the number of cores."
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
    raise Exception("VSN ERROR: Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

#
# Adjust the data
#

if args.method == 'linear_regression':
    # regress out variables (e.g.: total counts per cell and the percentage of mitochondrial genes expressed)
    sc.pp.regress_out(
        adata=adata,
        keys=args.vars_to_regress_out,
        n_jobs=args.n_jobs
    )
else:
    raise Exception("VSN ERROR: Method does not exist.")

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

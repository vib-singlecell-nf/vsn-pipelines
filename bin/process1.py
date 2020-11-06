#!/usr/bin/env python3

import argparse
import os
import scanpy as sc

parser = argparse.ArgumentParser(description='Template script')

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

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

################################################################################
# do some processing here...

print(adata)

################################################################################

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))


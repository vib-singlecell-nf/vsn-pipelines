#!/usr/bin/env python3

import argparse
import os
import scanpy as sc

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    nargs='+',
    type=argparse.FileType('r'),
    help='Input h5ad files.'
)

parser.add_argument(
    "-f", "--file-format",
    action="store",
    dest="format",
    default="h5ad",
    help="Concatenate the data. Choose one of : h5ad"
)

parser.add_argument(
    "-j", "--join",
    type=str,
    action="store",
    dest="join",
    default="inner",
    help="How to concatenate the multiple datasets. Choose one of : inner (intersect), outer (union)."
)

parser.add_argument(
    "-o", "--output",
    action="store",
    dest="output",
    default=None,
    help="Output file name."
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output)[0]

# I/O
files = []

if args.format == 'h5ad':
    for FILE_PATH_IN in args.input:
        try:
            FILE_PATH_IN = FILE_PATH_IN.name
            adata = sc.read_h5ad(filename=FILE_PATH_IN)
            files.append(adata)
        except IOError:
            raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(FILE_PATH_IN))

#
# Adjust the data
#

if args.format == 'h5ad':
    # Concatenate multiple h5ad files
    # Source: https://anndata.readthedocs.io/en/latest/anndata.AnnData.concatenate.html#anndata.AnnData.concatenate
    adata = files[0].concatenate(files[1:], join=args.join)
else:
    raise Exception("Concatenation of .{} files is not implemented.".format(args.format))

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

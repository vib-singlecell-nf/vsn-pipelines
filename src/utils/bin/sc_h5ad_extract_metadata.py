#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import scanpy as sc
import numpy as np

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='The path to the input h5ad file '
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='The path to the output containing cells IDs that will be used for applying the filter.'
)

parser.add_argument(
    '-a', '--axis',
    type=str,
    dest="axis",
    help='The axis defining the metadata which the given column_names will be extracted from. '
)

parser.add_argument(
    '-c', '--column-name',
    type=str,
    action="append",
    dest="column_names",
    help=""
)

args = parser.parse_args()

FILE_PATH_IN = args.input.name

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

#
# Extract the given column_names from the feature/observation-based metadata.
#

if args.axis == 'feature':
    metadata = adata.var[args.column_names]
elif args.axis == 'observation':
    raise Exception("VSN ERROR: Extracting the observation-based metadata is currently not implemented.")
else:
    raise Exception(f"Cannot extract from the {args.axis}-based metadata.")

# I/O
metadata.to_csv(
    path_or_buf=args.output,
    sep='\t',
    header=True,
    columns=args.column_names,
    index=False
)

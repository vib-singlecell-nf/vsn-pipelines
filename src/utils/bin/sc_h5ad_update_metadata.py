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
    '-m', "--additional-metadata",
    type=argparse.FileType('r'),
    dest="additional_metadata",
    required=True,
    help='The path the additional metadata used to update the metadata of the given input h5ad.'
)

parser.add_argument(
    '-a', '--axis',
    type=str,
    dest="axis",
    required=True,
    help='The axis defining the metadata which the given column_names will be extracted from. '
)

parser.add_argument(
    '-j', '--join-key',
    type=str,
    dest="join_key",
    required=True,
    help="The column name used to join the metadata with the given additional metadata."
)

parser.add_argument(
    '-i', '--index-column-name',
    type=str,
    dest="index_column_name",
    help="The column name to use as index for the metadata."
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

additional_metadata = pd.read_csv(
    filepath_or_buffer=args.additional_metadata,
    sep="\t",
    header=0
)

if args.axis == 'feature':
    adata.var = pd.merge(
        adata.var,
        additional_metadata,
        on=args.join_key
    )
    if args.index_column_name is not None:
        adata.var.set_index(args.index_column_name, inplace=True)
        adata.var.index.names = ['index']

elif args.axis == 'observation':
    raise Exception("VSN ERROR: Updating the observation-based metadata is currently not implemented.")

else:
    raise Exception(f"Cannot update the {args.axis}-based metadata.")


# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

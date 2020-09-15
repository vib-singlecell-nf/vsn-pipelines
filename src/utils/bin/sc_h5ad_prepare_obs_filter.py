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
    help='Either a .h5ad file with additional metadata (which can be used as filters) or .tsv file containing cell-based metadata.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='The path to the output containing cells IDs that will be used for applying the filter.'
)

parser.add_argument(
    '-m', '--method',
    type=str,
    dest="method",
    choices=['internal', 'external'],
    default='internal',
    help="The method to prepare the filters. Internal means, the input is expected to be a .h5ad otherwise it expects a .tsv."
)

parser.add_argument(
    '-i', '--sample-id',
    type=str,
    dest="sample_id",
    help="The sample ID to filter on."
)

parser.add_argument(
    '-s', '--sample-column-name',
    type=str,
    dest="sample_column_name",
    help="The column name containing the sample ID for each row in the cell meta data."
)

parser.add_argument(
    '-x', '--index-column-name',
    type=str,
    dest="index_column_name",
    help="The column name containing the index (unique identifier) for each row in the cell meta data."
)

parser.add_argument(
    '-f', '--filter-column-name',
    type=str,
    dest="filter_column_name",
    help="The column name that should be used to filter the given input anndata."
)

parser.add_argument(
    '-k', '--value-to-keep-from-filter-column',
    type=str,
    action="append",
    dest="values_to_keep_from_filter_column",
    help="The values to keep from the filter column."
)

args = parser.parse_args()

# I/O

#
# Extract cells that will be kept
#

if args.method == 'internal':
    # Expects h5ad file
    try:
        FILE_PATH_IN = args.input.name
        adata = sc.read_h5ad(filename=FILE_PATH_IN)
        metadata = adata.obs
    except IOError:
        raise Exception("VSN ERROR: Can only handle .h5ad files.")

elif args.method == 'external':
    metadata = pd.read_csv(
        filepath_or_buffer=args.input,
        header=0,
        index_col=args.index_column_name,
        sep="\t"
    )
else:
    raise Exception(f"VSN ERROR: The given method {args.method} is invalid.")

if args.sample_column_name not in metadata.columns:
    raise Exception(f"VSN ERROR: The meta data .tsv file expects a header with a required '{args.sample_column_name}' column.")
if args.filter_column_name not in metadata.columns:
    raise Exception(f"VSN ERROR: The meta data .tsv file expects a header with a required '{args.filter_column_name}' column.")


def is_bool(s):
    return s.lower() in ['true', 'false']


def str_to_bool(s):
    if s == 'True':
        return True
    elif s == 'False':
        return False
    else:
        raise ValueError


# Convert to boolean type if needed
values_to_keep_from_filter_column_formatted = [
    value_to_keep if not is_bool(s=value_to_keep) else str_to_bool(s=value_to_keep) for value_to_keep in args.values_to_keep_from_filter_column
]

filter_mask = np.logical_and(
    metadata[args.sample_column_name] == args.sample_id,
    metadata[args.filter_column_name].isin(values_to_keep_from_filter_column_formatted)
)
cells_to_keep = metadata.index[filter_mask]

# I/O
np.savetxt(args.output, cells_to_keep, fmt="%s")

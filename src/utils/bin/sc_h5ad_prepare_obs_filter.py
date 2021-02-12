#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import scanpy as sc
import numpy as np


def is_bool(s):
    return s.lower() in ['true', 'false']


def str_to_bool(s):
    if s == 'True':
        return True
    elif s == 'False':
        return False
    else:
        raise ValueError


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
    help="The method to prepare the filters. Internal means, the input is expected to be a .h5ad otherwise it expects a .tsv file."
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
    help="The column name containing the sample ID for each row in the cell metadata."
)

parser.add_argument(
    '-x', '--index-column-name',
    type=str,
    dest="index_column_name",
    help="The column name containing the index (unique identifier) for each row in the cell metadata."
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

    if args.index_column_name is None:
        raise Exception(f"VSN ERROR: Missing --index-column-name argument (indexColumnName param).")

    metadata = pd.read_csv(
        filepath_or_buffer=args.input,
        header=0,
        index_col=args.index_column_name,
        sep="\t"
    )
else:
    raise Exception(f"VSN ERROR: The given method {args.method} is invalid.")

filter_mask = None
values_to_keep_from_filter_column_formatted = None

if args.filter_column_name in metadata.columns:
    if args.values_to_keep_from_filter_column is None:
        raise Exception(f"VSN ERROR: Missing --value-to-keep-from-filter-column argument (valuesToKeepFromFilterColumn param).")
    # Convert to boolean type if needed
    values_to_keep_from_filter_column_formatted = [
        value_to_keep if not is_bool(s=value_to_keep) else str_to_bool(s=value_to_keep) for value_to_keep in args.values_to_keep_from_filter_column
    ]

if args.method == 'internal' or len(np.unique(metadata.index)) != len(metadata.index):

    # Check if index are all unique, in that case index only may be used
    if len(np.unique(metadata.index)) == metadata.shape[0]:
        # Check existence of 1-to-many relationships between adata.obs.index and metadata.index
        if np.sum(np.isin(adata.obs.index, metadata.index)) == adata.obs.shape[0]:
            print(f"VSN MSG: Creating a filter mask based only on '{args.filter_column_name}'...")
            filter_mask = metadata[args.filter_column_name].isin(values_to_keep_from_filter_column_formatted)
        else:
            raise Exception(f"VSN ERROR: There is not a 1-to-1 relationship between the index elements of the dataset {args.input.name} and index elements of the metadata {args.cell_meta_data_file_paths[0].name}.")
    else:
        if args.sample_column_name is None:
            raise Exception(f"VSN ERROR: Missing --sample-column-name argument (sampleColumnName param).")

        if args.sample_column_name not in metadata.columns:
            raise Exception(f"VSN ERROR: Missing '{args.sample_column_name}' column in obs slot of the given h5ad input file.")

        if values_to_keep_from_filter_column_formatted is None:
            raise Exception(f"VSN ERROR: Missing --filter-column-name argument (filterColumnName param) and/or --value-to-keep-from-filter-column (valuesToKeepFromFilterColumn param). These are required since the '{args.index_column_name}' index column does not contain unique values.")

        print(f"VSN MSG: Creating a filter mask based on '{args.sample_column_name}' and '{args.filter_column_name}'...")
        filter_mask = np.logical_and(
            metadata[args.sample_column_name] == args.sample_id,
            metadata[args.filter_column_name].isin(values_to_keep_from_filter_column_formatted)
        )
else:
    if values_to_keep_from_filter_column_formatted is not None:
        print(f"VSN MSG: Creating a filter mask based only on '{args.filter_column_name}' filter column...")
        filter_mask = metadata[args.filter_column_name].isin(values_to_keep_from_filter_column_formatted)
    else:
        print(f"VSN MSG: No filter mask: use all the cells from the given cell-based metadata .tsv filter file...")
        filter_mask = np.in1d(metadata.index, metadata.index)

cells_to_keep = metadata.index[filter_mask]

# I/O
np.savetxt(args.output, cells_to_keep, fmt="%s")

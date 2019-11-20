#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import scanpy as sc
import numpy as np

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "cell_meta_data_file_path",
    type=argparse.FileType('r'),
    help='A file path to meta data (TSV) for each cell where values from a column could be used as filter.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='The path to the output containing cells IDs that will be used for applying the filter.'
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
    help="The column name containing the sample ID for each cell entry in the cell meta data."
)

parser.add_argument(
    '-b', '--barcode-column-name',
    type=str,
    dest="barcode_column_name",
    help="The column name containing the barcode for each cell entry in the cell meta data."
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

metadata = pd.read_csv(
    filepath_or_buffer=args.cell_meta_data_file_path,
    header=0,
    sep="\t"
)

if args.sample_column_name not in metadata.columns:
    raise Exception(f"The meta data TSV file expects a header with a required '{args.sample_column_name}' column.")
if args.barcode_column_name not in metadata.columns:
    raise Exception(f"The meta data TSV file expects a header with a required '{args.barcode_column_name}' column.")
if args.filter_column_name not in metadata.columns:
    raise Exception(f"The meta data TSV file expects a header with a required '{args.filter_column_name}' column.")

sample_mask = np.logical_and(
    metadata[args.sample_column_name] == args.sample_id,
    metadata[args.filter_column_name].isin(args.values_to_keep_from_filter_column)
)
sample_cells = metadata[args.barcode_column_name][sample_mask]

# I/O
sample_cells.to_csv(path=args.output, index=False)

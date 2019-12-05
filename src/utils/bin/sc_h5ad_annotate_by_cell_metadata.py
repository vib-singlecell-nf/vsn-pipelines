#!/usr/bin/env python

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
    help='The h5ad file path which to annotate the cells from.'
)

parser.add_argument(
    "cell_meta_data_file_path",
    type=argparse.FileType('r'),
    help='The file path to meta data (TSV with header) for each cell where values from a column could be used to annotate the cells.'
)

parser.add_argument(
    "-o", "--output",
    type=argparse.FileType('w'),
    required=True,
    help='The h5ad file path to which the AnnData updated with annotated cells will be saved.'
)

parser.add_argument(
    '-i', '--sample-id',
    type=str,
    dest="sample_id",
    help="The sample ID to annotate the cells from."
)

parser.add_argument(
    '-s', '--sample-column-name',
    type=str,
    dest="sample_column_name",
    help="The column name containing the sample ID for each cell entry in the cell meta data."
)

parser.add_argument(
    '-b', '--index-column-name',
    type=str,
    required=True,
    dest="index_column_name",
    help="The column of the given cell_meta_data_file_path that should be used as index."
)

parser.add_argument(
    '-a', '--annotation-column-name',
    type=str,
    action="append",
    required=True,
    dest="annotation_column_names",
    help="The columns of the given cell_meta_data_file_path that will be used as annotation to annotate the cells."
)

args = parser.parse_args()

# Define the arguments properly
args = parser.parse_args()
FILE_PATH_IN = args.input.name
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except Exception:
    raise Exception("Can only handle .h5ad files.")

#
# Annotate the data
#

metadata = pd.read_csv(
    filepath_or_buffer=args.cell_meta_data_file_path,
    sep="\t",
    header=0
)

for annotation_column_name in args.annotation_column_names:
    # Extract the metadata for the cells from the adata
    metadata_subset = metadata[
        np.logical_and(
            metadata[args.sample_column_name] == args.sample_id,  # Taking sample into consideration is important here (barcode collision between samples might happen!)
            metadata[args.index_column_name].isin(adata.obs.index.values))
    ]
    # Annotate
    adata.obs[annotation_column_name] = None
    adata.obs[annotation_column_name][metadata_subset[args.index_column_name]] = metadata_subset[annotation_column_name]

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

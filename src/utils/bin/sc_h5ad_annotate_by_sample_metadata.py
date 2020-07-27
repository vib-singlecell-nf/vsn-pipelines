#!/usr/bin/env python3

import argparse
import os
import pandas as pd
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
    "-t", "--type",
    type=str,
    action="store",
    dest="type",
    default="sample",
    help="Are the entries of the meta data sample-based or cell-based? Choose one of: sample or cell"
)

parser.add_argument(
    "-m", "--meta-data-file-path",
    type=str,
    action="store",
    dest="meta_data_file_path",
    help="Path to the meta data. It expects a tabular separated file (.tsv) with header and a required 'id' column."
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]
SAMPLE_NAME = os.path.splitext(FILE_PATH_OUT_BASENAME)[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

#
# Annotate the data
#

metadata = pd.read_csv(
    filepath_or_buffer=args.meta_data_file_path,
    sep="\t"
)

if 'id' not in metadata.columns:
    raise Exception("VSN ERROR: The meta data TSV file expects a header with a required 'id' column.")

if args.type == "sample":
    sample_info = metadata[metadata.id == SAMPLE_NAME]

    if len(sample_info) == 0:
        raise Exception(f"VSN ERROR: The meta data TSV file does not contain sample ID '{SAMPLE_NAME}'.")
    elif len(sample_info) > 1:
        raise Exception(f"VSN ERROR: The meta data TSV file contains more than one line with sample ID '{SAMPLE_NAME}'.")

    for (column_name, column_data) in sample_info.iteritems():
        adata.obs[column_name] = column_data.values[0]
else:
    raise Exception("VSN ERROR: This meta data type {} is not implemented".format(args.type))

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import scanpy as sc

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help=''
)

parser.add_argument(
    "-o", "--output",
    type=argparse.FileType('w'),
    help=''
)

parser.add_argument(
    '-f', '--filter-file-path',
    type=argparse.FileType('r'),
    action="append",
    dest="filter_file_paths",
    help=""
)

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
# Subset the h5ad using the given cell IDs
#

obs_to_keep = []

for filter_file_path in args.filter_file_paths:
    obs_to_keep.extend(
        pd.read_csv(filepath_or_buffer=filter_file_path, header=None)[0].values
    )

print(f"Dimension of pre-filtered AnnData: {adata.shape}")
adata_filtered = adata[obs_to_keep, :]
print(f"Dimension of post-filtered AnnData: {adata_filtered.shape}")

# I/O
adata_filtered.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

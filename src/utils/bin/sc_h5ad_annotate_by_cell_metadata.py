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
    help='The h5ad file path which to annotate the cells from.'
)

parser.add_argument(
    "cell_meta_data_file_path",
    type=argparse.FileType('r'),
    help='The file path to metadata (.tsv with header) for each cell where values from a column could be used to annotate the cells.'
)

parser.add_argument(
    "-o", "--output",
    type=argparse.FileType('w'),
    required=True,
    help='The h5ad file path to which the AnnData updated with annotated cells will be saved.'
)

parser.add_argument(
    '-m', '--method',
    type=str,
    dest="method",
    choices=['aio', 'obo'],
    required=True,
    default='obo',
    help="The method to use to annotate the cells in the AnnData."
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
    required=False,
    help="The column name containing the sample ID for each cell entry in the cell metadata."
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
    required=False,
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
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

#
# Annotate the data
#

metadata = pd.read_csv(
    filepath_or_buffer=args.cell_meta_data_file_path,
    sep="\t",
    header=0,
    index_col=args.index_column_name
)

if args.method == 'obo':
    print("Annotating using the OBO method...")
    if "sample_id" in metadata.columns:
        # Drop the sample_id if already contained in adata.obs otherwise join needs a suffix
        # but we want to keep sample_id as column name w/o any suffix
        metadata.drop(
            columns=['sample_id'],
            inplace=True
        )
    if args.annotation_column_names is not None:
        metadata = metadata.filter(items=args.annotation_column_names)

    # Check if all elements from the metadata subset are present in the given input file (h5ad file)
    if np.sum(np.isin(adata.obs.index, metadata.index)) != len(adata.obs):
        raise Exception(f"VSN ERROR: Make sure the cell IDs from the given input h5ad {FILE_PATH_IN} exist in the column {args.index_column_name} of the following metadata file ({args.cell_meta_data_file_path.name}) you provided in params.sc.cell_annotate.cellMetaDataFilePath.")

    adata.obs = adata.obs.join(
        other=metadata
    )
elif args.method == 'aio':
    print("Annotating using the AIO method...")
    if args.sample_column_name is None:
        raise Exception("VSN ERROR: Missing the --sample-column-name (sampleColumnName param) which is required for the 'aio' method.")

    for annotation_column_name in args.annotation_column_names:
        # Extract the metadata for the cells from the adata
        metadata_subset = metadata[
            np.logical_and(
                metadata[args.sample_column_name] == args.sample_id,  # Taking sample into consideration is important here (barcode collision between samples might happen!)
                metadata.index.isin(adata.obs.index.values))
        ]
        # Check if all elements from the metadata subset are present in the given input file (h5ad file)
        if np.sum(np.isin(adata.obs.index, metadata_subset.index)) != len(adata.obs):
            raise Exception(f"VSN ERROR: Make sure the sample IDs inferred from the data files (e.g.: {args.sample_id}) exist in the column {args.sample_column_name} of the following metadata file ({args.cell_meta_data_file_path.name}) you provided in params.sc.cell_annotate.cellMetaDataFilePath.")

        # Annotate
        adata.obs[annotation_column_name] = None
        adata.obs.loc[metadata_subset.index, annotation_column_name] = metadata_subset[annotation_column_name]
else:
    raise Exception(f"VSN ERROR: The given method '{args.method}' is not valid.")

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

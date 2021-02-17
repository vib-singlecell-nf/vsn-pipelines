#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import scanpy as sc
import numpy as np
import warnings

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='The h5ad file path which to annotate the cells from.'
)

parser.add_argument(
    "cell_meta_data_file_paths",
    type=argparse.FileType('r'),
    nargs='+',
    help='The file path(s) to metadata (.tsv with header) for each cell where values from a column could be used to annotate the cells.'
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

if args.method == 'obo':
    print("Annotating using the OBO method...")
    for cell_meta_data_file_path in args.cell_meta_data_file_paths:
        metadata = pd.read_csv(
            filepath_or_buffer=cell_meta_data_file_path,
            sep="\t",
            header=0,
            index_col=args.index_column_name
        )
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
            warnings.warn(f"VSN WARNING: Incomplete join between given .h5ad ({FILE_PATH_IN}, num cells: {adata.obs.shape[0]}) and given metadata file ({cell_meta_data_file_path}, num cells: {metadata.shape[0]}).")
        adata.obs = adata.obs.join(
            other=metadata
        )
elif args.method == 'aio':
    if len(args.cell_meta_data_file_paths) > 1:
        raise Exception("VSN ERROR: Using multiple metadata files is currently not supported with the AIO method.")
    
    metadata = pd.read_csv(
        filepath_or_buffer=args.cell_meta_data_file_paths[0].name,
        sep="\t",
        header=0,
        index_col=args.index_column_name
    )
    print("Annotating using the AIO method...")

    if len(np.unique(adata.obs.index)) != adata.obs.shape[0]:
        raise Exception("VSN ERROR: Cells of the given h5ad are not unique.")

    for annotation_column_name in args.annotation_column_names:

        # Extract the metadata for the cells from the adata
        # Check if index are all unique, in that case index only may be used
        if len(np.unique(metadata.index)) == metadata.shape[0]:
            # Check existence of 1-to-many relationships between adata.obs.index and metadata.index
            num_matching_cells = np.sum(np.isin(adata.obs.index, metadata.index))
            num_matching_cells/adata.obs.shape[0]

            # proceed as long as we can annotate "most" of the cells in the metadata file:
            if num_matching_cells/adata.obs.shape[0] >= 0.8:
                print(f"VSN MSG: subsetting metadata based on only '{args.index_column_name}' column of the metadata {args.cell_meta_data_file_paths[0].name}.")
                print(f"VSN MSG: Annotating {num_matching_cells} ({round(num_matching_cells/adata.obs.shape[0],2)*100}%) of {args.sample_id} cells in the adata with information from the metadata.")
                metadata_subset = metadata[metadata.index.isin(adata.obs.index.values)]
            else:
                raise Exception(f"VSN ERROR: There is a dimension mismatch between the dataset {args.input.name} and the metadata {args.cell_meta_data_file_paths[0].name}: expected [80% of {len(adata.obs)}] but got {num_matching_cells} cells matching.")
        else:
            if args.sample_column_name is None:
                raise Exception("VSN ERROR: Missing the --sample-column-name (sampleColumnName param) which is required for the 'aio' method.")

            print(f"VSN MSG: subsetting metadata based on index on '{args.sample_column_name}' and '{args.index_column_name}' columns of the metadata {args.cell_meta_data_file_paths[0].name}.")
            metadata_subset = metadata[
                np.logical_and(
                    metadata[args.sample_column_name] == args.sample_id,  # Taking sample into consideration is important here (barcode collision between samples might happen!)
                    metadata.index.isin(adata.obs.index.values))
            ]
            # Check if all elements from the metadata subset are present in the given input file (h5ad file)
            num_matching_cells = np.sum(np.isin(adata.obs.index, metadata_subset.index))

            if num_matching_cells != len(adata.obs):
                raise Exception(f"VSN ERROR: Dimensions mismatch between {args.input.name} and {args.cell_meta_data_file_paths[0].name}: expected {len(adata.obs)} but got {num_matching_cells} cells matching. Make sur all cells from metadata file can be found in the data and/or make sure the sample IDs inferred from the data files (e.g.: {args.sample_id}) exist in the column {args.sample_column_name} of the following metadata file ({args.cell_meta_data_file_paths[0].name}) you provided in params.utils.cell_annotate.cellMetaDataFilePath.")

        # Annotate
        adata.obs[annotation_column_name] = None
        adata.obs.loc[metadata_subset.index, annotation_column_name] = metadata_subset[annotation_column_name]
else:
    raise Exception(f"VSN ERROR: The given method '{args.method}' is not valid.")

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

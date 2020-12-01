#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import scanpy as sc
from difflib import SequenceMatcher

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
    "-s", "--sample-id",
    type=str,
    action="store",
    dest="sample_id",
    help="The sample ID of the given h5ad file."
)

parser.add_argument(
    "-m", "--method",
    type=str,
    action="store",
    dest="method",
    default="sample",
    choices=["sample", "sample+"],
    help="Method used to do the sample-based annotation. Besides the sample ID, if another annotation column is required to annotate the cells, use 'sample+'."
)

parser.add_argument(
    "-f", "--metadata-file-path",
    type=str,
    action="store",
    dest="metadata_file_path",
    help="Path to the metadata. It expects a tabular separated file (.tsv) with header and a required 'id' column."
)

parser.add_argument(
    "-t", "--sample-column-name",
    type=str,
    action="store",
    dest="sample_column_name",
    help="The name of the column of the given metadata containing the sample ID."
)

parser.add_argument(
    '-i', '--adata-comp-index-column-name',
    type=str,
    action="append",
    dest="adata_comp_index_column_names",
    help="The names of complementary columns to use as index with the sample ID/name to annotate the cells. These column names should exist in the obs slot of the given h5ad input and have a 1-to-1 mapping with the -j/--metadata-comp-index-column-name."
)

parser.add_argument(
    '-j', '--metadata-comp-index-column-name',
    type=str,
    action="append",
    dest="metadata_comp_index_column_names",
    help="The names of complementary columns to use as index with the sample ID/name to annotate the cells.  These column names should exist in the obs slot of the given h5ad input and have a 1-to-1 mapping with the -i/--adata-comp-index-column-name"
)

parser.add_argument(
    '-a', '--annotation-column-name',
    type=str,
    action="append",
    dest="annotation_column_names",
    help="The names of columns from the given metadata to keep."
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]
SAMPLE_NAME = args.sample_id

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

#
# Annotate the data
#

if "sample_id" not in adata.obs.columns:
    raise Exception("VSN ERROR: Missing sample_id column in the obs slot of the AnnData of the given h5ad.")

if args.sample_column_name is None:
    raise Exception("VSN ERROR: Missing --sample-column-name argument (sampleColumnName param in sample_annotate config)")

metadata = pd.read_csv(
    filepath_or_buffer=args.metadata_file_path,
    sep="\t"
)


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


if args.sample_column_name not in metadata.columns:
    raise Exception(f"VSN ERROR: Missing '{args.sample_column_name}'' column in the given metadata file.")

sample_info = metadata[metadata[args.sample_column_name] == SAMPLE_NAME]
sample_scores = [similar(index_entry, SAMPLE_NAME) for index_entry in metadata[args.sample_column_name]]

if all(sample_score < 0.5 for sample_score in sample_scores):
    # Skip annotation for this sample
    print(f"Skipping annotation for {SAMPLE_NAME}.")
else:
    if len(sample_info) == 0:
        raise Exception(f"VSN ERROR: The metadata .tsv file does not contain sample ID '{SAMPLE_NAME}'.")
    elif args.method == "sample" and len(sample_info) > 1:
        raise Exception(f"VSN ERROR: The metadata .tsv file contains duplicate entries with the sample ID '{SAMPLE_NAME}'. Fix your metadata or use the 'sample+' method.")

    if args.method == "sample":
        for (column_name, column_data) in sample_info.iteritems():
            adata.obs[column_name] = column_data.values[0]
    elif args.method == "sample+":

        if args.adata_comp_index_column_names is None or args.metadata_comp_index_column_names is None:
            raise Exception("VSN ERROR: compIndexColumnNames param is missing in the sample_annotate config.")

        new_obs = pd.merge(
            adata.obs,
            sample_info,
            left_on=["sample_id"] + args.adata_comp_index_column_names,
            right_on=[args.sample_column_name] + args.metadata_comp_index_column_names
        )

        if new_obs.isnull().values.any():
            raise Exception("VSN ERROR: Merged adata.obs not complete, some NaN values detected.")

        # Update the obs slot of the AnnData
        adata.obs = new_obs
    else:
        raise Exception(f"VSN ERROR: Unrecognized method {args.method}.")

    if args.annotation_column_names is not None and len(args.annotation_column_names) > 0:
        adata.obs = adata.obs[args.annotation_column_names]

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

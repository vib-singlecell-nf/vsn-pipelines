#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import scanpy as sc
import numpy as np
import json
import re
from scipy.sparse import csr_matrix

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='The path to the input h5ad file '
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='The path to the updated h5ad output'
)

parser.add_argument(
    '-p', "--x-pca",
    type=argparse.FileType('r'),
    dest="x_pca",
    required=False,
    help='The path the (compressed) .tsv file containing the new PCA embeddings.'
)
parser.add_argument(
    '-r', "--empty-x",
    action="store_true",
    dest="remove_x",
    default=False,
    help="Empty X in the given h5ad file."
)
parser.add_argument(
    '-c', "--obs-column-mapper",
    type=str,
    action="store",
    dest="obs_column_mapper",
    help="Rename the columns in the obs slot of AnnData of the given input h5ad file using this mapper. This will be performed 1st."
)
parser.add_argument(
    '-v', "--obs-column-value-mapper",
    type=str,
    action="store",
    dest="obs_column_value_mapper",
    help="Rename the values in the obs slot of AnnData of the given input h5ad file using this mapper. This will be performed 2nd."
)
parser.add_argument(
    '-z', "--obs-column-to-remove",
    action="append",
    dest="obs_column_to_remove",
    help="Remove that column from the obs slot of the AnnData of the given input h5ad file. This will be performed 3rd."
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
# Update the feature/observation-based metadata
#

if args.remove_x:
    print("Emptying X slot of the given AnnData...")
    adata.X = csr_matrix(
        (adata.shape[0], adata.shape[1]),
        dtype=np.uint8
    )

if args.x_pca is not None:
    print("Updating X_pca slot of the given AnnData...")
    x_pca = pd.read_csv(
        filepath_or_buffer=args.x_pca,
        sep="\t",
        index_col=0
    ).values
    adata.obsm["X_pca"] = x_pca
    adata.uns["harmony"] = {
        "computed": True,
        "location": "AnnData.obsm['X_pca']"
    }

if args.obs_column_mapper is not None:
    print("Renaming the columns in the obs slot of the given AnnData...")
    # Renaming will be performed as in https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.rename.html
    # Mapper e.g.: {"ColumnA": "NewColumnAName", "ColumnB": "NewColumnBName"}
    obs_column_mapper = json.loads(args.obs_column_mapper)
    print(obs_column_mapper)
    adata.obs = adata.obs.rename(columns=obs_column_mapper)

if args.obs_column_value_mapper is not None:
    print("Renaming the values in the obs slot of the given AnnData...")
    # Renaming will be performed as in https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.replace.html
    # Mapper e.g.: {'ColumnA': {'OldValueToReplace_1': 'NewValue_1', 'OldValueToReplace_2': 'NewValue_2}
    obs_column_value_mapper = json.loads(args.obs_column_value_mapper)
    if type(obs_column_value_mapper) is dict:
        print("...using dict-like approach.")
        adata.obs = adata.obs.replace(to_replace=obs_column_value_mapper)
    else:
        print("...using RegExp-like approach.")
        for mapping in obs_column_value_mapper:
            adata.obs = adata.obs.replace(to_replace=rf'{mapping["from"]}', value=f'{mapping["to"]}', regex=True)

if args.obs_column_to_remove is not None:
    print("Removing some columns in the obs slot of the given AnnData...")
    for obs_column_pattern_to_remove in args.obs_column_to_remove:
        remove_columns_mask = [bool(re.match(obs_column_pattern_to_remove, obs_column)) for obs_column in adata.obs.columns]
        columns_to_remove = adata.obs.columns[remove_columns_mask]
        adata.obs = adata.obs.drop(columns_to_remove, axis=1)

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

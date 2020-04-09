#!/usr/bin/env python3

import argparse
import os
import re
import scanpy as sc
import numpy as np


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

in_formats = [
    '10x_cellranger_mex',
    '10x_cellranger_h5',
    'tsv',
    'csv'
]

out_formats = [
    'h5ad'
]

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=str,
    help='Input h5ad file.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output h5ad file.'
)

parser.add_argument(
    "-i", "--input-format",
    action="store",
    dest="input_format",
    default="",
    help="Input format of the file to be converted. Choose one of: {}.".format(', '.join(in_formats))
)

parser.add_argument(
    "-s", "--sample-id",
    type=str,
    dest="sample_id",
    default=None,
    action='store'
)

parser.add_argument(
    "-t", "--tag-cell-with-sample-id",
    type=str2bool,
    action="store",
    dest="tag_cell_with_sample_id",
    default=False,
    help="Tag each cell with the given sample_id."
)

parser.add_argument(
    "-u", "--make-var-index-unique",
    type=str2bool,
    action="store",
    dest="make_var_index_unique",
    default=False,
    help="Make the var index unique of the AnnData."
)

parser.add_argument(
    "-o", "--output-format",
    action="store",  # optional because action defaults to "store"
    dest="output_format",
    default="",
    help="Output format which the file should be converted to. Choose one of: {}.".format(', '.join(out_formats))
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]
INPUT_FORMAT = args.input_format
OUTPUT_FORMAT = args.output_format


def check_10x_cellranger_mex_path(path):
    # Sanity checks
    if not os.path.isdir(path):
        raise Exception(
            "Expecting a directory with an .mtx file when converting from 10x_cellranger_mex. {} does not seem to be one.".format(
                path
            )
        )
    if not os.path.exists(path):
        raise Exception("VSN ERROR: The given directory {} does not exist.".format(path))
    if not (
        not os.path.exists(os.path.join(path, "matrix.mtx")) or not os.path.exists(os.path.join(path, "matrix.mtx.gz"))
    ):
        raise Exception(
            "The given directory {} is not a proper 10xGenomics CellRanger folder. No .mtx[.gz] file found.".format(
                path
            )
        )


def add_sample_id(adata, args):
    # Annotate the file with the sample ID
    adata.obs["sample_id"] = args.sample_id
    return adata


if INPUT_FORMAT == '10x_cellranger_mex' and OUTPUT_FORMAT == 'h5ad':
    check_10x_cellranger_mex_path(path=FILE_PATH_IN)
    # Convert
    print("Reading 10x data from MEX format...")
    adata = sc.read_10x_mtx(
        FILE_PATH_IN,  # the directory with the `.mtx` file
        var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
        cache=False
    )
    adata = add_sample_id(
        adata=adata,
        args=args
    )
    # If is tag_cell_with_sample_id is given, add the sample ID as suffix
    if args.tag_cell_with_sample_id:
        adata.obs.index = map(lambda x: re.sub('-[0-9]+', f"-{args.sample_id}", x), adata.obs.index)
    adata.var.index = adata.var.index.astype(str)
    # Check if var index is unique
    if len(np.unique(adata.var.index)) < len(adata.var.index) and not args.make_var_index_unique:
        raise Exception("VSN ERROR: AnnData var index is not unique. This can be fixed by making it unique. To do so update the following param 'makeVarIndexUnique = true' (under params.sc.sc_file_converter) in your config.")
    if len(np.unique(adata.var.index)) < len(adata.var.index) and args.make_var_index_unique:
        adata.var_names_make_unique()
        print("Making AnnData var index unique...")
    # Sort var index
    adata = adata[:, np.sort(adata.var.index)]
    print("Writing 10x data to h5ad...")
    adata.write_h5ad(filename="{}.h5ad".format(FILE_PATH_OUT_BASENAME))

elif INPUT_FORMAT == '10x_cellranger_h5' and OUTPUT_FORMAT == 'h5ad':
    if not os.path.exists(FILE_PATH_IN):
        raise Exception("VSN ERROR: The given file {} does not exist.".format(FILE_PATH_IN))
    # Convert
    print("Reading 10x data from HDF5 format...")
    adata = sc.read_10x_h5(
        FILE_PATH_IN
    )
    adata = add_sample_id(
        adata=adata,
        args=args
    )
    # If is tag_cell_with_sample_id is given, add the sample ID as suffix
    if args.tag_cell_with_sample_id:
        adata.obs.index = map(lambda x: re.sub('-[0-9]+', f"-{args.sample_id}", x), adata.obs.index)
    adata.var.index = adata.var.index.astype(str)
    # Check if var index is unique
    if len(np.unique(adata.var.index)) < len(adata.var.index) and not args.make_var_index_unique:
        raise Exception("VSN ERROR: AnnData var index is not unique. This can be fixed by making it unique. To do so update the following param 'makeVarIndexUnique = true' (under params.sc.sc_file_converter) in your config.")
    if len(np.unique(adata.var.index)) < len(adata.var.index) and args.make_var_index_unique:
        adata.var_names_make_unique()
        print("Making AnnData var index unique...")
    # Sort var index
    adata = adata[:, np.sort(adata.var.index)]
    print("Writing 10x data to h5ad...")
    adata.write_h5ad(filename="{}.h5ad".format(FILE_PATH_OUT_BASENAME))

elif INPUT_FORMAT in ['tsv', 'csv'] and OUTPUT_FORMAT == 'h5ad':
    if INPUT_FORMAT == 'tsv':
        delim = '\t'
    elif INPUT_FORMAT == 'csv':
        delim = ','

    adata = sc.read_csv(
        FILE_PATH_IN,
        delimiter=delim,
        first_column_names=True
    ).T
    adata = add_sample_id(
        adata=adata,
        args=args
    )
    # If is tag_cell_with_sample_id is given, add the sample ID as suffix
    if args.tag_cell_with_sample_id:
        adata.obs.index = map(lambda x: re.sub('-[0-9]+', f"-{args.sample_id}", x), adata.obs.index)
    adata.var.index = adata.var.index.astype(str)
    # Check if var index is unique
    if len(np.unique(adata.var.index)) < len(adata.var.index) and not args.make_var_index_unique:
        raise Exception("VSN ERROR: AnnData var index is not unique. This can be fixed by making it unique. To do so update the following param 'makeVarIndexUnique = true' (under params.sc.sc_file_converter) in your config.")
    if len(np.unique(adata.var.index)) < len(adata.var.index) and args.make_var_index_unique:
        adata.var_names_make_unique()
        print("Making AnnData var index unique...")
    # Sort var index
    adata = adata[:, np.sort(adata.var.index)]
    adata.write_h5ad(filename="{}.h5ad".format(FILE_PATH_OUT_BASENAME))

elif INPUT_FORMAT == 'h5ad' and OUTPUT_FORMAT == 'h5ad':
    adata = sc.read_h5ad(
        FILE_PATH_IN
    )
    adata = add_sample_id(
        adata=adata,
        args=args
    )
    # If is tag_cell_with_sample_id is given, add the sample ID as suffix
    if args.tag_cell_with_sample_id:
        adata.obs.index = map(lambda x: re.sub('-[0-9]+', f"-{args.sample_id}", x), adata.obs.index)
    adata.var.index = adata.var.index.astype(str)
    # Check if var index is unique
    if len(np.unique(adata.var.index)) < len(adata.var.index) and not args.make_var_index_unique:
        raise Exception("VSN ERROR: AnnData var index is not unique. This can be fixed by making it unique. To do so update the following param 'makeVarIndexUnique = true' (under params.sc.sc_file_converter) in your config.")
    if len(np.unique(adata.var.index)) < len(adata.var.index) and args.make_var_index_unique:
        adata.var_names_make_unique()
        print("Making AnnData var index unique...")
    # Sort var index
    adata = adata[:, np.sort(adata.var.index)]
    adata.write_h5ad(filename="{}.h5ad".format(FILE_PATH_OUT_BASENAME))
else:
    raise Exception(
        "File format conversion {0} --> {1} hasn't been implemented yet.".format(INPUT_FORMAT, OUTPUT_FORMAT))

print("Done!")

#!/usr/bin/env python3

import argparse
import os
import scanpy as sc

in_formats = [
    '10x_mtx',
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

if INPUT_FORMAT == '10x_mtx' and OUTPUT_FORMAT == 'h5ad':
    # Sanity checks
    if not os.path.isdir(FILE_PATH_IN):
        raise Exception(
            "Expecting a directory with an .mtx file when converting from 10x_mtx. {} does not seem to be one.".format(
                FILE_PATH_IN
            )
        )
    if not os.path.exists(FILE_PATH_IN):
        raise Exception("The given directory {} does not exist.".format(FILE_PATH_IN))
    if not (
        not os.path.exists(os.path.join(FILE_PATH_IN, "matrix.mtx")) or not os.path.exists(os.path.join(FILE_PATH_IN, "matrix.mtx.gz"))
    ):
        raise Exception(
            "The given directory {} is not a proper 10xGenomics CellRanger folder. No .mtx[.gz] file found.".format(
                FILE_PATH_IN
            )
        )
    # Convert
    print("Reading 10x data...")
    adata = sc.read_10x_mtx(
        FILE_PATH_IN,  # the directory with the `.mtx` file
        var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
        cache=False
    )
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
    adata.write_h5ad(filename="{}.h5ad".format(FILE_PATH_OUT_BASENAME))


else:
    raise Exception(
        "File format conversion {0} --> {1} hasn't been implemented yet.".format(INPUT_FORMAT, OUTPUT_FORMAT))

print("Done!")

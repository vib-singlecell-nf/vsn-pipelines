#!/usr/bin/env python3

import argparse
import os
import scanpy as sc


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


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
    "-x", "--method",
    type=str,
    action="store",
    dest="method",
    default="mean_disp_plot",
    help="Method to choose top variable features. Choose one of : mean_disp_plot"
)

parser.add_argument(
    "-m", "--min-mean",
    type=float,
    action="store",
    dest="min_mean",
    default=None,
    help="Select genes with average mean greater than minimum mean."
)

parser.add_argument(
    "-M", "--max-mean",
    type=float,
    action="store",
    dest="max_mean",
    default=None,
    help="Select genes with average mean less than maximum mean."
)

parser.add_argument(
    "-d", "--min-dispersion",
    type=float,
    action="store",
    dest="min_disp",
    default=None,
    help="Select genes with dispersion greater than minimum dispersion."
)

parser.add_argument(
    "-D", "--max-dispersion",
    type=float,
    action="store",
    dest="max_disp",
    default=None,
    help="Select genes with dispersion less than maximum dispersion."
)

parser.add_argument(
    "-s", "--subset",
    type=str2bool,
    action="store",
    dest="subset",
    default=True,
    help="Subset the data with the highly variable features."
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

#
# Feature selection
#

if args.method == "mean_disp_plot":
    # Identify highly variable genes.
    # Expects logarithmized data: https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.api.pp.highly_variable_genes.html#scanpy.api.pp.highly_variable_genes
    sc.pp.highly_variable_genes(
        adata,
        min_mean=args.min_mean,
        max_mean=args.max_mean,
        min_disp=args.min_disp,
        max_disp=args.max_disp
    )
    if args.subset:
        # Keep only highly variable genes
        print("Subsetting highly variable features from the data...")
        adata = adata[:, adata.var['highly_variable']]
else:
    raise Exception("Method does not exist.")

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
import warnings

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
    default="wilcoxon",
    help="The default 'wilcoxon' uses Wilcoxon rank-sum, ‘t-test_overestim_var’ overestimates variance of each group,"
         " 't-test' uses t-test,  'logreg' uses logistic regression. See [Ntranos18], here and here, for why this is"
         " meaningful."
)
parser.add_argument(
    "-g", "--groupby",
    type=str,
    action="store",
    dest="groupby",
    default="louvain",
    help="The key of the observations grouping to consider."
)
parser.add_argument(
    "-n", "--ngenes",
    type=int,
    action="store",
    dest="ngenes",
    default=0,
    help="The number of genes that appear in the returned tables. Value of 0 will report all genes."
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
    raise Exception("Can only handle .h5ad files.")

if 'raw' not in dir(adata):
    warnings.warn(
        "There is no raw attribute in your anndata object. Differential analysis will be performed on the main (possibly normalised) matrix."
    )

if args.ngenes == 0:
    try:
        ngenes = adata.raw.shape[1]
    except AttributeError:
        ngenes = adata.shape[1]
        warnings.warn(
            "There is no raw attribute in your anndata object. Using the shape of the main matrix as the number of genes to report."
        )

sc.tl.rank_genes_groups(adata, groupby=args.groupby, method=args.method, n_genes=ngenes)

adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

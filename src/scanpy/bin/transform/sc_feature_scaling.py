#!/usr/bin/env python3

import os
from optparse import OptionParser

import scanpy as sc

parser = OptionParser(
    usage="usage: %prog [options] h5ad_file_path",
    version="%prog 1.0"
)
parser.add_option(
    "-x", "--method",
    type="string",
    action="store",
    dest="method",
    default="zscore_scale",
    help="Scale the data. Choose one of : zscore_scale"
)
parser.add_option(
    "-M", "--max-sd",
    type="float",
    action="store",
    dest="max_sd",
    default=-1,
    help="Clip values greater than maximum number of standard deviation."
)
(options, args) = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args[0]
FILE_PATH_OUT_BASENAME = os.path.splitext(args[1])[0]
SCALE__MAX_SD = options.max_sd if options.max_sd > 0 else None

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except IOError:
    raise Exception("Can only handle .h5ad files.")

#
# Transform the distribution of the data
#

if options.method == "zscore_scale":
    # scale each gene to unit variance, clip values exceeding SD max_sd.
    sc.pp.scale(adata, max_value=SCALE__MAX_SD)
else:
    raise Exception("Method does not exist.")

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

#!/usr/bin/env python
import os
from optparse import OptionParser

import scanpy as sc

parser = OptionParser(
    usage="usage: %prog [options] h5ad_file_path",
    version="%prog 1.0"
)
parser.add_option(
    "-x", "--method",
    action="store",
    dest="method",
    default="log1p",
    help="Transform the data."
)
(options, args) = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args[0]
FILE_PATH_OUT_BASENAME = os.path.splitext(args[1])[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

#
# Transform the distribution of the data
#

if options.method == "log1p":
    # log transform the data.
    sc.pp.log1p(adata)
else:
    raise Exception("Method does not exist.")

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

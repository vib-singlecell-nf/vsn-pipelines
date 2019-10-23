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
    type="string",
    action="store",
    dest="method",
    default="cpx",
    help="Normalize the data. Choose one of : cpx, regress_out"
)
parser.add_option(
    "-f", "--counts-per-cell-after",
    type="int",
    action="store",
    dest="counts_per_cell_after",
    default=1e4,
    help="Multiplying factor used when running 'cpx' method."
)
(options, args) = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args[0]
FILE_PATH_OUT_BASENAME = os.path.splitext(args[1])[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

#
# Normalize the data
#

adata.raw = adata

if options.method == 'cpx':
    # Total-count normalize (library-size correct) to '-r' reads/cell
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=options.counts_per_cell_after)
else:
    raise Exception("Method does not exist.")

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

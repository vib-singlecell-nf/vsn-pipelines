#!/usr/bin/env python
import os
from optparse import OptionParser
import scanpy as sc

parser = OptionParser(usage="usage: %prog [options] h5ad_file_path",
                      version="%prog 1.0")
parser.add_option("-x", "--method",
                  action="store",
                  dest="method",
                  default="linear_regression",
                  help="Normalize the data. Choose one of : regress_out")
parser.add_option("-r", "--variable-to-regress-out",
                  action="append",
                  dest="vars_to_regress_out",
                  default=None,
                  help="Variable to regress out. To regress multiple variables, add that many -v arguments. Used when running 'regress_out")
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
# Adjust the data
#

if options.method == 'linear_regression':
    # regress out variables (e.g.: total counts per cell and the percentage of mitochondrial genes expressed)
    sc.pp.regress_out(adata, options.vars_to_regress_out)
else:
    raise Exception("Method does not exist.")

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

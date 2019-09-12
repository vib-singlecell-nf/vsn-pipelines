#!/usr/bin/env python
import os
from optparse import OptionParser
import scanpy as sc

parser = OptionParser(usage="usage: %prog [options] h5ad_file_path",
                      version="%prog 1.0")
parser.add_option("-x", "--method",
                  type="string",
                  action="store",
                  dest="method",
                  default="mean_disp_plot",
                  help="Method to choose top variable features. Choose one of : mean_disp_plot")
parser.add_option("-m", "--min-mean",
                  type="float",
                  action="store",
                  dest="min_mean",
                  default=None,
                  help="Select genes with average mean greater than minimum mean.")
parser.add_option("-M", "--max-mean",
                  type="float",
                  action="store",
                  dest="max_mean",
                  default=None,
                  help="Select genes with average mean less than maximum mean.")
parser.add_option("-d", "--min-dispersion",
                  type="float",
                  action="store",
                  dest="min_disp",
                  default=None,
                  help="Select genes with dispersion greater than minimum dispersion.")
parser.add_option("-D", "--max-dispersion",
                  type="float",
                  action="store",
                  dest="max_disp",
                  default=None,
                  help="Select genes with dispersion less than maximum dispersion.")
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
# Feature selection
#

if options.method == "mean_disp_plot":
    # identify highly variable genes.
    # Expects logarithmized data: https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.api.pp.highly_variable_genes.html#scanpy.api.pp.highly_variable_genes
    sc.pp.highly_variable_genes(
        adata,
        min_mean=options.min_mean,
        max_mean=options.max_mean,
        min_disp=options.min_disp,
        max_disp=options.max_disp)
    # sc.pl.highly_variable_genes(adata)

    # keep only highly variable genes:
    # adata = adata[:, adata.var['highly_variable']]
else:
    raise Exception("Method does not exist.")

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

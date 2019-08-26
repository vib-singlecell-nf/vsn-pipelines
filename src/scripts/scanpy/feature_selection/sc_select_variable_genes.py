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
MEAN_DIS_POT__MIN_MEAN = options.min_mean if options.min_mean > 0 else None
MEAN_DIS_POT__MAX_MEAN = options.max_mean if options.max_mean > 0 else None
MEAN_DIS_POT__MIN_DISP = options.min_disp if options.min_disp > 0 else None
MEAN_DIS_POT__MAX_DISP = options.max_disp if options.max_disp > 0 else None

# I/O
# Expects h5ad file
try:
    adata=sc.read_h5ad(filename=FILE_PATH_IN)
except:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

#
# Feature selection
#

if options.method == "mean_disp_plot":
    # identify highly variable genes.
    sc.pp.highly_variable_genes(
        adata, 
        min_mean=MEAN_DIS_POT__MIN_MEAN, 
        max_mean=MEAN_DIS_POT__MAX_MEAN, 
        min_disp=MEAN_DIS_POT__MIN_DISP, 
        max_disp=MEAN_DIS_POT__MAX_DISP)
    # sc.pl.highly_variable_genes(adata)

    # keep only highly variable genes:
    # adata = adata[:, adata.var['highly_variable']]
else:
    raise Exception("Method does not exist.")

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

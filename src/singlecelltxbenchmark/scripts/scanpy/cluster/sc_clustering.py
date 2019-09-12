#!/usr/bin/env python
import os
from optparse import OptionParser
import scanpy as sc
import anndata as ad

from singlecelltxbenchmark.scripts.scanpy.utils import sc_neighbors

parser = OptionParser(usage="usage: %prog [options] h5ad_file_path",
                      version="%prog 1.0")
parser.add_option("-x", "--method",
                  type="string",
                  action="store",
                  dest="method",
                  default="Louvain",
                  help="Cluster cells using the Louvain algorithm [Blondel et al. (2008)] in the implementation of [Traag (2017)].")
parser.add_option("-r", "--resolution",
                  type="float",
                  action="store",
                  dest="resolution",
                  default=1.0,
                  help="[louvain], For the default flavor ('vtraag'), you can provide a resolution (higher resolution means finding more and smaller clusters) (Default: 1.0).")

# Add parser option for sc_neighbors
parser = sc_neighbors.add_options(parser=parser, method="Louvain")
(options, args) = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args[0]
FILE_PATH_OUT_BASENAME = os.path.splitext(args[1])[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except:
    raise Exception("Can only handle .h5ad files.")

#
# Transform the distribution of the data
#

if options.method == "Louvain":
    # Run Louvain clustering
    if "neighbors" not in adata.uns.keys():
        raise Exception("The neighborhood graph of observations has not been computed. Please do so before running {} clustering".format(options.method))
        # sc.pp.neighbors(adata, n_pcs=options.n_pcs, n_neighbors=options.n_neighbors)
    sc.tl.louvain(adata, resolution=options.resolution)
else:
    raise Exception("The dimensionality reduction method {} does not exist.".format(options.method))

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

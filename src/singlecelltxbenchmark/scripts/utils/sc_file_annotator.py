#!/usr/bin/env python
import os
from optparse import OptionParser
import scanpy as sc
import anndata as ad
import pandas as pd

parser = OptionParser(usage="usage: %prog [options] h5ad_file_path",
                      version="%prog 1.0")
parser.add_option("-t", "--type",
                  type="string",
                  action="store",
                  dest="type",
                  default="sample",
                  help="Are the entries of the meta data sample-based or cell-based? Choose one of: sample or cell")
parser.add_option("-m", "--meta-data-file-path",
                  type="string",
                  action="store",
                  dest="meta_data_file_path",
                  help="Path to the meta data. It expects a tabular separated file (.tsv) with header and a required 'id' column.")
(options, args) = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args[0]
FILE_PATH_OUT_BASENAME = os.path.splitext(args[1])[0]
SAMPLE_NAME = os.path.splitext(FILE_PATH_OUT_BASENAME)[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except:
    raise Exception("Can only handle .h5ad files.")

#
# Annotate the data
#

metadata = pd.read_csv(filepath_or_buffer=options.meta_data_file_path,
                       sep="\t")

if 'id' not in metadata.columns:
    raise Exception("The meta data TSV file expects a header with a required 'id' column.")

if options.type == "sample":
    for (column_name, column_data) in metadata[metadata.id == SAMPLE_NAME].iteritems():
        adata.obs[column_name] = column_data.values[0]
else:
    raise Exception("This meta data type {} is not implemented".format(options.type))

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

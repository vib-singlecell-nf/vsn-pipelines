#!/usr/bin/env python
import os
from optparse import OptionParser
import scanpy as sc

parser = OptionParser(usage="usage: %prog [options] h5ad_file_paths",
                      version="%prog 1.0")
parser.add_option("-f", "--file-format",
                  action="store",
                  dest="format",
                  default="h5ad",
                  help="Concatenate the data. Choose one of : h5ad")
parser.add_option("-j", "--join",
                  type="string",
                  action="store",
                  dest="join",
                  default="inner",
                  help="How to concatenate the multiple datasets. Choose one of : inner (intersect), outer (union).")
parser.add_option("-o", "--output",
                  action="store",
                  dest="output",
                  default=None,
                  help="Output file name.")
(options, args) = parser.parse_args()

# Define the arguments properly
FILE_PATH_OUT_BASENAME = os.path.splitext(options.output)[0]

# I/O
files = []

if options.format == 'h5ad':
    for FILE_PATH_IN in args:
        try:
            adata = sc.read_h5ad(filename=FILE_PATH_IN)
            files.append(adata)
        except:
            raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(FILE_PATH_IN))

#
# Adjust the data
#

if options.format == 'h5ad':
    # Concatenate multiple h5ad files
    # Source: https://anndata.readthedocs.io/en/latest/anndata.AnnData.concatenate.html#anndata.AnnData.concatenate
    adata = files[0].concatenate(files[1:], join=options.join)
else:
    raise Exception("Concatenation of .{} files is not implemented.".format(options.format))

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

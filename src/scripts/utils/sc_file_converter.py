#!/usr/bin/env python
import os
from optparse import OptionParser
import scanpy as sc

formats = ['10x_mtx', 'h5ad']

parser = OptionParser(usage="usage: %prog [options] datapath",
                        version="%prog 1.0")
parser.add_option("-i", "--input-format",
                    action="store",
                    dest="input_format",
                    default="",
                    help="Input format of the file to be converted. Choose one of: ".format(', '.join(formats)))
parser.add_option("-o", "--output-format",
                    action="store", # optional because action defaults to "store"
                    dest="output_format",
                    default="",
                    help="Output format which the file should be converted to. Choose one of: ".format(', '.join(formats)))
(options, args) = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args[0]
FILE_PATH_OUT_BASENAME = os.path.splitext(args[1])[0]
INPUT_FORMAT = options.input_format
OUTPUT_FORMAT = options.output_format

print(FILE_PATH_IN)

if INPUT_FORMAT == '10x_mtx' and OUTPUT_FORMAT == 'h5ad':
    # Sanity checks
    if not os.path.isdir(FILE_PATH_IN):
        raise Exception("Expecting a directory with an .mtx file when converting from 10x_mtx. {} does not seem to be one.".format(FILE_PATH_IN))
    if not os.path.exists(FILE_PATH_IN):
        raise Exception("The given directory {} does not exist.".format(FILE_PATH_IN))
    if not (not os.path.exists(os.path.join(FILE_PATH_IN, "matrix.mtx")) or not os.path.exists(os.path.join(FILE_PATH_IN, "matrix.mtx.gz"))):
        raise Exception("The given directory {} is not a proper 10xGenomics CellRanger folder. No .mtx[.gz] file found.".format(FILE_PATH_IN))    
    # Convert
    print("Reading 10x data...")
    adata = sc.read_10x_mtx(
        FILE_PATH_IN ,      # the directory with the `.mtx` file
        var_names='gene_symbols',   # use gene symbols for the variable names (variables-axis index)
        cache=False)
    print("Writing 10x data to h5ad...")
    adata.write_h5ad(filename="{}.h5ad".format(FILE_PATH_OUT_BASENAME))
else:
    raise Exception("File format conversion {0} --> {1} hasn't been implemented yet.".format(INPUT_FORMAT,OUTPUT_FORMAT))
print("Done!")
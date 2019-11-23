#!/usr/bin/env python
import os
from optparse import OptionParser
import pandas as pd

strand_options = {
    "no": 1,
    "forward": 2,
    "reverse": 3
}

parser = OptionParser(usage="usage: %prog [options] h5ad_file_paths",
                      version="%prog 1.0")
parser.add_option("-s", "--stranded",
                  action="store",
                  dest="stranded",
                  default="no",
                  help=f"Stranded nature of the library. Choose one of: {', '.join(strand_options.keys())}")
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

for FILE_PATH_IN in args:
    if not os.path.isfile(FILE_PATH_IN):
        raise Exception(f"Could not find file {FILE_PATH_IN}.")
    if not FILE_PATH_IN.endswith('ReadsPerGene.out.tab'):
        raise Exception(f"Expecting file ending with 'ReadsPerGene.out.tab', {os.path.basename(FILE_PATH_IN)} does not.")

    try:
        cell_name = os.path.basename(FILE_PATH_IN)[:-len("ReadsPerGene.out.tab")]
        counts = pd.read_csv(FILE_PATH_IN, sep='\t', index_col=0, skiprows=4, header=None)
        files.append((counts, cell_name))
    except IOError:
        raise Exception("Wrong input format. Expects .tab files, got .{}".format(FILE_PATH_IN))

#
# Adjust the data
#
try:
    all_counts = pd.DataFrame()
    for counts, cell_name in files:
        all_counts.loc[:, cell_name] = counts[strand_options[options.stranded]].astype(int)
except IOError:
    raise Exception("Concatenation failed.")

all_counts.to_csv(f"{FILE_PATH_OUT_BASENAME}.tsv", header=True, index=True, sep='\t')

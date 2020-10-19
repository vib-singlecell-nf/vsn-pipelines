#!/usr/bin/env python3

import argparse
import os
import pandas as pd

strand_options = {
    "no": 1,
    "forward": 2,
    "reverse": 3
}

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    nargs='+',
    type=argparse.FileType('r'),
    help='Input h5ad file.'
)

parser.add_argument(
    "-s", "--stranded",
    action="store",
    dest="stranded",
    default="no",
    help=f"Stranded nature of the library. Choose one of: {', '.join(strand_options.keys())}"
)

parser.add_argument(
    "-o", "--output",
    action="store",
    dest="output",
    default=None,
    help="Output file name."
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output)[0]

# I/O
files = []

for FILE_PATH_IN in args.input:
    FILE_PATH_IN = FILE_PATH_IN.name
    if not os.path.isfile(FILE_PATH_IN):
        raise Exception(f"Could not find file {FILE_PATH_IN}.")
    if not FILE_PATH_IN.endswith('ReadsPerGene.out.tab'):
        raise Exception(f"Expecting file ending with 'ReadsPerGene.out.tab', {os.path.basename(FILE_PATH_IN)} does not.")

    try:
        cell_name = os.path.basename(FILE_PATH_IN)[:-len("ReadsPerGene.out.tab")]
        counts = pd.read_csv(FILE_PATH_IN, sep='\t', index_col=0, skiprows=4, header=None)
        files.append((counts, cell_name))
    except IOError:
        raise Exception("VSN ERROR: Wrong input format. Expects .tab files, got .{}".format(FILE_PATH_IN))

#
# Adjust the data
#
try:
    all_counts = pd.DataFrame()
    for counts, cell_name in files:
        all_counts.loc[:, cell_name] = counts[strand_options[args.stranded]].astype(int)
except IOError:
    raise Exception("VSN ERROR: Concatenation failed.")

all_counts.to_csv(f"{FILE_PATH_OUT_BASENAME}.tsv", header=True, index=True, sep='\t')

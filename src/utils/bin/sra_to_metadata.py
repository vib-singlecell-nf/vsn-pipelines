#!/usr/bin/env python3

import os
import re
import argparse
from argparse import RawTextHelpFormatter
from pysradb import SRAdb
from pysradb.sraweb import SRAweb
import pandas as pd

parser = argparse.ArgumentParser(
    description='''
Convert a SRA ID to a metadata .tsv file with the following information
- experiment_accession, e.g.: SRX4084637
- experiment_title, e.g.: GSM3142622: w1118_1d_WholeBrain_Unstranded_RNA-seq; Drosophila melanogaster; RNA-Seq
- experiment_desc, e.g.: GSM3142622: w1118_1d_WholeBrain_Unstranded_RNA-seq; Drosophila melanogaster; RNA-Seq
- organism_taxid, e.g.: 7227
- organism_name, e.g.: Drosophila melanogaster
- library_strategy, e.g.: RNA-Seq
- library_source, e.g.: TRANSCRIPTOMIC
- library_selection, e.g.: cDNA
- sample_accession, e.g.: SRS3301695
- sample_title, e.g.: None
- study_accession, e.g.: SRP125768
- run_accession, e.g.: SRR7166639
- run_total_spots, e.g.: 3552575
- run_total_bases, e.g.: 176271295
- geo_accession, e.g.: GSM3142622
- sample_name, e.g.: w1118_1d_WholeBrain_Unstranded_RNA-seq

Example:
    python bin/sra_to_metadata.py \\
        -o test.tsv SRP125768 \\
        --sample-filter "DGRP-551_*d_r*" \\
        --sample-filter "w1118_*d_r*"

Using:
Zhu, Yuelin, Robert M. Stephens, Paul S. Meltzer, and Sean R. Davis. "SRAdb: query and use public next-generation sequencing data from within R." BMC bioinformatics 14, no. 1 (2013): 19.

''',
    formatter_class=RawTextHelpFormatter
)

parser.add_argument(
    "sra_project_id",
    type=str,
    help='The SRA project ID to get the metadata from.'
)

parser.add_argument(
    '-f', '--sample-filter',
    type=str,
    action="append",
    dest="sample_filters",
    help="A glob used as filter to keep only the matching samples."
)

parser.add_argument(
    "-d", "--sra-db",
    type=argparse.FileType('r'),
    required=False,
    dest="sra_db",
    help='The file path of the unzipped SQLite SRA database.'
)

parser.add_argument(
    "-o", "--output",
    type=argparse.FileType('w'),
    required=True,
    help='The .tsv file path that will stored the metadata for the given SRA Project ID.'
)

args = parser.parse_args()

#
# Get the metadata
#

if args.sra_db is not None:
    db = SRAdb(args.sra_db.name)
    print(f"Using local SRA SQLite database to query...")
else:
    print(f"Using NCBi's esearch and esummary interface to query...")
    db = SRAweb()

metadata = db.sra_metadata(
    args.sra_project_id,
    detailed=True,
    expand_sample_attributes=True,
    sample_attribute=True
)
# Drop any None columns
# pysradb does not lock the versions
# pandas 0.25.3 generates an additional None column compared to pandas 0.25.0
# Bug in 0.25.3 ?
metadata = metadata[metadata.columns.dropna()]

metadata = pd.concat(
    [
        metadata,
        metadata["experiment_title"].str.extract(
            r'^(.*): (.*); (.*); (.*)$', expand=True
        ).rename(
            columns={
                0: 'geo_accession',
                1: 'sample_name',
                2: 'species',
                3: 'library_type'
            }, inplace=False
        ).drop(
            ["species", "library_type"],
            axis=1,
            inplace=False
        )
    ],
    axis=1
)

# Filter the metadata based on the given ilters (if provided)
if args.sample_filters is not None:
    # Convert * (if not preceded by .) to .*
    def replace_bash_asterisk_wildcard(glob):
        return re.sub(r'(?<!\.)\*', '.*', glob)

    metadata = pd.concat(
        list(map(
            lambda glob: metadata[
                list(map(
                    lambda x: re.match(replace_bash_asterisk_wildcard(glob=glob), x) is not None,
                    metadata['sample_name'].values
                ))
            ].drop_duplicates(),
            args.sample_filters
        ))
    )

# I/O
metadata.to_csv(path_or_buf=args.output, index=False, header=True, sep="\t")

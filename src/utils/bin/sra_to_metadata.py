#!/usr/bin/env python3

import os
import argparse
from argparse import RawTextHelpFormatter
from pysradb.sraweb import SRAweb
import pandas as pd

parser = argparse.ArgumentParser(
    description='''
Convert a SRA ID to a meta data TSV file with the following information
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
- sample_title, e.g.: w1118_1d_WholeBrain_Unstranded_RNA-seq

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
    "-o", "--output",
    type=argparse.FileType('w'),
    required=True,
    help='The TSV file path that will stored the metadata for the given SRA Project ID.'
)

args = parser.parse_args()

#
# Get the metadata
#

db = SRAweb()
metadata = db.sra_metadata(srp=args.sra_project_id)
metadata = pd.concat(
    [
        metadata,
        metadata["experiment_title"].str.extract(
            r'^(.*): ([a-zA-Z0-9_-]*); (.*); (.*)$', expand=True
        ).rename(
            columns={
                0: 'geo_accession',
                1: 'sample_name',
                2: 'species',
                3: 'library_type'
            }, inplace=False
        ).drop(
            ["library_type", "species"], axis=1, inplace=False
        )
    ],
    axis=1
)

# I/O
metadata.to_csv(path_or_buf=args.output, index=False, header=True, sep="\t")

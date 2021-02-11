#!/usr/bin/env python3

import argparse
import pybiomart as pbm
import pandas as pd
import pickle

from pycisTopic.qc import compute_qc_stats

################################################################################

parser = argparse.ArgumentParser(description='Compute QC stats')

parser.add_argument(
    "--sampleId",
    type=str,
    required=True,
    help='Sample ID.'
)
parser.add_argument(
    "--fragments",
    type=str,
    required=True,
    help='Input fragments file.'
)
parser.add_argument(
    "--regions",
    type=str,
    required=True,
    help='Path to regions file.'
)
parser.add_argument(
    "--n_frag",
    type=int,
    required=True,
    default=50,
    help='Threshold on the number of fragments to keep for a barcode.'
)
parser.add_argument(
    "--output_metadata",
    type=str,
    help='Output file, tsv format.'
)
parser.add_argument(
    "--output_metadata_pkl",
    type=str,
    help='Metadata output, pickle format.'
)
parser.add_argument(
    "--output_profile_data_pkl",
    type=str,
    help='Profile data output, pickle format.'
)
parser.add_argument(
    "--threads",
    type=int,
    required=True,
    default=1,
    help='Number of threads to use.'
)

args = parser.parse_args()

################################################################################

dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://www.ensembl.org')
annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT', na=False)
annot = annot[~filter]
annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\S)', r'chr\1')
annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
annot = annot[annot.Transcript_type == 'protein_coding']

##################################################


fragments_dict = { 
    args.sampleId: args.fragments
    }
path_to_regions = {
    args.sampleId: args.regions
    }

metadata_bc_dict, profile_data_dict = compute_qc_stats(
                fragments_dict= fragments_dict,
                tss_annotation = annot,
                stats=['barcode_rank_plot', 'duplicate_rate', 'insert_size_distribution', 'profile_tss', 'frip'],
                label_list = None,
                path_to_regions = path_to_regions,
                n_cpu = args.threads,
                valid_bc = None,
                n_frag = args.n_frag,
                n_bc = None,
                tss_flank_window = 1000,
                tss_window = 50,
                tss_minimum_signal_window = 100,
                tss_rolling_window = 10,
                remove_duplicates = True)


### outputs:

metadata_bc_dict[args.sampleId].to_csv(args.output_metadata, sep='\t', index_label='barcode')

with open(args.output_metadata_pkl, 'wb') as f:
    pickle.dump(metadata_bc_dict, f)

with open(args.output_profile_data_pkl, 'wb') as f:
    pickle.dump(profile_data_dict, f)


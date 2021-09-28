#!/usr/bin/env python3

import argparse
import pickle

from pycisTopic.qc import plot_sample_metrics

################################################################################

parser = argparse.ArgumentParser(description='Compute QC stats')

parser.add_argument(
    "--sampleId",
    type=str,
    required=True,
    help='Sample ID.'
)
parser.add_argument(
    "--profile_data_pkl",
    type=str,
    help='Profile data, pickle format.'
)
parser.add_argument(
    "--output_pdf",
    type=str,
    help='Output plots, pdf format.'
)

args = parser.parse_args()

################################################################################

# Load sample metrics
infile = open(args.profile_data_pkl, 'rb')
profile_data_dict = pickle.load(infile)
infile.close()


# plot:
plot_sample_metrics(profile_data_dict,
           insert_size_distriubtion_xlim=[0,600],
           ncol=5,
           cmap='viridis',
           save=args.sampleId + "_qc_sample_metrics.pdf"
           )


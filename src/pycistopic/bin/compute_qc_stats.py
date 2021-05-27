#!/usr/bin/env python3

import argparse
import pybiomart as pbm
import pandas as pd
import pickle
import os

from pycisTopic.qc import compute_qc_stats

################################################################################

parser = argparse.ArgumentParser(description='Compute QC stats')

parser.add_argument(
    "--input_files",
    type=str,
    required=True,
    nargs='+',
    help='Input files in the form of [SampleId, path_to_fragments, path_to_peaks]. Multiple inputs are possible.'
)
parser.add_argument(
    "--sampleId",
    type=str,
    required=True,
    nargs='+',
    help='Sample ID.'
)
parser.add_argument(
    "--fragments",
    type=str,
    required=True,
    nargs='+',
    help='Input fragments file.'
)
parser.add_argument(
    "--regions",
    type=str,
    required=True,
    nargs='+',
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
    "--biomart_annot_pkl",
    type=str,
    help='Biomart annotations, pickle format.'
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

# Load biomart annotations:
infile = open(args.biomart_annot_pkl, 'rb')
annot = pickle.load(infile)
infile.close()


print(args.input_files)


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
                remove_duplicates = True,
                ### ray init args:
                include_dashboard=False,
                )


# load bap results to use for duplicate rate (if we are using bap output):
f_bap_qc = os.path.join(os.path.dirname(args.fragments),args.sampleId+'.QCstats.csv')
if os.path.isfile(f_bap_qc) and all(metadata_bc_dict[args.sampleId]['Dupl_rate_bap'] == 0):
    bapqc = pd.read_csv(f_bap_qc, index_col=0)
    metadata_bc_dict[args.sampleId]['Dupl_rate'] = bapqc['duplicateProportion']

### outputs:

metadata_bc_dict[args.sampleId].to_csv(args.output_metadata, sep='\t', index_label='barcode')

with open(args.output_metadata_pkl, 'wb') as f:
    pickle.dump(metadata_bc_dict, f)

with open(args.output_profile_data_pkl, 'wb') as f:
    pickle.dump(profile_data_dict, f)


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
    action='append',
    help='Input files in the form of "SampleId,fragments_filename,peaks_filename". Multiple inputs are possible.'
)
parser.add_argument(
    "--n_frag",
    type=int,
    required=True,
    default=50,
    help='Threshold on the number of fragments to keep for a barcode.'
)
parser.add_argument(
    "--tss_flank_window",
    type=int,
    required=True,
    default=2000,
    help='Flanking window around the TSS.'
)
parser.add_argument(
    "--tss_window",
    type=int,
    required=True,
    default=50,
    help='Window around the TSS used to count fragments in the TSS when calculating the TSS enrichment per barcode.'
)
parser.add_argument(
    "--tss_minimum_signal_window",
    type=int,
    required=True,
    default=100,
    help='Tail window use to normalize the TSS enrichment (average signal in the X bp in the extremes of the TSS window).'
)
parser.add_argument(
    "--tss_rolling_window",
    type=int,
    required=True,
    default=10,
    help='Rolling window used to smooth signal.'
)
parser.add_argument(
    "--min_norm",
    type=float,
    required=True,
    default=0.1,
    help='Minimum normalization score. If the average minimum signal value is below this value, this number is used to normalize the TSS signal. This approach penalizes cells with fewer reads.'
)

parser.add_argument(
    "--biomart_annot_pkl",
    type=str,
    help='Biomart annotations, pickle format.'
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


fragments_dict = { x[0].split(',')[0]: x[0].split(',')[1] for x in args.input_files }
path_to_regions = { x[0].split(',')[0]: x[0].split(',')[2] for x in args.input_files }

# Load biomart annotations:
infile = open(args.biomart_annot_pkl, 'rb')
annot = pickle.load(infile)
infile.close()



metadata_bc_dict, profile_data_dict = compute_qc_stats(
                fragments_dict=fragments_dict,
                tss_annotation=annot,
                stats=['barcode_rank_plot', 'duplicate_rate', 'insert_size_distribution', 'profile_tss', 'frip'],
                label_list=None,
                path_to_regions=path_to_regions,
                n_cpu=args.threads,
                valid_bc=None,
                n_frag=args.n_frag,
                n_bc=None,
                tss_flank_window=args.tss_flank_window,
                tss_window=args.tss_window,
                tss_minimum_signal_window=args.tss_minimum_signal_window,
                tss_rolling_window=args.tss_rolling_window,
                min_norm=args.min_norm,
                remove_duplicates = True,
                #_temp_dir=
                )


## load bap results to use for duplicate rate (if we are using bap output):
#f_bap_qc = os.path.join(os.path.dirname(args.fragments),args.sampleId+'.QCstats.csv')
#if os.path.isfile(f_bap_qc) and all(metadata_bc_dict[args.sampleId]['Dupl_rate_bap'] == 0):
#    bapqc = pd.read_csv(f_bap_qc, index_col=0)
#    metadata_bc_dict[args.sampleId]['Dupl_rate'] = bapqc['duplicateProportion']

### outputs:

with open(args.output_metadata_pkl, 'wb') as f:
    pickle.dump(metadata_bc_dict, f)

with open(args.output_profile_data_pkl, 'wb') as f:
    pickle.dump(profile_data_dict, f)


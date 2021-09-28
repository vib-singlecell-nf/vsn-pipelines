#!/usr/bin/env python3

import argparse
import pandas as pd
import pickle
import numpy as np

from pycisTopic.qc import plot_barcode_metrics

################################################################################

parser = argparse.ArgumentParser(description='Barcode level statistics')

parser.add_argument(
    "--sampleId",
    type=str,
    required=True,
    help='Sample ID.'
)
parser.add_argument(
    "--metadata_pkl",
    type=str,
    help='Metadata, pickle format.'
)
parser.add_argument(
    "--selected_barcodes",
    type=str,
    help='Output file containing selected barcodes.'
)

### fragments filters:
parser.add_argument(
    "--filter_frags_lower",
    type=float,
    required=False,
    default=3,
    help='Lower threshold on the number of fragments for keeping a barcode. Log10-scaled.'
)
parser.add_argument(
    "--filter_frags_upper",
    type=float,
    required=False,
    default=None,
    help='Upper threshold on the number of fragments for keeping a barcode. Log10-scaled.'
)

### TSS Encrichment filters:
parser.add_argument(
    "--filter_tss_lower",
    type=float,
    required=False,
    default=8,
    help='Lower threshold on the TSS Enrichment for keeping a barcode.'
)
parser.add_argument(
    "--filter_tss_upper",
    type=float,
    required=False,
    default=None,
    help='Upper threshold on the TSS Enrichment for keeping a barcode.'
)

### FRIP filters:
parser.add_argument(
    "--filter_frip_lower",
    type=float,
    required=False,
    default=None,
    help='Lower threshold on FRIP for keeping a barcode.'
)
parser.add_argument(
    "--filter_frip_upper",
    type=float,
    required=False,
    default=None,
    help='Upper threshold on FRIP for keeping a barcode.'
)

### Duplication rate filters:
parser.add_argument(
    "--filter_dup_rate_lower",
    type=float,
    required=False,
    default=None,
    help='Lower threshold on duplication rate for keeping a barcode.'
)
parser.add_argument(
    "--filter_dup_rate_upper",
    type=float,
    required=False,
    default=None,
    help='Upper threshold on duplication rate for keeping a barcode.'
)

args = parser.parse_args()

################################################################################

# the pycisTopic filter setting is in log10 scale
if args.filter_frags_lower is not None:
    args.filter_frags_lower = np.log10(args.filter_frags_lower)
if args.filter_frags_upper is not None:
    args.filter_frags_upper = np.log10(args.filter_frags_upper)


# Load barcode metrics
infile = open(args.metadata_pkl, 'rb')
metadata_bc = pickle.load(infile)
infile.close()


# Return figure to plot together with other metrics, and cells passing filters. Figure will be saved as pdf.
FRIP_NR_FRAG_fig, FRIP_NR_FRAG_filter=plot_barcode_metrics(metadata_bc[args.sampleId],
                                       var_x='Log_unique_nr_frag',
                                       var_y='FRIP',
                                       min_x=args.filter_frags_lower,
                                       max_x=args.filter_frags_upper,
                                       min_y=args.filter_frags_lower,
                                       max_y=args.filter_frags_upper,
                                       cmap='viridis',
                                       return_cells=True,
                                       return_fig=True,
                                       plot=False,
                                       save='./'+args.sampleId+'__FRIP-vs-nFrag.pdf')

# Return figure to plot together with other metrics, and cells passing filters
TSS_NR_FRAG_fig, TSS_NR_FRAG_filter=plot_barcode_metrics(metadata_bc[args.sampleId],
                                      var_x='Log_unique_nr_frag',
                                      var_y='TSS_enrichment',
                                      min_x=args.filter_frags_lower,
                                      max_x=args.filter_frags_upper,
                                      min_y=args.filter_tss_lower,
                                      max_y=args.filter_tss_upper,
                                      cmap='viridis',
                                      return_cells=True,
                                      return_fig=True,
                                      plot=False,
                                      save='./'+args.sampleId+'__TSS-vs-nFrag.pdf')

# Return figure to plot together with other metrics, but not returning cells (no filter applied for the duplication rate  per barcode)
DR_NR_FRAG_fig=plot_barcode_metrics(metadata_bc[args.sampleId],
                                      var_x='Log_unique_nr_frag',
                                      var_y='Dupl_rate',
                                      min_x=args.filter_frags_lower,
                                      max_x=args.filter_frags_upper,
                                      min_y=args.filter_dup_rate_lower,
                                      max_y=args.filter_dup_rate_upper,
                                      cmap='viridis',
                                      return_cells=False,
                                      return_fig=True,
                                      plot=False,
                                      save='./'+args.sampleId+'__duprate-vs-nFrag.pdf')

# intersection of barcodes to keep:
bc_passing_filters = list(set(FRIP_NR_FRAG_filter) & set(TSS_NR_FRAG_filter))


### outputs:

pd.DataFrame(bc_passing_filters).to_csv(args.selected_barcodes, sep='\t', index=False, header=False)


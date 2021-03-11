#!/usr/bin/env python3

import argparse
import pandas as pd
import pickle
import numpy as np

from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt

################################################################################

parser = argparse.ArgumentParser(description='Call cells from barcodes')

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

### fragments filters:
parser.add_argument(
    "--filter_frags_lower",
    type=float,
    required=False,
    default=3,
    help='Lower threshold on the number of fragments for keeping a barcode.'
)
parser.add_argument(
    "--filter_frags_upper",
    type=float,
    required=False,
    default=None,
    help='Upper threshold on the number of fragments for keeping a barcode.'
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


def plot_frag_qc(x, y,
                 ax,
                 x_thr_min=None, x_thr_max=None,
                 y_thr_min=None, y_thr_max=None,
                 ylab=None,
                 xlab="Number of (unique) fragments",
                 cmap='viridis',
                 density_overlay=False,
                 s=10,
                 marker='+',
                 c='#343434',
                 xlim=None,
                 ylim=None,
                 **kwargs
                ):
    assert all(x.index == y.index)
    barcodes = x.index.values
    if density_overlay:
        #pdf,axes = fastKDE.pdf(x.to_numpy(),y.to_numpy())
        xy = np.vstack([np.log(x),y])
        z = gaussian_kde(xy)(xy)
        idx = z.argsort()
        x, y, z, barcodes = x[idx], y[idx], z[idx], barcodes[idx]
    else:
        z=c

    barcodes_to_keep=[]
    sp=ax.scatter(x, y, c=z, s=s, edgecolors=None, marker=marker, cmap=cmap, **kwargs)
    #fig.colorbar(sp)
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
    if xlim is not None:
        ax.set_ylim(xlim[0], xlim[1])
    # thresholds:
    if x_thr_min is not None:
        ax.axvline(x=x_thr_min, color='r', linestyle='--')
        barcodes_to_keep.append(barcodes[x>x_thr_min])
    if x_thr_max is not None:
        ax.axvline(x=x_thr_max, color='r', linestyle='--')
        barcodes_to_keep.append(barcodes[x<x_thr_max])
    if y_thr_min is not None:
        ax.axhline(y=y_thr_min, color='r', linestyle='--')
        barcodes_to_keep.append(barcodes[y>y_thr_min])
    if y_thr_max is not None:
        ax.axhline(y=y_thr_max, color='r', linestyle='--')
        barcodes_to_keep.append(barcodes[y<y_thr_max])
    ax.set_xscale("log")
    ax.set_xmargin(0.01)
    ax.set_ymargin(0.01)
    ax.set_xlabel(xlab,fontsize=10)
    ax.set_ylabel(ylab,fontsize=10)

    if len(barcodes_to_keep)>0:
        return list(set.intersection(*map(set, barcodes_to_keep)))
    else:
        return barcodes


# Load barcode metrics
infile = open(args.metadata_pkl, 'rb')
metadata_bc_dict = pickle.load(infile)
infile.close()

fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(15,5), dpi=150 )
p1_cells = plot_frag_qc(
    x = metadata_bc_dict[args.sampleId]['Unique_nr_frag'],
    y = metadata_bc_dict[args.sampleId]['TSS_enrichment'],
    ylab = 'TSS Enrichment',
    x_thr_min=args.filter_frags_lower,
    x_thr_max=args.filter_frags_upper,
    y_thr_min=args.filter_tss_lower,
    y_thr_max=args.filter_tss_upper,
    density_overlay=True,
    ax=ax1
)
p2_cells = plot_frag_qc(
    x = metadata_bc_dict[args.sampleId]['Unique_nr_frag'],
    y = metadata_bc_dict[args.sampleId]['FRIP'],
    ylab = 'FRIP',
    ylim=[0,1],
    x_thr_min=args.filter_frags_lower,
    x_thr_max=args.filter_frags_upper,
    y_thr_min=args.filter_frip_lower,
    y_thr_max=args.filter_frip_upper,
    density_overlay=True,
    ax=ax2
)
p3_cells = plot_frag_qc(
    x = metadata_bc_dict[args.sampleId]['Unique_nr_frag'],
    y = metadata_bc_dict[args.sampleId]['Dupl_rate'],
    ylab = 'Duplicate rate per cell',
    ylim=[0,1],
    x_thr_min=args.filter_frags_lower,
    x_thr_max=args.filter_frags_upper,
    y_thr_min=args.filter_dup_rate_lower,
    y_thr_max=args.filter_dup_rate_upper,
    density_overlay=True,
    ax=ax3
)
fig.suptitle(args.sampleId)
plt.tight_layout()
plt.savefig(args.sampleId + '__fragments_qc.pdf', dpi=300, bbox_inches = 'tight')

# intersection of barcodes to keep:
bc_passing_filters = list(set(p1_cells) & set(p2_cells) & set(p3_cells))
metadata_bc_dict[args.sampleId]['Keep'] = [ 1 if x in bc_passing_filters else 0 for x in metadata_bc_dict[args.sampleId].index ]


### outputs:
pd.DataFrame(bc_passing_filters).to_csv(args.sampleId + '__selected_barcodes.txt', sep='\t', index=False, header=False)

with open(args.sampleId + '__metadata_with_calls.pickle', 'wb') as f:
    pickle.dump(metadata_bc_dict, f)


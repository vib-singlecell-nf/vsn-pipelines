#!/usr/bin/env python3

import argparse
import os
import numpy as np
import scrublet as scr
import scanpy as sc


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


parser = argparse.ArgumentParser(description='Template script')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input h5ad file containing the raw data'
)

parser.add_argument(
    "input_hvg",
    type=argparse.FileType('w'),
    help='Input h5ad containing the highly_variable slot.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output txt file containing the predicted doublet by Scrublet.'
)

parser.add_argument(
    "-s", "--synthetic-doublet-umi-subsampling",
    type=float,
    dest="synthetic_doublet_umi_subsampling",
    default=1.0,
    help='Rate for sampling UMIs when creating synthetic doublets.'
)

parser.add_argument(
    "-m", "--min-counts",
    type=int,
    dest="min_counts",
    default=3,
    help='Used for gene filtering prior to PCA. Genes expressed at fewer than `min_counts` in fewer than `min_cells` (see below) are excluded.'
)

parser.add_argument(
    "-n", "--min-cells",
    type=int,
    dest="min_cells",
    default=3,
    help='Used for gene filtering prior to PCA. Genes expressed at fewer than `min_counts` (see above) in fewer than `min_cells` are excluded.'
)

parser.add_argument(
    "-v", "--min-gene-variability-pctl",
    type=float,
    dest="min_gene_variability_pctl",
    default=0.85,
    help='''
        Used for gene filtering prior to PCA. Keep the most highly variable genes (in the top min_gene_variability_pctl percentile),
        as measured by the v-statistic [Klein et al., Cell 2015].
        '''
)

parser.add_argument(
    "-l", "--log-transform",
    type=str2bool,
    action="store",
    dest="log_transform",
    default=False,
    help='''
        If True, log-transform the counts matrix (log10(1+TPM)).
        `sklearn.decomposition.TruncatedSVD` will be used for dimensionality reduction, unless `mean_center` is True.
        '''
)

parser.add_argument(
    "-c", "--mean-center",
    type=str2bool,
    action="store",
    dest="mean_center",
    default=True,
    help='''
        If True, center the data such that each gene has a mean of 0.  
        `sklearn.decomposition.PCA` will be used for dimensionality
        '''
)

parser.add_argument(
    "-w", "--normalize-variance",
    type=str2bool,
    action="store",
    dest="normalize_variance",
    default=True,
    help='''
        If True, normalize the data such that each gene has a variance of 1.
        `sklearn.decomposition.TruncatedSVD` will be used for dimensionality
        reduction, unless `mean_center` is True.
        '''
)

parser.add_argument(
    "-p", "--n-prin-comps",
    type=int,
    dest="n_prin_comps",
    default=30,
    help='Number of principal components used to embed the transcriptomes prior to k-nearest-neighbor graph construction.'
)

parser.add_argument(
    "-t", "--technology",
    type=str,
    dest="technology",
    choices=["10xv2"],
    help='Single-cell technology used.'
)

args = parser.parse_args()


# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_IN_HVG = args.input_hvg
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

# I/O
# Expects h5ad file
try:
    adata_raw = sc.read_h5ad(filename=FILE_PATH_IN.name)
    hvg_adata = sc.read_h5ad(filename=FILE_PATH_IN_HVG.name)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

################################################################################
# Processing...

adata_raw_var = adata_raw.X[:, np.array(hvg_adata.var['highly_variable'])]
scrub = scr.Scrublet(adata_raw_var)
adata_raw.obs['doublet_scores'], adata_raw.obs['predicted_doublets'] = scrub.scrub_doublets(
    synthetic_doublet_umi_subsampling=args.synthetic_doublet_umi_subsampling,
    use_approx_neighbors=True,
    distance_metric='euclidean',
    get_doublet_neighbor_parents=False,
    min_counts=args.min_counts,
    min_cells=args.min_cells,
    min_gene_variability_pctl=args.min_gene_variability_pctl,
    log_transform=args.log_transform,
    mean_center=args.mean_center,
    normalize_variance=args.normalize_variance,
    n_prin_comps=args.n_prin_comps,
    verbose=True
)

if args.technology == "10xv2":
    # Take doublet cells based on expected doublets based on number of cells (10x Chromium)
    cells_recovered = len(adata_raw)
    doublet_rate = 0.0008 * cells_recovered + 0.0527
    expected_doublets = np.int(doublet_rate / 100 * cells_recovered)
    doublet_cells = adata_raw.obs['doublet_scores'].sort_values(
        ascending=False
    ).head(
        n=expected_doublets
    ).index
else:
    raise Exception(f"Doublet detection with Scrublet for the given technolog {args.technology} is not implemented")

################################################################################

# I/O
np.savetxt(f"{FILE_PATH_OUT_BASENAME}.txt", doublet_cells, fmt="%s")

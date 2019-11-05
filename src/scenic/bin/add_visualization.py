#!/usr/bin/env python3

import argparse
import base64
import json
import zlib
from multiprocessing import cpu_count
from shutil import copyfile

import loompy as lp
import numpy as np
import pandas as pd

import umap
from MulticoreTSNE import MulticoreTSNE as TSNE

import export_to_loom

################################################################################
################################################################################

parser = argparse.ArgumentParser(description='Create t-SNE and UMAP for SCENIC AUCell matrix')
parser.add_argument(
    '--loom_input',
    help='Loom file containing AUCell results',
    required=True,
    default='input.loom'
)
parser.add_argument(
    '--loom_output',
    help='Final loom file with dimensionality reductions embedded',
    required=True,
    default='output.loom'
)
parser.add_argument(
    '--num_workers',
    type=int,
    default=(cpu_count() - 1),
    help='The number of workers to use. (default: {}).'.format(cpu_count() - 1)
)
args = parser.parse_args()


def visualize_AUCell(args):

    ################################################################################
    # load data from loom
    ################################################################################

    with lp.connect(args.loom_input, mode='r', validate=False) as lf:
        auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)

    ################################################################################
    # Visualize AUC matrix:
    ################################################################################

    # UMAP
    run_umap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
    dr_umap = run_umap(auc_mtx.dropna())

    # tSNE
    tsne = TSNE(n_jobs=args.num_workers)
    dr_tsne = tsne.fit_transform(auc_mtx.dropna())

    ################################################################################
    # update scenic data
    ################################################################################

    scope_loom = export_to_loom.SCopeLoom.read_loom(filename=args.loom_input)
    scope_loom.add_embedding(embedding=dr_umap, embedding_name="SCENIC AUC UMAP", is_default=True)
    scope_loom.add_embedding(embedding=dr_tsne, embedding_name="SCENIC AUC t-SNE", is_default=False)
    scope_loom.export(out_fname=args.loom_output)


if __name__ == "__main__":
    visualize_AUCell(args)
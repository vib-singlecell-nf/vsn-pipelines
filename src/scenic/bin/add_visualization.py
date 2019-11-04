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


def df_to_named_matrix(df):
    arr_ip = [tuple(i) for i in df.as_matrix()]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr


def visualize_AUCell(args):

    ################################################################################
    # load data from loom
    ################################################################################

    # scenic motif output
    lf = lp.connect(args.loom_input, mode='r', validate=False)
    meta = json.loads(zlib.decompress(base64.b64decode(lf.attrs.MetaData)))
    auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
    regulons = pd.DataFrame(lf.ra.Regulons, index=lf.ra.Gene)
    lf.close()


    ################################################################################
    # Fix regulon objects to display properly in SCope:
    ################################################################################

    # add underscore for SCope compatibility:
    auc_mtx.columns = auc_mtx.columns.str.replace('\\(', '_(')

    # add underscore for SCope compatibility:
    regulons.columns = regulons.columns.str.replace('\\(', '_(')

    # Rename regulons in the thresholds object, motif
    rt = meta['regulonThresholds']
    for i, x in enumerate(rt):
        tmp = x.get('regulon').replace("(", "_(") # + '-motif'
        x.update({'regulon': tmp})


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
    # embeddings
    ################################################################################

    default_embedding = pd.DataFrame(dr_umap, columns=['_X', '_Y'], index=auc_mtx.dropna().index)

    embeddings_x = pd.DataFrame(dr_tsne, columns=['_X', '_Y'], index=auc_mtx.dropna().index)[['_X']].astype('float32')
    embeddings_y = pd.DataFrame(dr_tsne, columns=['_X', '_Y'], index=auc_mtx.dropna().index)[['_Y']].astype('float32')

    embeddings_x.columns = ['1']
    embeddings_y.columns = ['1']


    ################################################################################
    # copy loom
    ################################################################################

    copyfile(args.loom_input, args.loom_output)


    ################################################################################
    # update scenic data
    ################################################################################

    lf = lp.connect(args.loom_output, mode='r+', validate=False)

    # write regulon information:
    lf.ca['RegulonsAUC'] = df_to_named_matrix(auc_mtx)
    lf.ra['Regulons'] = df_to_named_matrix(regulons)

    # write embeddings:
    lf.ca['Embedding']    = df_to_named_matrix(default_embedding)
    lf.ca['Embeddings_X'] = df_to_named_matrix(embeddings_x)
    lf.ca['Embeddings_Y'] = df_to_named_matrix(embeddings_y)

    metaJson = {}
    metaJson['embeddings'] = [
        {
            "id": -1,
            "name": "SCENIC AUC UMAP"
        },
        {
            "id": 1,
            "name": "SCENIC AUC t-SNE"
        },
    ]

    lf.attrs['MetaData'] = json.dumps(metaJson)

    lf.close()


if __name__ == "__main__":
    visualize_AUCell(args)


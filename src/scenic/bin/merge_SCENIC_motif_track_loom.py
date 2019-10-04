#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from multiprocessing import cpu_count
from MulticoreTSNE import MulticoreTSNE as TSNE
import umap
import json
import zlib
import base64
import argparse

################################################################################
################################################################################

parser = argparse.ArgumentParser(description='Integrate output from pySCENIC motif- and track-based runs')
parser.add_argument('--loom_motif', help='Loom file from pySCENIC motif run', required=True, default='pyscenic_motif.loom')
parser.add_argument('--loom_track', help='Loom file from pySCENIC track run', required=True, default='pyscenic_motif.loom')
parser.add_argument('--loom_output', help='Final loom file with pySCENIC motif and track results integrated', required=True, default='pyscenic.loom')
parser.add_argument('--num_workers',
                    type=int, default=(cpu_count() - 1),
                    help='The number of workers to use. (default: {}).'.format(cpu_count() - 1))
args = parser.parse_args()

################################################################################
################################################################################

# args.loom_motif = "/ddn1/vol1/staging/leuven/stg_00002/lcb/cflerin/testruns/SingleCellTxBenchmark/src/singlecelltxbenchmark/pipelines/scenic/work/6d/50e6af1508cb58f88d1800bee20720/auc_mtf.loom"
# args.loom_track = "/ddn1/vol1/staging/leuven/stg_00002/lcb/cflerin/testruns/SingleCellTxBenchmark/src/singlecelltxbenchmark/pipelines/scenic/work/73/09b0e04f95b803b1a1393e19944df0/auc_trk.loom"
# args.loom_output = "/ddn1/vol1/staging/leuven/stg_00002/lcb/cflerin/testruns/SingleCellTxBenchmark/src/singlecelltxbenchmark/pipelines/scenic/work/testout.loom"


def dfToNamedMatrix(df):
    arr_ip = [tuple(i) for i in df.as_matrix()]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr


def integrateMotifTrack(args):
    ################################################################################
    # load data from loom
    ################################################################################

    # scenic motif output
    # scenic output
    lf = lp.connect(args.loom_motif, mode='r', validate=False)
    meta_mtf = json.loads(zlib.decompress(base64.b64decode(lf.attrs.MetaData)))
    auc_mtx_mtf = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
    regulons_mtf = pd.DataFrame(lf.ra.Regulons, index=lf.ra.Gene)
    lf.close()

    # scenic track output
    lf = lp.connect(args.loom_track, mode='r', validate=False)
    meta_trk = json.loads(zlib.decompress(base64.b64decode(lf.attrs.MetaData)))
    auc_mtx_trk = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
    regulons_trk = pd.DataFrame(lf.ra.Regulons, index=lf.ra.Gene)
    lf.close()

    ################################################################################
    #
    ################################################################################

    ################################################################################
    # Fix track auc mtx names:
    ################################################################################

    # auc_mtx_mtf.columns = auc_mtx_mtf.columns.str.replace('\(','_(')
    # auc_mtx_trk.columns = auc_mtx_trk.columns.str.replace('\(','_(')

    # relabel columns with suffix indicating the regulon source
    auc_mtx_trk.columns = auc_mtx_trk.columns + '-track'
    auc_mtx_mtf.columns = auc_mtx_mtf.columns + '-motif'

    # merge the AUC matrices:
    auc_mtx = pd.concat([auc_mtx_mtf, auc_mtx_trk], sort=False, axis=1, join='outer')

    # fill NAs (if any) with 0s:
    auc_mtx.fillna(0, inplace=True)

    # add underscore for SCope compatibility:
    auc_mtx.columns = auc_mtx.columns.str.replace('\(', '_(')

    ################################################################################
    # Fix regulon objects to display properly in SCope:
    ################################################################################

    # combine regulon assignment matrices:
    regulons_trk.columns = regulons_trk.columns + '-track'
    regulons_mtf.columns = regulons_mtf.columns + '-motif'

    # merge the regulon assignment matrices:
    regulons = pd.concat([regulons_mtf, regulons_trk], sort=False, axis=1, join='outer')

    # replace NAs with 0s:
    regulons.fillna(0, inplace=True)

    # add underscore for SCope compatibility:
    regulons.columns = regulons.columns.str.replace('\(', '_(')

    # Rename regulons in the thresholds object, motif
    rt_mtf = meta_mtf['regulonThresholds']
    for i, x in enumerate(rt_mtf):
        tmp = x.get('regulon').replace("(", "_(") + '-motif'
        x.update({'regulon': tmp})

    # Rename regulons in the thresholds object, track
    rt_trk = meta_trk['regulonThresholds']
    for i, x in enumerate(rt_trk):
        tmp = x.get('regulon').replace("(", "_(") + '-track'
        x.update({'regulon': tmp})
        # blank out the "motifData" field for track-based regulons:
        x.update({'mofitData': 'NA.png'})

    # merge regulon threshold dictionaries:
    rt = rt_mtf + rt_trk

    ################################################################################
    # Visualize AUC matrix:
    ################################################################################

    # UMAP
    runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
    dr_umap = runUmap(auc_mtx.dropna())
    # pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.dropna().index).to_csv("scenic_motif-track_umap.txt", sep='\t')
    # tSNE
    tsne = TSNE(n_jobs=args.num_workers)
    dr_tsne = tsne.fit_transform(auc_mtx.dropna())

    ################################################################################
    # embeddings
    ################################################################################

    defaultEmbedding = pd.DataFrame(dr_umap, columns=['_X', '_Y'], index=auc_mtx.dropna().index)

    Embeddings_X = pd.DataFrame(dr_tsne, columns=['_X', '_Y'], index=auc_mtx.dropna().index).iloc[:, 0].astype('float32')
    Embeddings_Y = pd.DataFrame(dr_tsne, columns=['_X', '_Y'], index=auc_mtx.dropna().index).iloc[:, 1].astype('float32')

    Embeddings_X.columns = ['1']
    Embeddings_Y.columns = ['1']

    ################################################################################
    # metadata
    ################################################################################

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

    metaJson["metrics"] = [
        {
            "name": "nUMI"
        }, {
            "name": "nGene"
        }, {
            "name": "Percent_mito"
        }
    ]

    metaJson["annotations"] = [
    ]

    # SCENIC regulon thresholds:
    metaJson["regulonThresholds"] = rt

    ################################################################################
    # row/column attributes
    ################################################################################

    # re-open the connection to the loom file to copy the original expression data
    lf = lp.connect(args.loom_motif, mode='r', validate=False)
    # exprMat = pd.DataFrame(lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
    # lf.close()

    col_attrs = {
        "CellID": lf.ca.CellID,  # np.array(adata.obs.index),
        "Embedding": dfToNamedMatrix(defaultEmbedding),
        # "Embeddings_X": dfToNamedMatrix(Embeddings_X),
        # "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
        "RegulonsAUC": dfToNamedMatrix(auc_mtx),
    }

    row_attrs = {
        "Gene": lf.ra.Gene,
        "Regulons": dfToNamedMatrix(regulons),
    }

    attrs = {
        "MetaData": json.dumps(metaJson),
    }

    attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')

    lp.create(
        filename=args.loom_output,
        layers=lf[:, :],
        row_attrs=row_attrs,
        col_attrs=col_attrs,
        file_attrs=attrs
    )
    lf.close()


if __name__ == "__main__":
    integrateMotifTrack(args)

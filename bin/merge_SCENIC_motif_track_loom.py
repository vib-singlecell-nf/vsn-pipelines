#!/usr/bin/env python3

import argparse
import base64
import json
import zlib
from multiprocessing import cpu_count

import loompy as lp
import numpy as np
import pandas as pd

################################################################################
################################################################################

parser = argparse.ArgumentParser(description='Integrate output from pySCENIC motif- and track-based runs')
parser.add_argument(
    '--loom_motif',
    help='Loom file from pySCENIC motif run',
    required=True,
    default='pyscenic_motif.loom'
)
parser.add_argument(
    '--loom_track',
    help='Loom file from pySCENIC track run',
    required=True,
    default='pyscenic_track.loom'
)
parser.add_argument(
    '--loom_output',
    help='Final loom file with pySCENIC motif and track results integrated',
    required=True,
    default='pyscenic.loom'
)
args = parser.parse_args()


################################################################################
################################################################################


def df_to_named_matrix(df):
    arr_ip = [tuple(i) for i in df.as_matrix()]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr


def integrate_motif_track(args):
    ################################################################################
    # load data from loom
    ################################################################################

    # scenic motif output
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
    # Fix track auc mtx names:
    ################################################################################

    # relabel columns with suffix indicating the regulon source
    auc_mtx_trk.columns = auc_mtx_trk.columns + '-track'
    auc_mtx_mtf.columns = auc_mtx_mtf.columns + '-motif'

    # merge the AUC matrices:
    auc_mtx = pd.concat([auc_mtx_mtf, auc_mtx_trk], sort=False, axis=1, join='outer')

    # fill NAs (if any) with 0s:
    auc_mtx.fillna(0, inplace=True)

    ################################################################################
    # Combine track and motif
    ################################################################################

    # combine regulon assignment matrices:
    regulons_trk.columns = regulons_trk.columns + '-track'
    regulons_mtf.columns = regulons_mtf.columns + '-motif'

    # merge the regulon assignment matrices:
    regulons = pd.concat([regulons_mtf, regulons_trk], sort=False, axis=1, join='outer')

    # replace NAs with 0s:
    regulons.fillna(0, inplace=True)

    # Rename regulons in the thresholds object, motif
    rt_mtf = meta_mtf['regulonThresholds']
    for i, x in enumerate(rt_mtf):
        tmp = x.get('regulon') + '-motif'
        x.update({'regulon': tmp})

    # Rename regulons in the thresholds object, track
    rt_trk = meta_trk['regulonThresholds']
    for i, x in enumerate(rt_trk):
        tmp = x.get('regulon') + '-track'
        x.update({'regulon': tmp})
        # blank out the "motifData" field for track-based regulons:
        x.update({'mofitData': 'NA.png'})

    # merge regulon threshold dictionaries:
    rt = rt_mtf + rt_trk

    ################################################################################
    # metadata
    ################################################################################

    metadata = {}

    metadata["metrics"] = [
        {
            "name": "nUMI"
        }, {
            "name": "nGene"
        }, {
            "name": "Percent_mito"
        }
    ]

    # SCENIC regulon thresholds:
    metadata["regulonThresholds"] = rt

    ################################################################################
    # row/column attributes
    ################################################################################

    # re-open the connection to the loom file to copy the original expression data
    lf = lp.connect(args.loom_motif, mode='r', validate=False)

    col_attrs = {
        "CellID": lf.ca.CellID,  # np.array(adata.obs.index),
        "RegulonsAUC": df_to_named_matrix(auc_mtx),
    }

    row_attrs = {
        "Gene": lf.ra.Gene,
        "Regulons": df_to_named_matrix(regulons),
    }

    attrs = {
        "MetaData": json.dumps(metadata),
    }

    if "SCopeTreeL1" in attrs.keys():
        attrs['SCopeTreeL1'] = lf.attrs.SCopeTreeL1
    if "SCopeTreeL2" in attrs.keys():
        attrs['SCopeTreeL2'] = lf.attrs.SCopeTreeL2
    if "SCopeTreeL3" in attrs.keys():
        attrs['SCopeTreeL3'] = lf.attrs.SCopeTreeL3

    lp.create(
        filename=args.loom_output,
        layers=lf[:, :],
        row_attrs=row_attrs,
        col_attrs=col_attrs,
        file_attrs=attrs
    )
    lf.close()


if __name__ == "__main__":
    integrate_motif_track(args)


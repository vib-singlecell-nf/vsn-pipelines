#!/usr/bin/env python3

import argparse
import base64
import json
import sys
import zlib
from copy import deepcopy
from shutil import copyfile

import loompy as lp
import numpy as np
import pandas as pd

################################################################################
################################################################################

parser = argparse.ArgumentParser(
    description='Integrate output from pySCENIC with SCope loom from pipeline best-practices path'
)
parser.add_argument(
    '--loom_scope',
    help='Loom file for SCope from pipeline best-practices path',
    required=True,
    default='.loom'
)
parser.add_argument(
    '--loom_scenic',
    help='Loom file from pySCENIC run',
    required=True,
    default='pyscenic_motif.loom'
)
parser.add_argument(
    '--loom_output',
    help='Final loom file with pySCENIC results integrated',
    required=True,
    default='pyscenic.loom'
)
args = parser.parse_args()


################################################################################
################################################################################

def dfToNamedMatrix(df):
    arr_ip = [tuple(i) for i in df.as_matrix()]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr


def integrateSCENICloom(args):
    ################################################################################
    # load data from scenic loom
    ################################################################################

    # scenic output
    lf = lp.connect(args.loom_scenic, mode='r', validate=False)
    meta_scenic = json.loads(zlib.decompress(base64.b64decode(lf.attrs.MetaData)))
    genes_scenic = lf.ra.Gene
    cellid_scenic = lf.ca.CellID
    regulonsAUC = lf.ca.RegulonsAUC
    regulons = lf.ra.Regulons

    ### capture embeddings:
    dr = [
        pd.DataFrame(lf.ca.Embedding, index=lf.ca.CellID)
    ]
    dr_names = [
        meta_scenic['embeddings'][0]['name'].replace(" ", "_")
    ]

    # add other embeddings
    drx = pd.DataFrame(lf.ca.Embeddings_X, index=lf.ca.CellID)
    dry = pd.DataFrame(lf.ca.Embeddings_Y, index=lf.ca.CellID)

    for i in range(len(drx.columns)):
        dr.append(pd.concat([drx.iloc[:, i], dry.iloc[:, i]], sort=False, axis=1, join='outer'))
        dr_names.append(meta_scenic['embeddings'][i + 1]['name'].replace(" ", "_").replace('/', '-'))

    # rename columns:
    for i, x in enumerate(dr):
        x.columns = ['X', 'Y']

    lf.close()

    ################################################################################
    # copy loom
    ################################################################################

    copyfile(args.loom_scope, args.loom_output)

    ################################################################################
    # add scenic data
    ################################################################################

    lf = lp.connect(args.loom_output, mode='r+', validate=False)

    if not all(lf.ca.CellID == cellid_scenic):
        sys.exit(f"ERROR: Column attribute CellIDs does not match between {args.loom_scope} and {args.loom_scenic}")
    if not all(lf.ra.Gene == genes_scenic):
        sys.exit(f"ERROR: Row attribute 'Gene' does not match between {args.loom_scope} and {args.loom_scenic}")

    # write regulon information:
    lf.ca['RegulonsAUC'] = regulonsAUC
    lf.ra['Regulons'] = regulons

    # get existing metadata:
    meta = json.loads(zlib.decompress(base64.b64decode(lf.attrs.MetaData)))

    # append regulon thresholds:
    meta['regulonThresholds'] = meta_scenic['regulonThresholds']

    embeddings_start_index = 1
    if len(meta['embeddings']) > 1:
        embeddings_start_index = int(meta['embeddings'][-1]['id'])

    # append embedding labels:
    id = deepcopy(embeddings_start_index)
    for i, x in enumerate(dr_names):
        id += 1
        meta['embeddings'].append({'id': id, 'name': x})

    # get existing embeddings:
    Embeddings_X = pd.DataFrame(lf.ca.Embeddings_X, index=lf.ca.CellID)
    Embeddings_Y = pd.DataFrame(lf.ca.Embeddings_Y, index=lf.ca.CellID)

    # append scenic embeddings:
    id = deepcopy(embeddings_start_index)
    for i, x in enumerate(dr):
        id += 1
        Embeddings_X[str(id)] = dr[i].iloc[:, 0]
        Embeddings_Y[str(id)] = dr[i].iloc[:, 1]

    lf.ca.Embeddings_X = dfToNamedMatrix(Embeddings_X)
    lf.ca.Embeddings_Y = dfToNamedMatrix(Embeddings_Y)

    lf.attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(meta).encode('ascii'))).decode('ascii')

    lf.close()


if __name__ == "__main__":
    integrateSCENICloom(args)

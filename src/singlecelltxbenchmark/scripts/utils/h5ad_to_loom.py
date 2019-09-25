#!/usr/bin/env python
import os
from optparse import OptionParser
import scanpy as sc
import loompy as lp
import pandas as pd
import zlib
import json
import base64
import numpy as np

parser = OptionParser(usage="usage: %prog [options] h5ad_file_path",
                      version="%prog 1.0")
parser.add_option("-x", "--method",
                  action="store",
                  dest="method",
                  default="linear_regression",
                  help="Normalize the data. Choose one of : regress_out")
parser.add_option("-r", "--variable-to-regress-out",
                  action="append",
                  dest="vars_to_regress_out",
                  default=None,
                  help="Variable to regress out. To regress multiple variables, add that many -v arguments. Used when running 'regress_out")
(options, args) = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args[0]
FILE_PATH_OUT_BASENAME = os.path.splitext(args[1])[0]


def dfToNamedMatrix(df):
    arr_ip = [tuple(i) for i in df.as_matrix()]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr

try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))


# ClusterMarkers_0 = pd.DataFrame(index=adata.raw.var.index, columns=[str(x) for x in range(max(set([int(x) for x in adata.obs['louvain']])) + 1)])

# ClusterMarkers_0_avg_logFC = pd.DataFrame(index=adata.raw.var.index, columns=[str(x) for x in range(max(set([int(x) for x in adata.obs['louvain']])) + 1)])

# ClusterMarkers_0_pval = pd.DataFrame(index=adata.raw.var.index, columns=[str(x) for x in range(max(set([int(x) for x in adata.obs['louvain']])) + 1)])

# ClusterMarkers_0.fillna(0, inplace=True)
# ClusterMarkers_0_avg_logFC.fillna(0, inplace=True)
# ClusterMarkers_0_pval.fillna(0, inplace=True)

# for i in range(max(set([int(x) for x in adata.obs['louvain']])) + 1):
#     i = str(i)
#     tot_genes = len(adata.uns['rank_genes_groups']['pvals_adj'][i])
#     sigGenes = adata.uns['rank_genes_groups']['pvals_adj'][i] < 0.05
#     deGenes = np.logical_and(np.logical_or(adata.uns['rank_genes_groups']['logfoldchanges'][i] >= 1.5, adata.uns['rank_genes_groups']['logfoldchanges'][i] <= -1.5), np.isfinite(adata.uns['rank_genes_groups']['logfoldchanges'][i]))
#     sigAndDE = np.logical_and(sigGenes, deGenes)
#     print(f'Filtering {sum(sigAndDE)} sig and de genes of {tot_genes}')

#     names = adata.uns['rank_genes_groups']['names'][i][sigAndDE]
#     ClusterMarkers_0.loc[names, i] = 1
#     ClusterMarkers_0_avg_logFC.loc[names, i] = np.around(adata.uns['rank_genes_groups']['logfoldchanges'][i][sigAndDE], decimals=6)
#     ClusterMarkers_0_pval.loc[names, i] = np.around(adata.uns['rank_genes_groups']['pvals_adj'][i][sigAndDE], decimals=6)

metaJson = {}
metaJson["metrics"] = []
metaJson["annotations"] = []

main_dr = pd.DataFrame(adata.obsm['X_umap'], columns=['_X', '_Y'])

metaJson['embeddings'] = [
    {
        "id": -1,
        "name": f"HVG UMAP"
    }]

Embeddings_X = pd.DataFrame()
Embeddings_Y = pd.DataFrame()

embeddings_id = 1

if 'X_tsne' in adata.obsm.keys():
    Embeddings_X[str(embeddings_id)] = pd.DataFrame(adata.obsm['X_tsne'])[0]
    Embeddings_Y[str(embeddings_id)] = pd.DataFrame(adata.obsm['X_tsne'])[1]
    embeddings_id += 1
    metaJson['embeddings'].append(
        {
            "id": embeddings_id,
            "name": f"HVG t-SNE"
        }
    )

Embeddings_X[str(embeddings_id)] = pd.DataFrame(adata.obsm['X_pca'])[0]
Embeddings_Y[str(embeddings_id)] = pd.DataFrame(adata.obsm['X_pca'])[1]
metaJson['embeddings'].append(
    {
        "id": embeddings_id,
        "name": f"HVG PC1/PC2"
    }
)

metaJson["clusterings"] = [{
            "id": 0,
            "group": "Louvain",
            "name": "Louvain default resolution",
            "clusters": [],                
        }]

for i in range(max(set([int(x) for x in adata.obs['louvain']])) + 1):
    clustDict = {}
    clustDict['id'] = i
    clustDict['description'] = f'Unannotated Cluster {i}'
    metaJson['clusterings'][0]['clusters'].append(clustDict)

clusterings = pd.DataFrame()

clusterings["0"] = adata.obs['louvain'].values.astype(np.int64)

col_attrs = {"CellID": np.array(adata.obs.index),
             "Embedding": dfToNamedMatrix(main_dr),
             "Embeddings_X": dfToNamedMatrix(Embeddings_X),
             "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
             "Clusterings": dfToNamedMatrix(clusterings),
             "ClusterID": np.array(adata.obs['louvain'].values)}


for col in adata.obs.keys():
    if type(adata.obs[col].dtype) == pd.core.dtypes.dtypes.CategoricalDtype:
        metaJson["annotations"].append(
            {
                "name": col,
                "values": list(set(adata.obs[col].values))
            }
        )
    else:
        metaJson["metrics"].append(
            {
                "name": col
            }
        )
    col_attrs[col] = np.array(adata.obs[col].values)


row_attrs = {"Gene": np.array(adata.raw.var.index)}

attrs = {"MetaData": json.dumps(metaJson)}

attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')

lp.create(filename=f"{FILE_PATH_OUT_BASENAME}.loom", layers=(adata.raw.X).T.toarray(), row_attrs=row_attrs, col_attrs=col_attrs, file_attrs=attrs)

#!/usr/bin/env python3

import argparse
import base64
import json
import loompy as lp
import numpy as np
import os
import pandas as pd
import scanpy as sc
import zlib

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "raw_filtered_data",
    type=argparse.FileType('r'),
    help='Input h5ad file containing the raw filtered data.'
)

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input h5ad file.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output .loom file.'
)

parser.add_argument(
    '--nomenclature',
    type=str,
    dest="nomenclature",
    help='The name of the genome.'
)

parser.add_argument(
    '--scope-tree-level-1',
    type=str,
    dest="scope_tree_level_1",
    help='The name of the first level of the SCope tree.'
)

parser.add_argument(
    '--scope-tree-level-2',
    type=str,
    dest="scope_tree_level_2",
    help='The name of the second level of the SCope tree.'
)

parser.add_argument(
    '--scope-tree-level-3',
    type=str,
    dest="scope_tree_level_3",
    help='The name of the third level of the SCope tree.'
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]


def df_to_named_matrix(df):
    arr_ip = [tuple(i) for i in df.as_matrix()]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr


try:
    raw_filtered_adata = sc.read_h5ad(filename=args.raw_filtered_data.name)
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

#####################
# Global Attributes #
#####################

# Initialize

attrs = {}
attrs_metadata = {}

attrs["Genome"] = '' if args.nomenclature is None else args.nomenclature
attrs["SCopeTreeL1"] = 'Unknown' if args.scope_tree_level_1 is None else args.scope_tree_level_1
attrs["SCopeTreeL2"] = '' if args.scope_tree_level_2 is None else args.scope_tree_level_2
attrs["SCopeTreeL3"] = '' if args.scope_tree_level_3 is None else args.scope_tree_level_3


#####################
# Column Attributes #
#####################

col_attrs = {
    "CellID": np.array(adata.obs.index)
}

# ANNOTATIONS & METRICS

attrs_metadata["metrics"] = []
attrs_metadata["annotations"] = []

# Populate
for column_attr_key in adata.obs.keys():
    if type(adata.obs[column_attr_key].dtype) == pd.core.dtypes.dtypes.CategoricalDtype:
        attrs_metadata["annotations"].append(
            {
                "name": column_attr_key,
                "values": list(set(adata.obs[column_attr_key].values))
            }
        )
    else:
        attrs_metadata["metrics"].append(
            {
                "name": column_attr_key
            }
        )
    col_attrs[column_attr_key] = np.array(adata.obs[column_attr_key].values)

# EMBEDDINGS

default_embedding = pd.DataFrame(
    index=adata.raw.obs_names
)
embeddings_x = pd.DataFrame(
    index=adata.raw.obs_names
)
embeddings_y = pd.DataFrame(
    index=adata.raw.obs_names
)

default_embedding["_X"] = adata.obsm['X_umap'][:, 0]
default_embedding["_Y"] = adata.obsm['X_umap'][:, 1]

embeddings_x["-1"] = adata.obsm['X_umap'][:, 0]
embeddings_y["-1"] = adata.obsm['X_umap'][:, 1]

attrs_metadata['embeddings'] = [
    {
        "id": -1,
        "name": f"HVG UMAP"
    }
]

if 'X_tsne' in adata.obsm.keys():
    embedding_id = embeddings_x.shape[1] - 1
    embeddings_x[str(embedding_id)] = adata.obsm['X_tsne'][:, 0]
    embeddings_y[str(embedding_id)] = adata.obsm['X_tsne'][:, 1]
    attrs_metadata['embeddings'].append(
        {
            "id": embedding_id,
            "name": f"HVG t-SNE"
        }
    )

if 'X_pca' in adata.obsm.keys():
    embedding_id = embeddings_x.shape[1] - 1
    embeddings_x[str(embedding_id)] = adata.obsm['X_pca'][:, 0]
    embeddings_y[str(embedding_id)] = adata.obsm['X_pca'][:, 1]
    attrs_metadata['embeddings'].append(
        {
            "id": embedding_id,
            "name": f"HVG PC1/PC2"
        }
    )

# Update column attribute Dict
col_attrs_embeddings = {
    "Embedding": df_to_named_matrix(default_embedding),
    "Embeddings_X": df_to_named_matrix(embeddings_x),
    "Embeddings_Y": df_to_named_matrix(embeddings_y)
}
col_attrs = {**col_attrs, **col_attrs_embeddings}

# CLUSTERINGS

clustering_id = 0

clustering_algorithm = adata.uns["rank_genes_groups"]["params"]["groupby"]
clustering_resolution = adata.uns[clustering_algorithm]["params"]["resolution"]
cluster_marker_method = adata.uns['rank_genes_groups']["params"]["method"]

num_clusters = int(max(adata.obs[clustering_algorithm])) + 1

clusterings = pd.DataFrame(
    index=adata.raw.obs_names
)
clusterings[str(clustering_id)] = adata.obs['leiden'].values.astype(np.int64)

attrs_metadata["clusterings"] = [
    {
        "id": clustering_id,
        "group": clustering_algorithm.capitalize(),
        "name": f"{clustering_algorithm.capitalize()} resolution {clustering_resolution}",
        "clusters": [],
        "clusterMarkerMetrics": [
            {
                "accessor": "avg_logFC",
                "name": "Avg. logFC",
                "description": f"Average log fold change from {cluster_marker_method.capitalize()} test"
            }, {
                "accessor": "pval",
                "name": "Adjusted P-Value",
                "description": f"Adjusted P-Value from {cluster_marker_method.capitalize()} test"
            }
        ]
    }
]

for i in range(num_clusters):
    cluster = {}
    cluster['id'] = i
    cluster['description'] = f'Unannotated Cluster {i}'
    attrs_metadata['clusterings'][clustering_id]['clusters'].append(cluster)


##################
# Row Attributes #
##################

row_attrs = {
    "Gene": np.array(adata.raw.var.index)
}

# CLUSTER MARKERS

# Initialize
cluster_markers = pd.DataFrame(
    index=adata.raw.var.index,
    columns=[str(x) for x in np.arange(num_clusters)]
).fillna(0, inplace=False)
cluster_markers_avg_logfc = pd.DataFrame(
    index=adata.raw.var.index,
    columns=[str(x) for x in np.arange(num_clusters)]
).fillna(0, inplace=False)
cluster_markers_pval = pd.DataFrame(
    index=adata.raw.var.index,
    columns=[str(x) for x in np.arange(num_clusters)]
).fillna(0, inplace=False)

# Populate
for i in range(num_clusters):
    i = str(i)
    num_genes = len(adata.uns['rank_genes_groups']['pvals_adj'][i])
    sig_genes_mask = adata.uns['rank_genes_groups']['pvals_adj'][i] < 0.05
    deg_genes_mask = np.logical_and(
        np.logical_or(
            adata.uns['rank_genes_groups']['logfoldchanges'][i] >= 1.5,
            adata.uns['rank_genes_groups']['logfoldchanges'][i] <= -1.5
        ),
        np.isfinite(
            adata.uns['rank_genes_groups']['logfoldchanges'][i]
        )
    )
    sig_and_deg_genes_mask = np.logical_and(
        sig_genes_mask,
        deg_genes_mask
    )
    gene_names = adata.uns['rank_genes_groups']['names'][i][sig_and_deg_genes_mask]
    cluster_markers.loc[gene_names, i] = 1
    cluster_markers_avg_logfc.loc[gene_names, i] = np.around(
        adata.uns['rank_genes_groups']['logfoldchanges'][i][sig_and_deg_genes_mask],
        decimals=6
    )
    cluster_markers_pval.loc[gene_names, i] = np.around(
        adata.uns['rank_genes_groups']['pvals_adj'][i][sig_and_deg_genes_mask],
        decimals=6
    )


# Update row attribute Dict
row_attrs_cluster_markers = {
    f"ClusterMarkers_{clustering_id}": df_to_named_matrix(cluster_markers),
    f"ClusterMarkers_{clustering_id}_avg_logFC": df_to_named_matrix(cluster_markers_avg_logfc),
    f"ClusterMarkers_{clustering_id}_pval": df_to_named_matrix(cluster_markers_pval)
}
row_attrs = {**row_attrs, **row_attrs_cluster_markers}

# Update global attribute Dict
attrs["MetaData"] = json.dumps(attrs_metadata)
# attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')

##################
# Build the Loom #
##################

lp.create(
    filename=f"{FILE_PATH_OUT_BASENAME}.loom",
    layers=(raw_filtered_adata.X).T,
    row_attrs=row_attrs,
    col_attrs=col_attrs,
    file_attrs=attrs
)

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
    "input",
    type=argparse.FileType('r'),
    nargs='+',
    help='Input h5ad file.'
)

parser.add_argument(
    "raw_filtered_data",
    type=argparse.FileType('r'),
    help='Input h5ad file containing the raw filtered data.'
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

parser.add_argument(
    '--markers-log-fc-threshold',
    type=float,
    default=0,
    dest="markers_log_fc_threshold",
    help='The name of the genome.'
)

parser.add_argument(
    '--markers-fdr-threshold',
    type=float,
    default=0.05,
    dest="markers_fdr_threshold",
    help='The name of the genome.'
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATHS_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]


def df_to_named_matrix(df):
    arr_ip = [tuple(i) for i in df.as_matrix()]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr


def read_h5ad(file_path, backed='r'):
    return sc.read_h5ad(filename=file_path, backed=backed)

# Loading
# - the raw filtered H5AD
# - the first of H5AD input files
# - all H5AD input files
#
# From the first H5AD input file one will get:
# - Annotations and metrics
# - Embeddings
# Assumption:
# - All H5AD input files differ only by their clustering and cluster markers


adatas = []

try:
    raw_filtered_adata = read_h5ad(file_path=args.raw_filtered_data.name, backed=None)
    adata = read_h5ad(file_path=FILE_PATHS_IN[0].name)
    for adata_idx in range(0, len(FILE_PATHS_IN)):
        adatas = adatas + [read_h5ad(file_path=FILE_PATHS_IN[adata_idx].name)]
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATHS_IN[0])[0]))

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
    # Don't store the clustering as annotation
    if type(adata.obs[column_attr_key].dtype) == pd.core.dtypes.dtypes.CategoricalDtype and column_attr_key != adata.uns["rank_genes_groups"]["params"]["groupby"]:
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

clusterings = pd.DataFrame(
    index=adata.raw.obs_names
)
attrs_metadata["clusterings"] = []

for adata_idx in range(0, len(FILE_PATHS_IN)):

    clustering_id = adata_idx
    clustering_algorithm = adatas[adata_idx].uns["rank_genes_groups"]["params"]["groupby"]
    clustering_resolution = adatas[adata_idx].uns[clustering_algorithm]["params"]["resolution"]
    cluster_marker_method = adatas[adata_idx].uns['rank_genes_groups']["params"]["method"]
    num_clusters = len(np.unique(adatas[adata_idx].obs[clustering_algorithm]))

    # Data
    clusterings[str(clustering_id)] = adatas[adata_idx].obs[clustering_algorithm].values.astype(np.int64)

    # Metadata
    attrs_metadata["clusterings"] = attrs_metadata["clusterings"] + [{
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
    }]

    for i in range(0, num_clusters):
        cluster = {}
        cluster['id'] = i
        cluster['description'] = f'Unannotated Cluster {i}'
        attrs_metadata['clusterings'][clustering_id]['clusters'].append(cluster)

# Update column attribute Dict
col_attrs_clusterings = {
    "ClusterID": clusterings["0"].values,  # Pick the first one as default clustering (this is purely arbitrary)
    "Clusterings": df_to_named_matrix(clusterings)
}
col_attrs = {**col_attrs, **col_attrs_clusterings}

##################
# Row Attributes #
##################

row_attrs = {
    "Gene": np.array(raw_filtered_adata.var.index)
}

# CLUSTER MARKERS

for adata_idx in range(0, len(FILE_PATHS_IN)):
    clustering_id = attrs_metadata['clusterings'][adata_idx]["id"]
    num_clusters = len(attrs_metadata['clusterings'][adata_idx]["clusters"])

    # Initialize
    cluster_markers = pd.DataFrame(
        index=raw_filtered_adata.var.index,
        columns=[str(x) for x in np.arange(num_clusters)]
    ).fillna(0, inplace=False)
    cluster_markers_avg_logfc = pd.DataFrame(
        index=raw_filtered_adata.var.index,
        columns=[str(x) for x in np.arange(num_clusters)]
    ).fillna(0, inplace=False)
    cluster_markers_pval = pd.DataFrame(
        index=raw_filtered_adata.var.index,
        columns=[str(x) for x in np.arange(num_clusters)]
    ).fillna(0, inplace=False)

    # Populate
    for i in range(0, num_clusters):
        i = str(i)
        gene_names = adatas[adata_idx].uns['rank_genes_groups']['names'][i]
        pvals_adj = adatas[adata_idx].uns['rank_genes_groups']['pvals_adj'][i]
        logfoldchanges = adatas[adata_idx].uns['rank_genes_groups']['logfoldchanges'][i]
        num_genes = len(gene_names)
        sig_genes_mask = pvals_adj < args.markers_fdr_threshold
        deg_genes_mask = np.logical_and(
            np.logical_or(
                logfoldchanges >= args.markers_log_fc_threshold,
                logfoldchanges <= -args.markers_log_fc_threshold
            ),
            np.isfinite(
                logfoldchanges
            )
        )
        sig_and_deg_genes_mask = np.logical_and(
            sig_genes_mask,
            deg_genes_mask
        )
        marker_names = gene_names[sig_and_deg_genes_mask]

        marker_genes_along_raw_adata_mask = np.in1d(
            raw_filtered_adata.var.index,
            marker_names
        )
        marker_genes_along_raw_adata = cluster_markers.index[marker_genes_along_raw_adata_mask]

        # Populate the marker mask
        markers_df = pd.DataFrame(
            1,
            index=marker_names,
            columns=["is_marker"]
        )
        cluster_markers.loc[
            marker_genes_along_raw_adata_mask,
            i
        ] = np.around(
            markers_df["is_marker"][marker_genes_along_raw_adata],
            decimals=6
        )

        # Populate the marker gene log fold changes
        logfoldchanges_df = pd.DataFrame(
            logfoldchanges[sig_and_deg_genes_mask],
            index=marker_names,
            columns=["logfc"]
        )
        cluster_markers_avg_logfc.loc[
            marker_genes_along_raw_adata_mask,
            i
        ] = logfoldchanges_df["logfc"][marker_genes_along_raw_adata]

        # Populate the marker gene false discovery rates
        pvals_adj_df = pd.DataFrame(
            pvals_adj[sig_and_deg_genes_mask],
            index=marker_names,
            columns=["fdr"]
        )
        cluster_markers_pval.loc[
            marker_genes_along_raw_adata_mask,
            i
        ] = pvals_adj_df["fdr"][marker_genes_along_raw_adata]

    # Update row attribute Dict
    row_attrs_cluster_markers = {
        f"ClusterMarkers_{str(clustering_id)}": df_to_named_matrix(cluster_markers),
        f"ClusterMarkers_{str(clustering_id)}_avg_logFC": df_to_named_matrix(cluster_markers_avg_logfc),
        f"ClusterMarkers_{str(clustering_id)}_pval": df_to_named_matrix(cluster_markers_pval)
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

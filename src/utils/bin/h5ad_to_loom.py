#!/usr/bin/env python3

import argparse
import base64
import json
import loompy as lp
import numpy as np
import os
import pandas as pd
from pandas.core.dtypes.dtypes import CategoricalDtypeType
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
    '--matrix-slot',
    type=str,
    dest="matrix_slot",
    choices=['X', 'raw_X'],
    default='X',
    help='The slot name where to find the raw matrix in raw_filtered_data. X refers to adata.X while X_raw refers to adata.raw.X.'
)

parser.add_argument(
    '--nomenclature',
    type=str,
    dest="nomenclature",
    help='The name of the genome annotation the genes come from. e.g.: Flybase r6.31.'
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
    help='Threshold on the log fold change for the markers to not to be saved in the loom.'
)

parser.add_argument(
    '--markers-fdr-threshold',
    type=float,
    default=0.05,
    dest="markers_fdr_threshold",
    help='Threshold on the false discovery rate (FDR) for the markers to not to be saved in the loom.'
)

args = parser.parse_args()

SCANPY__CLUSTER_MARKER_DATA__ANNDATA_UNS_KEY = "rank_genes_groups"

# Define the arguments properly
FILE_PATHS_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]


def df_to_named_matrix(df: pd.DataFrame):
    # Create meta-data structure.
    # Create a numpy structured array
    return np.array([tuple(row) for row in df.values],
                    dtype=np.dtype(list(zip(df.columns, df.dtypes))))


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

    # Check cell overlap between input and raw_filtered_adata
    overlap = np.in1d(
        adata.obs.index.values.astype(str),
        raw_filtered_adata.obs.index.values.astype(str)
    )
    num_overlap = np.sum(overlap)
    if num_overlap != adata.shape[0]:
        raise Exception("VSN ERROR: Some cells from input is not present in raw_filtered_data.")
    if num_overlap < raw_filtered_adata.shape[0]:
        print("VSN MSG: Subsetting the raw_filtered_data .h5ad file by the cells found in the input .h5ad file.")
        raw_filtered_adata = raw_filtered_adata[adata.obs.index, :]

    for adata_idx in range(0, len(FILE_PATHS_IN)):
        adatas = adatas + [read_h5ad(file_path=FILE_PATHS_IN[adata_idx].name)]
except IOError:
    raise Exception("VSN ERROR: Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATHS_IN[0])[0]))

CELL_IDS = adata.obs.sample_id.index

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


ANNOTATION_MAX_UNIQUE_VALUES = 1024


def is_metric(column_attr_key, values):
    uniq_values = np.unique(values)
    if column_attr_key.startswith('n_') or column_attr_key.startswith('num_') or column_attr_key.startswith('pct_') or column_attr_key.startswith('percent_'):
        return True
    if issubclass(values.dtype.type, np.integer):
        if len(uniq_values) > ANNOTATION_MAX_UNIQUE_VALUES:
            return True
    if issubclass(values.dtype.type, np.floating):
        return True
    return False


def is_annotation(column_attr_key, values):
    uniq_values = np.unique(values)
    if issubclass(values.dtype.type, CategoricalDtypeType):
        return True
    if issubclass(values.dtype.type, np.bool_):
        return True
    if issubclass(values.dtype.type, np.integer):
        if len(uniq_values) < ANNOTATION_MAX_UNIQUE_VALUES:
            return True
        return False
    if values.dtype == np.object and len(values.dtype) == 0 and len(uniq_values) < ANNOTATION_MAX_UNIQUE_VALUES:
        return True
    return False


# Populate
for column_attr_key in adata.obs.keys():
    # Don't store the clustering as annotation
    if column_attr_key.lower().startswith('cell_id') or column_attr_key.lower().startswith('cellid'):
        continue

    vals = adata.obs[column_attr_key].values
    uniq_vals = np.unique(vals)

    if is_metric(column_attr_key, vals):
        attrs_metadata["metrics"] = attrs_metadata["metrics"] + [
            {
                "name": column_attr_key
            }
        ]
    elif is_annotation(column_attr_key, vals):
        # Convert boolean to [0,1] range so that it matches with values in column attribute
        _uniq_vals = uniq_vals.astype(int) if issubclass(vals.dtype.type, np.bool_) else uniq_vals
        attrs_metadata["annotations"] = attrs_metadata["annotations"] + [
            {
                "name": column_attr_key,
                "values": list(map(lambda x: str(x), _uniq_vals.tolist()))
            }
        ]

    col_attrs[column_attr_key] = np.array(vals)

# EMBEDDINGS

default_embedding = pd.DataFrame(
    index=CELL_IDS
)
embeddings_x = pd.DataFrame(
    index=CELL_IDS
)
embeddings_y = pd.DataFrame(
    index=CELL_IDS
)

if 'X_umap' in adata.obsm.keys():
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
    index=CELL_IDS
)
attrs_metadata["clusterings"] = []

for adata_idx in range(0, len(FILE_PATHS_IN)):

    # Only add clustering information if rank_genes_group has been computed

    if SCANPY__CLUSTER_MARKER_DATA__ANNDATA_UNS_KEY not in adatas[adata_idx].uns:
        continue

    clustering_id = adata_idx
    clustering_algorithm = adatas[adata_idx].uns[SCANPY__CLUSTER_MARKER_DATA__ANNDATA_UNS_KEY]["params"]["groupby"]

    # Check if from AnnData 0.7.x
    if isinstance(clustering_algorithm, np.ndarray):
        if len(clustering_algorithm) > 1:
            raise Exception("VSN ERROR: Currently there is no support for conversion of h5ad containing multiple clusterings.")
        else:
            clustering_algorithm = clustering_algorithm[0]
    clustering_resolution = adatas[adata_idx].uns[clustering_algorithm]["params"]["resolution"]
    cluster_marker_method = adatas[adata_idx].uns[SCANPY__CLUSTER_MARKER_DATA__ANNDATA_UNS_KEY]["params"]["method"]

    # Check if from AnnData 0.7.x
    if isinstance(cluster_marker_method, np.ndarray):
        if len(cluster_marker_method) > 1:
            raise Exception("VSN ERROR: Currently there is no support for conversion of h5ad containing multiple differential expression results.")
        else:
            cluster_marker_method = cluster_marker_method[0]
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

    if SCANPY__CLUSTER_MARKER_DATA__ANNDATA_UNS_KEY not in adatas[adata_idx].uns:
        continue

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
        gene_names = adatas[adata_idx].uns[SCANPY__CLUSTER_MARKER_DATA__ANNDATA_UNS_KEY]['names'][i]
        pvals_adj = adatas[adata_idx].uns[SCANPY__CLUSTER_MARKER_DATA__ANNDATA_UNS_KEY]['pvals_adj'][i]
        logfoldchanges = adatas[adata_idx].uns[SCANPY__CLUSTER_MARKER_DATA__ANNDATA_UNS_KEY]['logfoldchanges'][i]
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
        ] = markers_df["is_marker"][marker_genes_along_raw_adata]

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

if args.matrix_slot == "X":
    print("Using X slot in AnnData...")
    raw_filtered_mat = (raw_filtered_adata.X).T
elif args.matrix_slot == "raw_X":
    print("Using raw.X slot in AnnData...")
    raw_filtered_mat = (raw_filtered_adata.raw.X).T
else:
    raise Exception("VSN ERROR: Invalid --matrix-slot. Choose X or X_raw.")

lp.create(
    filename=f"{FILE_PATH_OUT_BASENAME}.loom",
    layers=raw_filtered_mat,
    row_attrs=row_attrs,
    col_attrs=col_attrs,
    file_attrs=attrs
)

#!/usr/bin/env python3

import os
import re
import argparse
from shutil import copyfile
import json
import hdbscan
import pandas as pd
import loompy as lp
import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score

parser = argparse.ArgumentParser(description='Template script')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input Loom (following the SCope standards) file.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output Loom (following the SCope standards) file.'
)

parser.add_argument(
    "-e", "--cell-embeddings-index",
    type=int,
    dest="cell_embeddings_index",
    default=0,
    help='The index of the cell embeddings to use for the density-based spatial clustering.'
)

parser.add_argument(
    "-f", "--from-min-cluster-size",
    type=int,
    dest="from_min_cluster_size",
    default=5,
    help='The start value of the min_cluster_size to test'
)

parser.add_argument(
    "-t", "--to-min-cluster-size",
    type=int,
    dest="to_min_cluster_size",
    default=100,
    help='The end value of the min_cluster_size to test'
)

parser.add_argument(
    "-b", "--by-min-cluster-size",
    type=int,
    dest="by_min_cluster_size",
    default=5,
    help='The end value of the min_cluster_size to test'
)

parser.add_argument(
    "-x", "--from-min-samples",
    type=int,
    dest="from_min_samples",
    default=5,
    help='The start value of the min_samples to test'
)

parser.add_argument(
    "-y", "--to-min-samples",
    type=int,
    dest="to_min_samples",
    default=100,
    help='The end value of the min_samples to test'
)

parser.add_argument(
    "-z", "--by-min-samples",
    type=int,
    dest="by_min_samples",
    default=5,
    help='The end value of the min_samples to test'
)

args = parser.parse_args()


# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]
FILE_PATH_OUT = "{}.loom".format(FILE_PATH_OUT_BASENAME)
# I/O
# Expects loom file

# Copy loom which will use to update the metadata and set the default clustering
copyfile(FILE_PATH_IN.name, FILE_PATH_OUT)

try:
    loom = lp.connect(FILE_PATH_OUT, validate=False)
except IOError:
    raise Exception("Wrong input format. Expects .loom files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

################################################################################
# Process the data...

CELL_EMBEDDINGS_INDEX = str(args.cell_embeddings_index)

if CELL_EMBEDDINGS_INDEX not in loom.ca.Embeddings_X.dtype.names or CELL_EMBEDDINGS_INDEX not in loom.ca.Embeddings_Y.dtype.names:
    raise Exception(f"The given cell embeddings index {args.cell_embeddings_index} does not exist in the embeddings of the given loom file {FILE_PATH_IN}.")

# Get the cell embeddings
md = json.loads(loom.attrs.MetaData)
cell_embeddings = pd.DataFrame({
    "_X": loom.ca.Embeddings_X[CELL_EMBEDDINGS_INDEX],
    "_Y": loom.ca.Embeddings_Y[CELL_EMBEDDINGS_INDEX]
}, index=loom.ca.CellID)

# Get the clusterings array
clusterings = loom.ca.Clusterings

grid_res = []
for min_cluster_size in range(args.from_min_cluster_size, args.to_min_cluster_size + args.by_min_cluster_size, args.by_min_cluster_size):
    for min_samples in range(args.from_min_samples, args.to_min_samples + args.by_min_samples, args.by_min_samples):
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            gen_min_span_tree=True
        )
        clusterer.fit(cell_embeddings)
        res = []
        for i in loom.ca.Clusterings.dtype.names:
            clustering_name = list(filter(lambda x: x["id"] == int(i), md["clusterings"]))[0]['name']
            resolution = re.sub(r'.* resolution ', '', clustering_name)
            ars = adjusted_rand_score(clusterer.labels_, clusterings[i])
            res = res + [(clustering_name, i, resolution, ars)]
        opt_res = sorted(res, key=lambda tup: tup[3], reverse=True)[0]
        grid_res = grid_res + [(min_cluster_size, min_samples) + opt_res]
        print(f"Optimal clustering resolution for HDBSCAN with min_cluster_size={min_cluster_size}, min_samples={min_samples} => {opt_res}")

# Aggregating the results
grid_res_df = pd.DataFrame(
    grid_res,
    columns=[
        "min_cluster_size",
        "min_samples",
        "clustering_name",
        "clustering_id",
        "opt_clustering_resolution",
        "opt_clustering_adjusted_rand_score"
    ]
)

# Get optimal clustering
vc = grid_res_df['clustering_id'].value_counts()
print(f"")
print(f"============================")
print(f"Clustering Occurrence Table.")
print(f"============================")
print(f"")
print(f"Clustering Name => Occurrence")
for idx in range(0, len(vc.index)):
    clustering_id = vc.index[idx]
    clustering_name = md["clusterings"][int(clustering_id)]["name"]
    occurrence = vc.values[idx]
    print(f"{clustering_name} => {occurrence}")
opt_clustering_id = vc.head(1).index[0]
opt_clustering_name = md["clusterings"][int(opt_clustering_id)]["name"]
print(f"Optimal clustering: {opt_clustering_name}")


md_clusterings_sorted_by_resolution = sorted(
    md["clusterings"],
    key=lambda k: float(re.sub(r'.* resolution ', '', k['name']))
)


# Update the clustering metadata
def update_md_clustering(md_clustering, idx, occurrence):
    md_clustering_copy = md_clustering.copy()
    if idx == 0:
        md_clustering_copy["name"] = md_clustering_copy["name"] + f" (default, {occurrence})"
    else:
        md_clustering_copy["name"] = md_clustering_copy["name"] + f" ({occurrence})"
    return md_clustering_copy


# Update the clusterings metadata
def update_md_clusterings(md_clusterings, idx, clustering_id, occurrence):
    md_clusterings_copy = md_clusterings.copy()
    new_md_clusterings = list(map(
        lambda x: update_md_clustering(md_clustering=x, idx=idx, occurrence=occurrence) if x["id"] == int(clustering_id) else x,
        md_clusterings_copy
    ))
    return new_md_clusterings


new_md_clusterings = md_clusterings_sorted_by_resolution.copy()

for idx in range(0, len(vc.index)):
    new_md_clusterings = update_md_clusterings(
        md_clusterings=new_md_clusterings,
        idx=idx,
        clustering_id=vc.index[idx],
        occurrence=vc.values[idx]
    )

md["clusterings"] = new_md_clusterings
loom.attrs.MetaData = json.dumps(md)

################################################################################

# I/O
loom.close()

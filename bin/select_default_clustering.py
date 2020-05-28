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

# Get the cell embeddings
md = json.loads(loom.attrs.MetaData)
cell_embeddings = pd.DataFrame({
    "_X": loom.ca.Embeddings_X['0'],
    "_Y": loom.ca.Embeddings_Y['0']
}, index=loom.ca.CellID)

# Get the clusterings array
clusterings = loom.ca.Clusterings

grid_res = []
for min_cluster_size in range(args.from_min_cluster_size, args.to_min_cluster_size + args.by_min_cluster_size, args.by_min_cluster_size):
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        gen_min_span_tree=True
    )
    clusterer.fit(cell_embeddings)
    res = []
    for i in loom.ca.Clusterings.dtype.names:
        clustering_name = list(filter(lambda x: x["id"] == int(i), md["clusterings"]))[0]['name']
        resolution = re.sub(r'.* resolution ', '', clustering_name)
        ars = adjusted_rand_score(clusterer.labels_, clusterings[i])
        res = res + [(i, resolution, ars)]
    opt_res = sorted(res, key=lambda tup: tup[2], reverse=True)[0]
    grid_res = grid_res + [(min_cluster_size, ) + opt_res]
    print(f"Optimal clustering resolution for HDBSCAN with min_cluster_size={min_cluster_size} => {opt_res}")

# Aggregating the results
grid_res_df = pd.DataFrame(
    grid_res,
    columns=[
        "min_cluster_size",
        "clustering_id",
        "opt_clustering_resolution",
        "opt_clustering_adjusted_rand_score"
    ]
)

# Get optimal clustering
vc = grid_res_df['clustering_id'].value_counts()
print(vc)
opt_clustering_id = vc.head(1).index[0]


# Update the clusterings metadata
def set_default_md_clustering(md_clustering):
    mc_clustering_copy = md_clustering.copy()
    mc_clustering_copy["name"] = mc_clustering_copy["name"] + " (default)"
    return mc_clustering_copy


md_clusterings_sorted_by_resolution = sorted(
    md["clusterings"],
    key=lambda k: float(re.sub(r'.* resolution ', '', k['name']))
)
new_md_clusterings = list(map(
    lambda x: set_default_md_clustering(md_clustering=x) if x["id"] == int(opt_clustering_id) else x,
    md_clusterings_sorted_by_resolution
))
md["clusterings"] = new_md_clusterings
loom.attrs.MetaData = json.dumps(md)

################################################################################

# I/O
loom.close()

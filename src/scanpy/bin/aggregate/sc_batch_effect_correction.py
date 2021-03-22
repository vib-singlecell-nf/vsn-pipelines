#!/usr/bin/env python3

import argparse
import numpy as np
import os
import scanpy as sc

print('''Correct the data from batch effects. Choose one of:
- bbknn:
    Batch balanced kNN alters the kNN procedure to identify each cell's top neighbours in each batch separately instead of the entire cell pool with no accounting for batch.
    Details at: https://icb-scanpy.readthedocs-hosted.com/en/stable/external/scanpy.external.pp.bbknn.html
- mnncorrect:
    Correct batch effects by matching mutual nearest neighbors.
    Details at: see https://icb-scanpy.readthedocs-hosted.com/en/stable/external/scanpy.external.pp.mnncorrect.html
''')

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    nargs='+',
    type=argparse.FileType('r'),
    help='Input h5ad files.'
)

parser.add_argument(
    "-x", "--method",
    type=str,
    action="store",
    dest="method",
    default="bbknn",
    help="Correct the data from batch effects. Choose one of: bkknn, mnncorrect."
)

parser.add_argument(
    "-o", "--output-file",
    type=str,
    action="store",
    dest="output_file",
    default='',
    help="Output file name."
)

parser.add_argument(
    "-b", "--batch-key",
    type=str,
    action="store",
    dest="batch_key",
    default='batch',
    help="[bbknn, mnncorrect], adata.obs column name discriminating between your batches."
)

parser.add_argument(
    "-K", "--key",
    type=str,
    action="store",
    dest="key",
    default='batch',
    help="[combat], adata.obs column name discriminating between your batches."
)

parser.add_argument(
    "-c", "--n-pcs",
    type=int,
    action="store",
    dest="n_pcs",
    default=50,
    help="[bbknn]. How many principal components to use in the analysis."
)

parser.add_argument(
    "-k", "--k",
    type=int,
    action="store",
    dest="k",
    default=20,
    help="[mnncorrect]. Number of mutual nearest neighbors."
)

parser.add_argument(
    "-i", "--var-index",
    type=str,
    action="store",
    dest="var_index",
    default=None,
    help="[mnncorrect]. The index (list of str) of vars (genes). Necessary when using only a subset of vars to perform MNN correction, and should be supplied with var_subset."
)

parser.add_argument(
    "-s", "--var-subset",
    action="store",
    dest="var_subset",
    default=None,
    help="[mnncorrect]. The subset of vars (list of str) to be used when performing MNN correction. Typically, a list of highly variable genes (HVGs). When set to None, uses all vars."
)

parser.add_argument(
    "-j", "--n-jobs",
    type=int,
    action="store",
    dest="n_jobs",
    default=None,
    help="[bbknn, mnncorrect], The number of jobs. When set to None, automatically uses the number of cores."
)

parser.add_argument(
    "-t", "--trim",
    type=int,
    action="store",
    dest="trim",
    default=None,
    help="[bbknn], Trim the neighbours of each cell to these many top connectivities. May help with population independence and improve the tidiness of clustering."
)

parser.add_argument(
    "-n", "--neighbors-within-batch",
    type=int,
    action="store",
    dest="neighbors_within_batch",
    default=3,
    help="[bbknn], How many top neighbours to report for each batch; total number of neighbours will be this number times the number of batches."
)

args = parser.parse_args()

# Define the arguments properly
DATA_FILE_PATH_BASENAME = os.path.splitext(args.output_file)[0]

# I/O
# Expects h5ad file
# Get all input files into a list
adatas = []
for FILE_PATH_IN in args.input:
    try:
        FILE_PATH_IN = FILE_PATH_IN.name
        adata = sc.read_h5ad(filename=FILE_PATH_IN)
        adatas.append(adata)
    except IOError:
        raise Exception("VSN ERROR: Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

#
# Aggregate the data
#

if args.method == 'combat':
    sc.pp.combat(
        adatas[0],
        key=args.key
    )
elif args.method == 'bbknn':
    # Expects:
    # - the PCA to have been computed and stored in adata.obsm['X_pca']
    if 'X_pca' not in adata.obsm.keys():
        raise Exception("VSN ERROR: Expect the PCA to have been computed and stored in adata.obsm['X_pca']")
    # Run BBKNN
    sc.external.pp.bbknn(
        adatas[0],
        batch_key=args.batch_key,
        n_pcs=args.n_pcs,
        neighbors_within_batch=args.neighbors_within_batch,
        trim=args.trim
    )
elif args.method == 'mnncorrect':
    # Run MNN_CORRECT (mnnpy)
    # GitHub: https://github.com/chriscainx/mnnpy/tree/master
    # Open issue: https://github.com/chriscainx/mnnpy/issues/24
    #  Use numba==0.43.1
    from mnnpy import mnn_correct

    _adatas = []
    cell_ids = []
    print(f"Splitting the AnnData into batches defined by {args.batch_key}...")
    for batch in np.unique(adata.obs[args.batch_key]):
        _adata_by_batch = adata[adata.obs[args.batch_key] == batch, :]
        cell_ids.extend(_adata_by_batch.obs.index.values)
        _adatas.append(
            _adata_by_batch
        )

    # Get all HVG into a list
    # Union
    print(f"Getting all highly variable genes...")
    hvgs = []
    for ADATA in _adatas:
        hvgs = np.union1d(hvgs, ADATA.var.index[ADATA.var.highly_variable])

    index_unique = None

    if len(cell_ids) != len(np.unique(cell_ids)):
        print("Non-unique cell index detected!")
        print("Make the index unique by joining the existing index names with the batch category, using index_unique='-'")
        index_unique = '-'

    # [mnn_correct] Subset only HVG otherwise "ValueError: Lengths must match to compare"
    print(f"Perform mnnCorrect...")
    corrected = mnn_correct(
        *_adatas,
        var_index=args.var_index,
        var_subset=hvgs if args.var_subset is None else args.var_subset,
        batch_key=args.batch_key,
        index_unique=index_unique,
        k=args.k,
        n_jobs=args.n_jobs)
    adata = corrected[0]
    # Run MNN_CORRECT (mnnpy)
    # Open GitHub issue: https://github.com/theislab/scanpy/issues/757
    # sc.external.pp.mnn_correct(adata,
    #     var_index=options.var_index,
    #     var_subset=options.var_subset,
    #     batch_key=options.batch_key,
    #     k=options.k)
else:
    raise Exception(f"The given batch effect correction method {args.method} is not implemented.")

# I/O
adata.write_h5ad("{}.h5ad".format(DATA_FILE_PATH_BASENAME))

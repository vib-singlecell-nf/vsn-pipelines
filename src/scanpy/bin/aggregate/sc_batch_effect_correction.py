#!/usr/bin/env python
import os
from optparse import OptionParser

import numpy as np

import scanpy as sc

print('''Correct the data from batch effects. Choose one of:
- bbknn:
    Batch balanced kNN alters the kNN procedure to identify each cell's top neighbours in each batch separately instead of the entire cell pool with no accounting for batch.
    Details at: https://icb-scanpy.readthedocs-hosted.com/en/stable/external/scanpy.external.pp.bbknn.html
- mnn_correct:
    Correct batch effects by matching mutual nearest neighbors.
    Details at: see https://icb-scanpy.readthedocs-hosted.com/en/stable/external/scanpy.external.pp.mnn_correct.html
''')

parser = OptionParser(
    usage="usage: %prog [options] h5ad_file_path",
    version="%prog 1.0"
)
parser.add_option(
    "-x", "--method",
    type="string",
    action="store",
    dest="method",
    default="bbknn",
    help="Correct the data from batch effects. Choose one of: bkknn, mnn_correct."
)
parser.add_option(
    "-o", "--output-file",
    type="string",
    action="store",
    dest="output_file",
    default='',
    help="Output file name."
)
parser.add_option(
    "-b", "--batch-key",
    type="string",
    action="store",
    dest="batch_key",
    default='batch',
    help="[bbknn, mnn_correct], adata.obs column name discriminating between your batches."
)
parser.add_option(
    "-K", "--key",
    type="string",
    action="store",
    dest="key",
    default='batch',
    help="[combat], adata.obs column name discriminating between your batches."
)
parser.add_option(
    "-c", "--n-pcs",
    type="int",
    action="store",
    dest="n_pcs",
    default=50,
    help="[bbknn]. How many principal components to use in the analysis."
)
parser.add_option(
    "-k", "--k",
    type="int",
    action="store",
    dest="k",
    default=20,
    help="[mnn_correct]. Number of mutual nearest neighbors."
)
parser.add_option(
    "-i", "--var-index",
    type="string",
    action="store",
    dest="var_index",
    default=None,
    help="[mnn_correct]. The index (list of str) of vars (genes). Necessary when using only a subset of vars to perform MNN correction, and should be supplied with var_subset."
)
parser.add_option(
    "-s", "--var-subset",
    action="store",
    dest="var_subset",
    default=None,
    help="[mnn_correct]. The subset of vars (list of str) to be used when performing MNN correction. Typically, a list of highly variable genes (HVGs). When set to None, uses all vars."
)
parser.add_option(
    "-j", "--n-jobs",
    type="int",
    action="store",
    dest="n_jobs",
    default=None,
    help="[bbknn, mnn_correct], The number of jobs. When set to None, automatically uses the number of cores."
)
parser.add_option(
    "-t", "--trim",
    type="int",
    action="store",
    dest="trim",
    default=None,
    help="[bbknn], Trim the neighbours of each cell to these many top connectivities. May help with population independence and improve the tidiness of clustering."
)
parser.add_option(
    "-n", "--neighbors-within-batch",
    type="int",
    action="store",
    dest="neighbors_within_batch",
    default=3,
    help="[bbknn], How many top neighbours to report for each batch; total number of neighbours will be this number times the number of batches."
)

(options, args) = parser.parse_args()

# Define the arguments properly
DATA_FILE_PATH_BASENAME = os.path.splitext(options.output_file)[0]

# I/O
# Expects h5ad file
# Get all input files into a list
adatas = []
for FILE_PATH_IN in args:
    try:
        adata = sc.read_h5ad(filename=FILE_PATH_IN)
        adatas.append(adata)
    except Exception:
        raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

# Get all HVG into a list
# Union
hvgs = []
for ADATA in adatas:
    hvgs = np.union1d(hvgs, ADATA.var.index[ADATA.var.highly_variable])

#
# Aggregate the data
#

if options.method == 'combat':
    sc.pp.combat(
        adatas[0],
        key=options.key)
elif options.method == 'bbknn':
    # Expects:
    # - the PCA to have been computed and stored in adata.obsm['X_pca']
    if 'X_pca' not in adata.obsm.keys():
        raise Exception("Expect the PCA to have been computed and stored in adata.obsm['X_pca']")
    # Run BBKNN
    sc.external.pp.bbknn(
        adatas[0],
        batch_key=options.batch_key,
        n_pcs=options.n_pcs,
        neighbors_within_batch=options.neighbors_within_batch,
        trim=options.trim)
elif options.method == 'mnn':
    # Run MNN_CORRECT (mnnpy)
    # GitHub: https://github.com/chriscainx/mnnpy/tree/master
    # Open issue: https://github.com/chriscainx/mnnpy/issues/24
    #  Use numba==0.43.1
    from mnnpy import mnn_correct

    # [mnn_correct] Subset only HVG otherwise "ValueError: Lengths must match to compare"
    corrected = mnn_correct(
        *list(map(lambda adata: adata[:, adata.var.index.isin(hvgs)], adatas)),
        var_index=options.var_index,
        var_subset=hvgs if options.var_subset is None else options.var_subset,
        batch_key=options.batch_key,
        k=options.k,
        n_jobs=options.n_jobs)
    adata = corrected[0]
    # Run MNN_CORRECT (mnnpy)
    # Open GitHub issue: https://github.com/theislab/scanpy/issues/757
    # sc.external.pp.mnn_correct(adata,
    #     var_index=options.var_index,
    #     var_subset=options.var_subset,
    #     batch_key=options.batch_key,
    #     k=options.k)
else:
    raise Exception("Method does not exist.")

# I/O
adata.write_h5ad("{}.h5ad".format(DATA_FILE_PATH_BASENAME))

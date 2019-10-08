#!/usr/bin/env python
import os
from optparse import OptionParser
import scanpy as sc
import anndata as ad
import warnings

parser = OptionParser(usage="usage: %prog [options] h5ad_file_path",
                      version="%prog 1.0")
parser.add_option("-x", "--method",
                  type="string",
                  action="store",
                  dest="method",
                  default="wilcoxon",
                  help="The default 'wilcoxon' uses Wilcoxon rank-sum, ‘t-test_overestim_var’ overestimates variance of each group, 't-test' uses t-test,  'logreg' uses logistic regression. See [Ntranos18], here and here, for why this is meaningful.")
parser.add_option("-g", "--groupby",
                  type="string",
                  action="store",
                  dest="groupby",
                  default="louvain",
                  help="The key of the observations grouping to consider.")
parser.add_option("-n", "--ngenes",
                  type="int",
                  action="store",
                  dest="ngenes",
                  default=0,
                  help="The number of genes that appear in the returned tables. Value of 0 will report all genes.")

(options, args) = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args[0]
FILE_PATH_OUT_BASENAME = os.path.splitext(args[1])[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN)
except:
    raise Exception("Can only handle .h5ad files.")

if 'raw' not in dir(adata):
    warnings.warn("There is no raw attribute in your anndata object. Differential analysis will be performed on the main (possibly normalised) matrix.")

if options.ngenes == 0:
    try:
        ngenes = adata.raw.shape[1]
    except AttributeError:
        ngenes = adata.shape[1]
        warnings.warn("There is no raw attribute in your anndata object. Using the shape of the main matrix as the number of genes to report.")


sc.tl.rank_genes_groups(adata, groupby=options.groupby, method=options.method, n_genes=ngenes)

adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

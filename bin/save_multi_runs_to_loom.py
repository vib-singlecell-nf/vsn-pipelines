#!/usr/bin/env python3

import argparse
import re
import sys
import pandas as pd
import pickle
import gzip
import time
import utils
import export_to_loom
import warnings
import os

parser = argparse.ArgumentParser(description='Run AUCell on gene signatures saved as .tsv in folder.')

parser.add_argument(
    'expression_mtx_fname',
    type=argparse.FileType('r'),
    help='The name of the file that contains the expression matrix for the single cell experiment.'
         ' Two file formats are supported: csv (rows=cells x columns=genes) or loom (rows=genes x columns=cells).'
)

parser.add_argument(
    'aggregated_regulons_fname',
    type=argparse.FileType('r'),
    help='The name of the file (.pkl.gz aka compressed pickle format) that contains the regulons resulting from the aggregated motif enrichment table.'
)

parser.add_argument(
    'auc_mtx_fname',
    type=argparse.FileType('r'),
    help='The name of the file that contains the AUCell matrix.'
)

parser.add_argument(
    '-o', '--output',
    type=argparse.FileType('w'),
    default=sys.stdout,
    help='Output file/stream, i.e. a table of TF-target genes (CSV).'
)

parser.add_argument(
    '--min-genes-regulon',
    type=int,
    default=5,
    dest="min_genes_regulon",
    help='The threshold used for filtering the regulons based on the number of targets (default: {}).'.format(5)
)

parser.add_argument(
    '--min-regulon-gene-occurrence',
    type=int,
    default=5,
    dest="min_regulon_gene_occurrence",
    help='The threshold used for filtering the genes bases on their occurrence (default: {}).'.format(5)
)

parser.add_argument(
    '--cell-id-attribute',
    type=str,
    default='CellID',
    dest="cell_id_attribute",
    help='The name of the column attribute that specifies the identifiers of the cells in the loom file.'
)

parser.add_argument(
    '--gene-attribute',
    type=str,
    default='Gene',
    dest="gene_attribute",
    help='The name of the row attribute that specifies the gene symbols in the loom file.'
)

parser.add_argument(
    '--title',
    type=str,
    dest="title",
    help='The title for this loom file. If None than the basename of the filename is used as the title.'
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

print(f"Extracting the matrix form the loom...", flush=True)
start = time.time()
ex_matrix_df = utils.get_matrix(
    loom_file_path=args.expression_mtx_fname.name,
    gene_attribute=args.gene_attribute,
    cell_id_attribute=args.cell_id_attribute
)
print(f"... took {time.time() - start} seconds", flush=True)

print(f"Reading the aggregated regulons...", flush=True)
start = time.time()
with gzip.open(args.aggregated_regulons_fname.name, 'rb') as file_handler:
    regulons = pickle.load(file_handler)
print(f"... took {time.time() - start} seconds to run.", flush=True)

print(f"Reading AUCell matrix...", flush=True)
start = time.time()
# Read the regulons AUCell matrix
auc_mtx = pd.read_csv(args.auc_mtx_fname.name, sep='\t', header=0, index_col=0)
auc_mtx.columns.name = "Regulon"
print(f"... took {time.time() - start} seconds to run.", flush=True)

# Check whether the cell index between AUC matrix and expression matrix are the same
if len(auc_mtx.index.difference(ex_matrix_df.index)) > 0:
    warnings.warn("Difference detected in the cell indexes between the AUCell matrix and the expression matrix. Subsetting cells from AUC matrix...", Warning)
    auc_mtx = auc_mtx[auc_mtx.index.isin(ex_matrix_df.index)]

# Filter regulons based on the ones used in AUCell matrix
# In case of a lot of regulons (>800), this is required to avoid header limitation set by h5py
regulons = list(filter(lambda x: x.name in auc_mtx.columns, regulons))

# Create loom
print(f"Exporting to loom...", flush=True)
start = time.time()
# Create the basic loom
scope_loom = export_to_loom.SCopeLoom(
    ex_mtx=ex_matrix_df,
    regulons=regulons,
    title=args.title,
    nomenclature=args.nomenclature,
    auc_mtx=auc_mtx,
    tree_structure=[
        args.scope_tree_level_1,
        args.scope_tree_level_2,
        args.scope_tree_level_3
    ],
    compress=True,
    save_additional_regulon_meta_data=True
)
scope_loom.set_generic_loom()
# Add additional stuff specific to multi-runs SCENIC
scope_loom.add_row_attr_regulon_gene_weights()
scope_loom.add_row_attr_regulon_gene_occurrences()
scope_loom.set_scenic_min_genes_regulon(min_genes_regulon=args.min_genes_regulon)
scope_loom.set_scenic_min_regulon_gene_occurrence(min_regulon_gene_occurrence=args.min_regulon_gene_occurrence)
scope_loom.export(out_fname=args.output.name, save_embeddings=False)

print(f"... took {time.time() - start} seconds to run.", flush=True)
print(f"Done.", flush=True)

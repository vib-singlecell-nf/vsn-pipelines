#!/usr/bin/env python3

import argparse
import glob
import loompy as lp
import ntpath
import os
import pandas as pd
from pyscenic.genesig import GeneSignature
from pyscenic.aucell import aucell, derive_auc_threshold, enrichment, create_rankings
import sys
from typing import List
import warnings

################################################################################
################################################################################

parser_grn = argparse.ArgumentParser(description='Run AUCell on gene signatures saved as TSV in folder.')

parser_grn.add_argument('expression_mtx_fname',
                        type=argparse.FileType('r'),
                        help='The name of the file that contains the expression matrix for the single cell experiment.'
                        ' Two file formats are supported: csv (rows=cells x columns=genes) or loom (rows=genes x columns=cells).')
parser_grn.add_argument('signatures_folder',
                        help='The folder containing the signatures as TSV files.')
parser_grn.add_argument('-o', '--output',
                        type=argparse.FileType('w'), default=sys.stdout,
                        help='Output file/stream, i.e. a table of TF-target genes (CSV).')
parser_grn.add_argument('--regulon-type',
                        type=str,
                        dest="regulon_type",
                        help='The type of the database ( mtf or trk) used for cisTarget pruning step.')
parser_grn.add_argument('--auc-threshold',
                        type=float, default=0.05,
                        dest="auc_threshold",
                        help='The auc threshold used for calculating the AUC of a feature as fraction of ranked genes (default: {}).'.format(0.05))
parser_grn.add_argument('--percentile-threshold',
                        type=float,
                        dest="percentile_threshold",
                        help='The percentile to get the auc threshold used for calculating the AUC of a feature as fraction of ranked genes (default: {}).'.format(0.05))
parser_grn.add_argument('--gene-occurence-threshold',
                        type=int, default=5,
                        dest="gene_occurence_threshold",
                        help='The threshold used for filtering the genes bases on their occurence (default: {}).'.format(5))
parser_grn.add_argument('--num-workers',
                        type=int, default=1,
                        dest="num_workers",
                        help='The number of workers to use. (default: {}).'.format(1))
parser_grn.add_argument('--cell-id-attribute',
                        type=str, default='CellID',
                        dest="cell_id_attribute",
                        help='The name of the column attribute that specifies the identifiers of the cells in the loom file.')
parser_grn.add_argument('--gene-attribute',
                        type=str, default='Gene',
                        dest="gene_attribute",
                        help='The name of the row attribute that specifies the gene symbols in the loom file.')

args = parser_grn.parse_args()

# Do stuff


def read_signatures_from_tsv_dir(dpath: str, noweights=False, weight_threshold=0, show_warnings=False) -> List['GeneSignature']:
    """
    Load gene signatures from a list of TSV files in directory. Requires TSV with 1 or 2 columns. First column should be genes, second (optional) are weight for genes.
    :param dpath: The filepath to directory.
    :return: A list of signatures.
    """
    # https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
    assert os.path.exists(dpath), "{} does not exist.".format(dpath)
    gene_sig_file_paths = glob.glob(os.path.join(dpath, "*.tsv"))

    def signatures():
        for gene_sig_file_path in gene_sig_file_paths:
            fname = ntpath.basename(gene_sig_file_path)
            regulon = os.path.splitext(fname)[0]
            gene_sig = pd.read_csv(gene_sig_file_path, sep='\t', header=None, index_col=None)
            # Do some sanity checks
            if len(gene_sig.columns) == 0:
                assert os.path.exists(gene_sig_file_path), "{} has 0 columns. Requires TSV with 1 or 2 columns. First column should be genes (required), second (optional) are weight for the given genes.".format(gene_sig_file_path)
            if len(gene_sig.columns) > 2:
                assert os.path.exists(gene_sig_file_path), "{} has more than 2 columns. Requires TSV with 1 or 2 columns. First column should be genes, second (optional) are weight for the given genes.".format(gene_sig_file_path)
            if len(gene_sig.columns) == 1 or noweights:
                gene2weight = gene_sig[0]
            if len(gene_sig.columns) == 2 and not noweights:
                # Filter the genes based on the given weight_threshold
                # print(weight_threshold)
                gene_sig = gene_sig[gene_sig[1] > weight_threshold]
                if len(gene_sig.index) == 0:
                    if show_warnings:
                        warnings.warn("{0} is empty after apply filter with weight_threshold > {1}".format(regulon, weight_threshold))
                    continue
                gene2weight = [tuple(x) for x in gene_sig.values]
            yield GeneSignature(name=regulon, gene2weight=gene2weight)
    signatures = list(signatures())
    print("Signatures passed filtering {0} out of {1}".format(len(signatures), len(gene_sig_file_paths)))
    return signatures


def get_matrix(loom_file_path):
    with lp.connect(loom_file_path, mode='r', validate=False) as loom:  # Read in r mode otherwise concurrency problem
        ex_matrix = loom[:, :]
        ex_matrix_df = pd.DataFrame(data=ex_matrix[:, :], index=loom.ra[args.gene_attribute], columns=loom.ca[args.cell_id_attribute]).T  # Gene expression as (cell, gene) - matrix.
    return ex_matrix_df

ex_matrix_df = get_matrix(loom_file_path=args.expression_mtx_fname.name)
signatures = read_signatures_from_tsv_dir(dpath=args.signatures_folder, noweights=False, weight_threshold=args.gene_occurence_threshold)
auc_threshold = args.auc_threshold

if args.percentile_threshold is not None:
    percentiles = derive_auc_threshold(ex_matrix_df)
    auc_threshold = percentiles[args.percentile_threshold]

aucs_mtx = aucell(ex_matrix_df, signatures, auc_threshold=auc_threshold, num_workers=args.num_workers)
aucs_mtx.to_csv(path_or_buf=args.output, index=True, sep='\t')

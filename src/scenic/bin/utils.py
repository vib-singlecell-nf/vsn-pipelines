#!/usr/bin/env python3

import glob
import ntpath
import os
import warnings
from typing import List
import gzip
import loompy as lp
import pandas as pd
from pyscenic.genesig import GeneSignature
from pyscenic.transform import COLUMN_NAME_CONTEXT, COLUMN_NAME_TARGET_GENES


def read_signatures_from_tsv_dir(dpath: str, noweights=False, weight_threshold=0, min_genes=0, show_warnings=False) -> List['GeneSignature']:
    """
    Load gene signatures from a list of .tsv files in directory. Requires .tsv with 1 or 2 columns. First column should be genes, second (optional) are weight for genes.
    :param dpath: The filepath to directory.
    :return: A list of signatures.
    """
    # https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
    assert os.path.exists(dpath), "{} does not exist.".format(dpath)
    gene_sig_file_paths = glob.glob(os.path.join(dpath, "*.tsv"))

    def signatures():
        for gene_sig_file_path in gene_sig_file_paths:
            gene_sig = pd.read_csv(gene_sig_file_path, sep='\t', header=None, index_col=None)
            fname = ntpath.basename(gene_sig_file_path)
            regulon = os.path.splitext(fname)[0]

            # Check if the file is the regulon frequency file
            if regulon == 'regulons':
                continue

            # Do some sanity checks
            if len(gene_sig.columns) == 0:
                raise Exception(f"{gene_sig_file_path} has 0 columns. Requires .tsv with 1 or 2 columns. First column should be genes (required), second (optional) are weight for the given genes.")
            if len(gene_sig.columns) > 2:
                raise Exception(f"{gene_sig_file_path} has more than 2 columns. Requires .tsv with 1 or 2 columns. First column should be genes, second (optional) are weight for the given genes.")
            if len(gene_sig.columns) == 1 or noweights:
                gene2weight = gene_sig[0]
            if len(gene_sig.columns) == 2 and not noweights:
                # Filter the genes based on the given weight_threshold
                # 1st column: genes
                # 2nd column: weights
                gene_sig = gene_sig[gene_sig[1] > weight_threshold]
                if len(gene_sig.index) == 0:
                    if show_warnings:
                        warnings.warn(
                            "{0} is empty after apply filter with weight_threshold > {1}".format(
                                regulon,
                                weight_threshold)
                        )
                    continue
                gene2weight = [tuple(x) for x in gene_sig.values]
            yield GeneSignature(name=regulon, gene2weight=gene2weight)

    signatures = list(signatures())
    # Filter regulons with less than min_genes (>= min_genes)
    signatures = list(filter(lambda x: len(x.gene2weight) >= min_genes, signatures))
    # Subtract 1 because regulons.tsv should not be counted
    print("Signatures passed filtering {0} out of {1}".format(len(signatures), len(gene_sig_file_paths) - 1))
    return signatures


def read_feature_enrichment_table(fname, sep, chunksize=None, raw=False):
    converters = None

    if not raw and chunksize is None:
        with gzip.open(fname, 'r') as fh:
            fh.readline()
            header_line2 = str(fh.readline()).rstrip('\n').split(',')
            # Get the index column of COLUMN_NAME_CONTEXT and COLUMN_NAME_TARGET_GENES
            column_name_context_idx = header_line2.index(COLUMN_NAME_CONTEXT)
            column_name_target_genes_idx = header_line2.index(COLUMN_NAME_TARGET_GENES)
            # Use converters to evaluate the expression a the given columns
            converters = {
                column_name_context_idx: eval,
                column_name_target_genes_idx: eval,
            }
    return pd.read_csv(
        fname,
        sep=sep,
        index_col=[0, 1],
        header=[0, 1],
        skipinitialspace=True,
        chunksize=chunksize,
        converters=converters,
        engine='c'
    )


def get_matrix(loom_file_path, gene_attribute, cell_id_attribute):
    with lp.connect(loom_file_path, mode='r', validate=False) as loom:  # Read in r mode otherwise concurrency problem
        ex_matrix = loom[:, :]
        ex_matrix_df = pd.DataFrame(
            data=ex_matrix[:, :],
            index=loom.ra[gene_attribute],
            columns=loom.ca[cell_id_attribute]
        ).T  # Gene expression as (cell, gene) - matrix.
    return ex_matrix_df

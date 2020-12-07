#!/usr/bin/env python3

import argparse
import base64
import json
import zlib

import os
import loompy as lp
import numpy as np
import pandas as pd

import export_to_loom

################################################################################
################################################################################

parser = argparse.ArgumentParser(description='Integrate output from pySCENIC motif- and track-based runs')

parser.add_argument(
    '--loom_motif',
    help='Loom file from pySCENIC motif run',
    required=True,
    default='pyscenic_motif.loom'
)

parser.add_argument(
    '--loom_track',
    help='Loom file from pySCENIC track run',
    required=True,
    default='pyscenic_track.loom'
)

parser.add_argument(
    '--motif-regulons-folder',
    type=str,
    dest="motif_regulons_folder",
    required=False,
    help='Folder containing the track-based regulon .tsv files.'
)

parser.add_argument(
    '--track-regulons-folder',
    type=str,
    dest="track_regulons_folder",
    required=False,
    help='Folder containing the track-based regulon .tsv files.'
)

parser.add_argument(
    '--min-regulon-occurrence',
    type=int,
    dest="min_regulon_occurrence",
    required=False,
    help='Filter out regulons that have occurred strictly less than the given threshold.',
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
    '--loom_output',
    help='Final loom file with pySCENIC motif and track results integrated',
    required=True,
    default='pyscenic.loom'
)

args = parser.parse_args()


################################################################################
################################################################################


def integrate_motif_track(args):

    ################################################################################
    # load data from loom
    ################################################################################

    # scenic motif output
    mtf_scope_loom = export_to_loom.SCopeLoom.read_loom(filename=args.loom_motif, tag='motif')
    # scenic track output
    trk_scope_loom = export_to_loom.SCopeLoom.read_loom(filename=args.loom_track, tag='track')

    # if too many (h5py header limitation), filter the regulons
    if args.motif_regulons_folder is not None and args.min_regulon_occurrence is not None:
        print(f"Filtering motif-based regulons that occur less than {args.min_regulon_occurrence} times...")
        regulon_count_table = pd.read_csv(
            os.path.join(args.motif_regulons_folder, "regulons.tsv"),
            sep="\t"
        )
        regulons = regulon_count_table[regulon_count_table["count"] >= args.min_regulon_occurrence]["regulon"].str.replace("(", "_(")
        mtf_scope_loom.set_regulon_filter(regulons=regulons)

    if args.track_regulons_folder is not None and args.min_regulon_occurrence is not None:
        print(f"Filtering track-based regulons that occur less than {args.min_regulon_occurrence} times...")
        regulon_count_table = pd.read_csv(
            os.path.join(args.track_regulons_folder, "regulons.tsv"),
            sep="\t"
        )
        regulons = regulon_count_table[regulon_count_table["count"] >= args.min_regulon_occurrence]["regulon"].str.replace("(", "_(")
        trk_scope_loom.set_regulon_filter(regulons=regulons)

    ################################################################################
    # merge scenic outputs
    ################################################################################

    # Merge regulon data from scenic track output into scenic motif output
    mtf_scope_loom.merge_regulon_data(scope_loom=trk_scope_loom)
    mtf_scope_loom.set_scenic_min_genes_regulon(min_genes_regulon=args.min_genes_regulon)
    mtf_scope_loom.set_scenic_min_regulon_gene_occurrence(min_regulon_gene_occurrence=args.min_regulon_gene_occurrence)
    mtf_scope_loom.add_metrics(metrics=["nUMI", "nGene", "Percent_mito"])
    mtf_scope_loom.export(out_fname=args.loom_output, save_embeddings=False)


if __name__ == "__main__":
    integrate_motif_track(args)

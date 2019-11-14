#!/usr/bin/env python3

import argparse
import base64
import json
import zlib

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

    ################################################################################
    # merge scenic outputs
    ################################################################################

    # Merge regulon data from scenic track output into scenic motif output
    mtf_scope_loom.merge_regulon_data(scope_loom=trk_scope_loom)
    mtf_scope_loom.add_metrics(metrics=["nUMI", "nGene", "Percent_mito"])
    mtf_scope_loom.export(out_fname=args.loom_output, save_embeddings=False)


if __name__ == "__main__":
    integrate_motif_track(args)

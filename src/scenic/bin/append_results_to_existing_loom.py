#!/usr/bin/env python3

import argparse
import base64
import json
import sys
import zlib
from copy import deepcopy
from shutil import copyfile

import loompy as lp
import numpy as np
import pandas as pd
import export_to_loom

################################################################################
################################################################################

parser = argparse.ArgumentParser(
    description='Integrate output from pySCENIC with SCope loom from pipeline best-practices path'
)
parser.add_argument(
    '--loom_scope',
    help='Loom file for SCope from pipeline best-practices path',
    required=True,
    default='.loom'
)
parser.add_argument(
    '--loom_scenic',
    help='Loom file from pySCENIC run',
    required=True,
    default='pyscenic_motif.loom'
)
parser.add_argument(
    '--loom_output',
    help='Final loom file with pySCENIC results integrated',
    required=True,
    default='pyscenic.loom'
)
args = parser.parse_args()


################################################################################
################################################################################

def append_to_existing_loom(args):

    ################################################################################
    # Merge looms
    ################################################################################

    scope_loom = export_to_loom.SCopeLoom.read_loom(filename=args.loom_scope)
    scenic_loom = export_to_loom.SCopeLoom.read_loom(filename=args.loom_scenic)
    # Make sure that both loom files have the same ordered feature space
    scope_loom.sort(axis=0, by="Gene")
    scenic_loom.sort(axis=0, by="Gene")
    scope_loom.merge(loom=scenic_loom)
    scope_loom.export(out_fname=args.loom_output)


if __name__ == "__main__":
    append_to_existing_loom(args)

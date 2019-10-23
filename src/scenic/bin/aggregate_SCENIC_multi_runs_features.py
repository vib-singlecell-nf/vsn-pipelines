#!/usr/bin/env python3

import argparse
import re
import sys

import pandas as pd
from pyscenic.transform import COLUMN_NAME_TYPE, COLUMN_NAME_CONTEXT
from pyscenic.utils import ACTIVATING_MODULE, REPRESSING_MODULE

import utils

################################################################################
################################################################################


parser_grn = argparse.ArgumentParser(
    usage="usage: %prog [options] feature_enrichment_table_fnames",
    description='Aggregate feature (motif or track) enrichment tables from multiple SCENIC runs. '
)

parser_grn.add_argument(
    'feature_enrichment_table_fnames',
    nargs='+',
    help='The looms resulting from the SCENIC AUCell step.'
)
parser_grn.add_argument(
    '-o', '--output',
    type=argparse.FileType('w'),
    default=sys.stdout,
    help='Output file/stream, i.e. a table of aggregated feature (motif or track) enrichment table (CSV).'
)

args = parser_grn.parse_args()


# Do stuff

def stack_feature_enrichment_tables(feature_enrichment_table_fnames):
    multi_runs_feature_enrichment_table = None

    def get_type(row):
        ctx = row[('Enrichment', COLUMN_NAME_CONTEXT)]
        # Activating is the default!
        return REPRESSING_MODULE if REPRESSING_MODULE in ctx else ACTIVATING_MODULE

    for i in range(0, len(feature_enrichment_table_fnames)):
        fname = feature_enrichment_table_fnames[i]
        feature_enrichment_table = utils.read_feature_enrichment_table(fname=fname, sep=",")
        # Add the run ID
        feature_enrichment_table[('', 'RunID')] = int(
            re.search(r'run_([0-9]+)+__reg_(mtf|trk)\.csv$', str(fname)).group(1)
        )
        feature_enrichment_table[('', COLUMN_NAME_TYPE)] = feature_enrichment_table.apply(get_type, axis=1)
        if multi_runs_feature_enrichment_table is None:
            multi_runs_feature_enrichment_table = feature_enrichment_table
        else:
            multi_runs_feature_enrichment_table = pd.concat(
                [multi_runs_feature_enrichment_table, feature_enrichment_table]
            )
        # print("Reading {0}... Done! ({1})".format(reg_path, multi_runs_feature_enrichment_table.shape))
    return multi_runs_feature_enrichment_table


multi_runs_feature_enrichment_table = stack_feature_enrichment_tables(args.feature_enrichment_table_fnames)
multi_runs_feature_enrichment_table.to_csv(path_or_buf=args.output)

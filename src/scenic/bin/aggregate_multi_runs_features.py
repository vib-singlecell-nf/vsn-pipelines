#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import pickle
from pyscenic.transform import COLUMN_NAME_TYPE, COLUMN_NAME_CONTEXT
from pyscenic.utils import ACTIVATING_MODULE, REPRESSING_MODULE
from pyscenic.transform import COLUMN_NAME_CONTEXT, COLUMN_NAME_TARGET_GENES
import re
import sys
import time
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
    help='Output file/stream, i.e. a table of aggregated feature (motif or track) enrichment table.'
)
parser_grn.add_argument(
    '-f', '--output-format',
    dest='output_format',
    type=str,
    default='csv',
    choices=['pickle', 'parquet', 'csv'],
    help='Output format of the aggregated feature (motif or track) enrichment table.'
)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


parser_grn.add_argument(
    "-b", '--use-chunking',
    type=str2bool,
    nargs='?',
    dest='use_chunking',
    const=True,
    default=True,
    help="Use the chunking method (low memory usage). Output format is fixed to compressed csv (.csv.gz). If this is set to False, you will require a server w/ a lot of RAM! (> 60gb)")

args = parser_grn.parse_args()

# Utils


def get_type(row):
    ctx = row[('Enrichment', COLUMN_NAME_CONTEXT)]
    # Activating is the default!
    return REPRESSING_MODULE if REPRESSING_MODULE in ctx else ACTIVATING_MODULE


# Do stuff


def save_aggregated_feature_enrichment_table(feature_enrichment_table_fnames, out_fname, chunksize=None):

    def to_csv(fname, df, count):
        mode = 'a' if count > 0 else 'w'
        header = False if count > 0 else True
        print(f"Saving {f_name} {df.shape}... Done!", flush=True)
        df.to_csv(fname, mode=mode, header=header)

    for f_name_idx in range(0, len(feature_enrichment_table_fnames)):
        f_name = feature_enrichment_table_fnames[f_name_idx]
        df = utils.read_feature_enrichment_table(
            fname=f_name,
            sep=",",
            chunksize=chunksize
        )

        for df_chunk in df:
            to_csv(
                fname=out_fname,
                df=add_run_id_to_feature_enrichment_table(
                    df=df_chunk,
                    fname=f_name
                ),
                count=f_name_idx
            )


def add_run_id_to_feature_enrichment_table(df, fname):
    # Add the run ID
    df[('', 'RunID')] = int(
        re.search(r'run_([0-9]+)+__reg_(mtf|trk)\.csv.gz$', str(fname)).group(1)
    )
    df[('', COLUMN_NAME_TYPE)] = df.apply(get_type, axis=1)
    return df


def read_and_add_run_id_to_feature_enrichment_table(fname):
    # Check if the file is empty
    if os.stat(fname).st_size == 0:
        raise Exception(f"The given file {fname} is empty. Please rerun the GRNboost and/or the cisTarget for this run.", flush=True)

    print("Reading {0}...".format(fname), flush=True)
    return add_run_id_to_feature_enrichment_table(
        df=utils.read_feature_enrichment_table(fname=fname, sep=","),
        fname=fname
    )


def stack_feature_enrichment_tables(feature_enrichment_table_fnames):
    return pd.concat((read_and_add_run_id_to_feature_enrichment_table(fname=f) for f in feature_enrichment_table_fnames))


# Stack all the enrichment tables from the different SCENIC runs
print("Stacking all the enrichment tables from the different runs...", flush=True)
start = time.time()

# Use this method if memory is an issue
if args.use_chunking:
    save_aggregated_feature_enrichment_table(
        feature_enrichment_table_fnames=args.feature_enrichment_table_fnames,
        out_fname=args.output.name,
        chunksize=50000
    )
else:
    multi_runs_feature_enrichment_table = stack_feature_enrichment_tables(args.feature_enrichment_table_fnames)
    print("Final aggregated enrichment table ({0})".format(multi_runs_feature_enrichment_table.shape), flush=True)
    print(f"... took {time.time() - start} seconds to run.", flush=True)

    # Save
    print("Saving...", flush=True)
    start = time.time()
    if args.output_format == 'pickle':
        print("to pickle (pandas)...", flush=True)
        multi_runs_feature_enrichment_table.to_pickle(path=args.output.name)
    # Do not support MultiIndex
    elif args.output_format == 'parquet':
        print("to parquet...", flush=True)
        multi_runs_feature_enrichment_table.to_parquet(fname=args.output.name)
    elif args.output_format == 'csv':
        print("to csv...", flush=True)
        multi_runs_feature_enrichment_table.to_csv(path_or_buf=args.output)
    else:
        raise Exception(f"The given output_format {args.output_format} has not been implemented")
print(f"... took {time.time() - start} seconds to run.", flush=True)
print("Done.", flush=True)

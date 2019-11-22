#!/usr/bin/env python3

import argparse
import os

import loompy as lp
import numpy as np
import pandas as pd

################################################################################
################################################################################

parser_grn = argparse.ArgumentParser(
    usage="usage: %prog [options] auc_looms",
    description='Aggregate genes by regulon from multiple SCENIC runs.'
)

parser_grn.add_argument(
    'auc_looms',
    nargs='+',
    help='The looms resulting from the SCENIC AUCell step.'
)
parser_grn.add_argument(
    '-o', '--output',
    help='Output file/stream, i.e. a folder containing regulons (CSV).'
)

args = parser_grn.parse_args()


# Do stuff


def from_loom_connection_regulon_incidence_matrix_to_long_df(loom):
    df = pd.DataFrame(loom.ra.Regulons)
    df.index = loom.ra.Gene
    df.index.name = 'Gene'
    df.columns.name = 'Regulon'
    df_ = df.stack().reset_index()
    df_.columns = ["gene", "regulon", "is_target"]
    df_ = df_.loc[(df_['is_target'] == 1)]
    return df_


def make_regulon_count_table(auc_looms):
    regulon_names_all_runs = None
    for i in range(0, len(auc_looms)):
        loom = lp.connect(filename=auc_looms[i], validate=False)
        df = from_loom_connection_regulon_incidence_matrix_to_long_df(loom=loom)
        df = pd.DataFrame(df['regulon'].unique(), columns=['regulon'])
        if regulon_names_all_runs is None:
            regulon_names_all_runs = df
        else:
            regulon_names_all_runs = pd.concat([regulon_names_all_runs, df], axis=0)
    return regulon_names_all_runs


def stack_regulons(auc_looms):
    all_runs_regulons_stacked = None

    for i in range(0, len(auc_looms)):
        loom = lp.connect(filename=auc_looms[i], validate=False)
        df = from_loom_connection_regulon_incidence_matrix_to_long_df(loom=loom)
        if all_runs_regulons_stacked is None:
            all_runs_regulons_stacked = df
        else:
            all_runs_regulons_stacked = pd.concat([all_runs_regulons_stacked, df], axis=0)
    return all_runs_regulons_stacked


def aggregate_genes_by_regulons(all_runs_regulons_stacked):
    all_runs_regulons_stacked_aggr = all_runs_regulons_stacked.groupby(['gene', 'regulon']).agg(['count']).sort_values(
        by=['regulon'],
        ascending=False
    )
    all_runs_regulons_stacked_aggr.columns = all_runs_regulons_stacked_aggr.columns.droplevel()
    return all_runs_regulons_stacked_aggr.reset_index()


def save_aggregated_regulons(all_runs_regulons_aggregated, output_dir):
    os.mkdir(output_dir)
    for regulon_name in np.unique(all_runs_regulons_aggregated["regulon"]):
        regulon_target_gene_occurrence_df = all_runs_regulons_aggregated[
            (all_runs_regulons_aggregated.regulon == regulon_name)
        ].sort_values(
            by=['count'],
            ascending=False
        )[["gene", "count"]]
        regulon_target_gene_occurrence_df.to_csv(
            path_or_buf=os.path.join(output_dir, regulon_name + ".tsv"),
            header=False,
            sep="\t",
            index=False
        )


# Aggregate all target genes for each regulon and save it
all_runs_regulons_stacked = stack_regulons(auc_looms=args.auc_looms)
all_runs_regulons_aggregated = aggregate_genes_by_regulons(all_runs_regulons_stacked=all_runs_regulons_stacked)
save_aggregated_regulons(all_runs_regulons_aggregated=all_runs_regulons_aggregated, output_dir=args.output)

# Make the regulon count table and save to regulons.tsv
all_runs_regulon_count_table = make_regulon_count_table(auc_looms=args.auc_looms)
all_runs_regulon_count_table['regulon'].value_counts().to_frame(
    name="count"
).to_csv(
    path_or_buf=os.path.join(args.output, "regulons.tsv"),
    sep="\t",
    header=True,
    index=True,
    index_label="regulon"
)

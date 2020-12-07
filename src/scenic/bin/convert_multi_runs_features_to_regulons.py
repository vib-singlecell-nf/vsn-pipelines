#!/usr/bin/env python3

import argparse
import re
import sys
import pandas as pd
import pickle
import gzip
from pyscenic import transform
from pyscenic.transform import COLUMN_NAME_NES
from pyscenic.utils import COLUMN_NAME_MOTIF_SIMILARITY_QVALUE, COLUMN_NAME_ORTHOLOGOUS_IDENTITY, \
    COLUMN_NAME_ANNOTATION
import time
import utils

################################################################################
# TODO:
# This implementation should be optimized:
# It's taking several hours to run
################################################################################

parser_grn = argparse.ArgumentParser(description='Transform aggregated motif enrichment table to regulons')

parser_grn.add_argument(
    'motif_enrichment_table_fname',
    type=argparse.FileType('r'),
    help='The name of the file that contains the motif enrichments.'
)
parser_grn.add_argument(
    'signatures_fname',
    help='The name of the folder containing the signatures as .tsv files.'
)

parser_grn.add_argument(
    '-o', '--output',
    type=argparse.FileType('w'),
    default=sys.stdout,
    help='Output file/stream, i.e. a pickle file containing regulons inferred from motif enrichment table.'
)
# parser_grn.add_argument(
#     '--min-genes-regulon',
#     type=int,
#     default=5,
#     dest="min_genes_regulon",
#     help='The threshold used for filtering the regulons based on the number of targets (default: {}).'.format(5)
# )
# parser_grn.add_argument(
#     '--min-regulon-gene-occurrence',
#     type=int,
#     default=5,
#     dest="min_regulon_gene_occurrence",
#     help='The threshold used for filtering the genes bases on their occurrence (default: {}).'.format(5)
# )

args = parser_grn.parse_args()

# Transform motif enrichment table (generated from the cisTarget step) to regulons
print(f"Reading aggregated motif enrichment table...", flush=True)
start = time.time()
f = args.motif_enrichment_table_fname.name
if f.endswith('.pkl') or f.endswith('.pkl.gz') or f.endswith('.pickle') or f.endswith('.pickle.gz'):
    motif_enrichment_table = pd.read_pickle(path=f)
elif f.endswith('.csv') or f.endswith('.csv.gz'):
    motif_enrichment_table = utils.read_feature_enrichment_table(fname=args.motif_enrichment_table_fname.name, sep=",")
else:
    raise Exception("VSN ERROR: The aggregated feature enrichment table is in the wrong format. Expecting .pickle or .csv formats.")
print(f"... took {time.time() - start} seconds to run.", flush=True)

print(f"Making the regulons...", flush=True)
start = time.time()
regulons = transform.df2regulons(
    df=motif_enrichment_table,
    save_columns=[
        COLUMN_NAME_NES,
        COLUMN_NAME_ORTHOLOGOUS_IDENTITY,
        COLUMN_NAME_MOTIF_SIMILARITY_QVALUE,
        COLUMN_NAME_ANNOTATION
    ]
)
print(f"{len(regulons)} regulons from df2regulons.")

# Read the signatures saved in out/multi_runs_regulons_[mtf|trk]
# Keep all regulons and targets (so that all can be visualized in SCope)
signatures = utils.read_signatures_from_tsv_dir(
    dpath=args.signatures_fname,
    noweights=False,
    weight_threshold=0,
    min_genes=0
)
print(f"{len(signatures)} all regulons from out/multi_runs_regulons_[mtf|trk].")

# Filter regulons (regulons from motifs enrichment table) by the filtered signatures
regulons = list(filter(lambda x: x.name in list(map(lambda x: x.name, signatures)), regulons))
# Add gene2occurrence from filtered signatures to regulons
regulons = list(
    map(
        lambda x:
        x.copy(gene2occurrence=list(filter(lambda y: y.name == x.name, signatures))[0].gene2weight), regulons
    )
)
print(f"{len(regulons)} final regulons.")
print(f"Saving...")
with gzip.open(args.output.name, 'wb') as file_handler:
    pickle.dump(regulons, file_handler)
print(f"... took {time.time() - start} seconds to run.", flush=True)

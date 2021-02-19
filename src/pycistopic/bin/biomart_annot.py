#!/usr/bin/env python3

import pybiomart as pbm
import argparse
import pickle


parser = argparse.ArgumentParser(description='Biomart annotation download')

parser.add_argument(
    "--biomart_dataset_name",
    type=str,
    required=True,
    help='Biomart dataset name, e.g. "hsapiens_gene_ensembl".'
)

args = parser.parse_args()

################################################################################

dataset = pbm.Dataset(name=args.biomart_dataset_name,  host='http://www.ensembl.org')
annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT', na=False)
annot = annot[~filter]
annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\S)', r'chr\1')
annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
annot = annot[annot.Transcript_type == 'protein_coding']

with open('biomart_annot.pickle', 'wb') as f:
    pickle.dump(annot, f)


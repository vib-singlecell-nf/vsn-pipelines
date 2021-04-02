#!/usr/bin/env python3

import argparse
import pybiomart as pbm

from pycisTopic.qc import computeQCStats

#import os
#import scanpy as sc

parser = argparse.ArgumentParser(description='Template script')

parser.add_argument(
    "--sampleId",
    type=str,
    required=True,
    help='Sample ID.'
)
parser.add_argument(
    "--fragments",
    type=argparse.FileType('r'),
    required=True,
    help='Input fragments file.'
)
parser.add_argument(
    "--regions",
    type=argparse.FileType('r'),
    required=True,
    help='Path to regions file.'
)
parser.add_argument(
    "--output",
    type=argparse.FileType('w'),
    help='Output file.'
)
parser.add_argument(
    "-s", "--seed",
    type=int,
    action="store",
    dest="seed",
    default=0,
    help="Use this integer seed for reproducibility."
)

args = parser.parse_args()


singularity run -B /ddn1 /ddn1/vol1/site_scratch/leuven/325/vsc32528/sif/vibsinglecellnf-pycistopic-0.1.img ipython

args.sampleId = 'sample_test'
args.fragments = '/ddn1/vol1/staging/leuven/stg_00002/lcb/cflerin/analysis/pbmc_atac/analysis/nextflow/test/out/fragments/VIB_1.sinto.fragments.tsv.gz'

dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://www.ensembl.org')
annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT')
annot = annot[~filter]
annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\S)', r'chr\1')
annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
annot = annot[annot.Transcript_type == 'protein_coding']

##################################################

fragments_dict = { 
    args.sampleId: args.fragments
    }
path_to_regions = {'Run_1':'/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/MLV/10x/exp/ih/20190425_NextSeq500_10x_scATAC/MLV__4aa2e0__Mouse_liver_ctrl/outs/peaks.bed',
                  'Run_2':'/staging/leuven/stg_00002/lcb/lcb_projects/MLV/cellranger_atac/NovaSeq6000_20200730/MLV__0d3236__liver_fresh_07_07_2020/outs/peaks.bed'}

metadata_bc_dict, profile_data_dict = computeQCStats(fragments_dict= fragments_dict,
                tss_annotation = annot,
                stats=['barcodeRankPlot', 'insertSizeDistribution', 'profileTSS', 'FRIP'],
                label_list = None,
                path_to_regions = path_to_regions,
                n_cpu = 5,
                valid_bc = None,
                n_frag = 50,
                n_bc = None,
                tss_flank_window = 2000,
                tss_window = 50,
                tss_minimum_signal_window = 100,
                tss_rolling_window = 10,
                remove_duplicates = True)







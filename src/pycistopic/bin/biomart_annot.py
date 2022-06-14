#!/usr/bin/env python3

import argparse
import os
import pickle
import sys

import pybiomart as pbm


def main():
    parser = argparse.ArgumentParser(description="Biomart gene annotation download.")

    parser.add_argument(
        "--biomart_dataset_name",
        type=str,
        required=True,
        help='Biomart dataset name, e.g. "hsapiens_gene_ensembl", '
        '"mmusculus_gene_ensembl", "dmelanogaster_gene_ensembl", ... .',
    )
    parser.add_argument(
        "--biomart_host",
        type=str,
        required=True,
        help='Biomart host address, e.g. "http://www.ensembl.org", '
        '"http://nov2020.archive.ensembl.org/", ... .',
    )

    args = parser.parse_args()

    # Skip retrieving annotation from biomart, if it was already done.
    if os.path.exists("biomart_annot.pickle"):
        sys.exit(0)

    dataset = pbm.Dataset(name=args.biomart_dataset_name, host=args.biomart_host)
    annot = dataset.query(
        attributes=[
            "chromosome_name",
            "transcription_start_site",
            "strand",
            "external_gene_name",
            "transcript_biotype",
        ]
    )

    # Rename columns.
    annot.columns = ["Chromosome", "Start", "Strand", "Gene", "Transcript_type"]

    # Convert objects in chromosome column to strings.
    annot["Chromosome"] = annot["Chromosome"].astype(str)

    # Only keep protein coding genes.
    annot = annot[annot.Transcript_type == "protein_coding"]

    # Only keep genes on normal chromosomes: (1-99, X, Y, 2L, 2R, 2L, 3R).
    filter_chroms = annot["Chromosome"].str.contains("^[0-9]{1,2}$|^[XY]$|^[23][LR]$")
    annot = annot[(filter_chroms)]

    # Add "chr" to the beginning of the chromosome names to make them UCSC compatible.
    annot["Chromosome"] = annot["Chromosome"].str.replace(r"(\b\S)", r"chr\1")

    with open("biomart_annot.pickle", "wb") as fh:
        pickle.dump(annot, fh)


if __name__ == "__main__":
    main()

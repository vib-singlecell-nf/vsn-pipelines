#!/usr/bin/env python3

import argparse
import os
import pickle
import sys

import pandas as pd
import pybiomart as pbm
from pycisTopic.qc import compute_qc_stats


def main():
    parser = argparse.ArgumentParser(description="Compute QC stats.")

    parser.add_argument(
        "--input_files",
        type=str,
        required=True,
        nargs="+",
        action="append",
        help='Input files in the form of "SampleId,fragments_filename,peaks_filename". '
        "Multiple inputs are possible.",
    )
    parser.add_argument(
        "--n_frag",
        type=int,
        default=50,
        help="Threshold on the number of fragments to keep for a barcode.",
    )
    parser.add_argument(
        "--tss_flank_window",
        type=int,
        default=2000,
        help="Flanking window around the TSS.",
    )
    parser.add_argument(
        "--tss_window",
        type=int,
        default=50,
        help="Window around the TSS used to count fragments in the TSS when "
        "calculating the TSS enrichment per barcode.",
    )
    parser.add_argument(
        "--tss_minimum_signal_window",
        type=int,
        default=100,
        help="Tail window use to normalize the TSS enrichment (average signal "
        "in the X bp in the extremes of the TSS window).",
    )
    parser.add_argument(
        "--tss_rolling_window",
        type=int,
        default=10,
        help="Rolling window used to smooth signal.",
    )
    parser.add_argument(
        "--min_norm",
        type=float,
        default=0.1,
        help="Minimum normalization score. If the average minimum signal value "
        "is below this value, this number is used to normalize the TSS signal. "
        "This approach penalizes cells with fewer reads.",
    )

    parser.add_argument(
        "--biomart_annot_pkl",
        type=str,
        required=True,
        help="Biomart annotations in pickle format.",
    )
    parser.add_argument(
        "--output_metadata_dir",
        type=str,
        required=True,
        help="Metadata output dir for metadata pickle files per sample.",
    )
    parser.add_argument(
        "--output_profile_data_dir",
        type=str,
        required=True,
        help="Profile data output dir for profile data picle files per sample.",
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="Number of threads to use."
    )

    args = parser.parse_args()

    fragments_dict = {x[0].split(",")[0]: x[0].split(",")[1] for x in args.input_files}
    path_to_regions = {x[0].split(",")[0]: x[0].split(",")[2] for x in args.input_files}

    # Load biomart annotations:
    with open(args.biomart_annot_pkl, "rb") as fh:
        annot = pickle.load(fh)

    try:
        metadata_bc_dict, profile_data_dict = compute_qc_stats(
            fragments_dict=fragments_dict,
            tss_annotation=annot,
            stats=[
                "barcode_rank_plot",
                "duplicate_rate",
                "insert_size_distribution",
                "profile_tss",
                "frip",
            ],
            label_list=None,
            path_to_regions=path_to_regions,
            n_cpu=args.threads,
            valid_bc=None,
            n_frag=args.n_frag,
            n_bc=None,
            tss_flank_window=args.tss_flank_window,
            tss_window=args.tss_window,
            tss_minimum_signal_window=args.tss_minimum_signal_window,
            tss_rolling_window=args.tss_rolling_window,
            min_norm=args.min_norm,
            remove_duplicates=True,
            # _temp_dir=
        )
    except FileExistsError as e:
        print(e)
        print(
            "For errors with /tmp/ray, use an alternate temporary directory and "
            "map this to '/tmp' within the container. For example in Singularity: "
            "'-B /alt/tmp/path:/tmp', or with Docker: '-v /alt/tmp/path:/tmp'. "
            "Set these mappings in the Singularity/Docker 'runOptions' parameter "
            "within your config file."
        )
        sys.exit(1)

    ## load bap results to use for duplicate rate (if we are using bap output):
    # f_bap_qc = os.path.join(os.path.dirname(args.fragments),args.sampleId+'.QCstats.csv')
    # if os.path.isfile(f_bap_qc) and all(metadata_bc_dict[args.sampleId]['Dupl_rate_bap'] == 0):
    #    bapqc = pd.read_csv(f_bap_qc, index_col=0)
    #    metadata_bc_dict[args.sampleId]['Dupl_rate'] = bapqc['duplicateProportion']

    # Write output for each sample to separate pickle files.

    for sample, metadata_bc_df in metadata_bc_dict.items():
        with open(f"{args.output_metadata_dir}/{sample}.metadata.pkl", "wb") as fh:
            pickle.dump(metadata_bc_df, fh)

    for sample, profile_data_df in profile_data_dict.items():
        with open(
            f"{args.output_profile_data_dir}/{sample}.profile_data.pkl", "wb"
        ) as fh:
            pickle.dump(profile_data_df, fh)


if __name__ == "__main__":
    main()

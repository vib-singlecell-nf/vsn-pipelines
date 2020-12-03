#!/usr/bin/env python3

import os
import argparse
import numpy as np
import scanpy as sc


def add_IO_arguments(parser):
    parser.add_argument(
        "input",
        type=argparse.FileType('r'),
        help='The path to the anndata file.'
    )
    parser.add_argument(
        "output",
        type=argparse.FileType('w'),
        help='The path to the anndata output file.',
    )
    return parser


def add_cell_filters(parser):
    parser.add_argument(
        "-s", "--cell-filter-strategy",
        type=str,
        action="store",
        dest="cell_filter_strategy",
        default="fixedthresholds",
        choices=["fixedthresholds", "adaptivethresholds"],
        help="""
Strategy used to filter out cells. The 'fixedthresholds' strategy will only apply the given thresholds on the QC metrics (n_counts, n_genes, percent_mito).
The 'adaptivethresholds' strategy will remove outliers based on the median absolute deviation (MAD) from the median value of each metric (n_counts, n_genes) across all cells.
A value is considered an outlier if it is more than (+/-)3 MADs away (both directions) from the median of a given QC metric (n_counts, n_genes).
This filter strategy is applied in the log space of the QC metrics (n_counts, n_genes).
This filter strategy will NOT be applied on the percent_mito metric. A simple filter defined by max_percent_mito threshold will be applied.
If min_n_counts and and min_n_genes are given, they will overrule the lower bound of each metric defined by (-)3MADs away from the median.
             """
    )
    parser.add_argument(
        "-c", "--min-n-counts",
        type=int,
        action="store",
        dest="min_n_counts",
        default=-1,
        help="Filter out cells with less than the minimum n of counts."
    )
    parser.add_argument(
        "-C", "--max-n-counts",
        type=int,
        action="store",
        dest="max_n_counts",
        default=-1,
        help="Filter out cells with more than the maximum n of counts."
    )
    parser.add_argument(
        "-g", "--min-n-genes",
        type=int,
        action="store",
        dest="min_n_genes",
        default=-1,
        help="Filter out cells with less than the minimum n of genes expressed."
    )
    parser.add_argument(
        "-G", "--max-n-genes",
        type=int,
        action="store",
        dest="max_n_genes",
        default=-1,
        help="Filter out cells with more than the maximum n of genes expressed."
    )
    parser.add_argument(
        "-M", "--max-percent-mito",
        type=float,
        action="store",  # optional because action defaults to "store"
        dest="max_percent_mito",
        default=-1,
        help="Filter out cells with more than the maximum percentage of mitochondrial genes expressed."
    )
    return parser


def add_gene_filters(parser):
    parser.add_argument(
        "--min-number-cells",
        type=int,
        action="store",
        dest="min_number_cells",
        default=-1,
        help="Filter out genes that are detected in less than the minimum number of cells."
    )
    return parser


################################################################################


def compute_qc_stats(adata, args):
    #
    # Compute number of counts
    #
    adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))
    #
    # Compute number of genes
    #
    # simply compute the number of genes per cell (computes 'n_genes' column)
    sc.pp.filter_cells(adata, min_genes=0)
    #
    # Compute fraction of mitochondrial genes
    #
    # mito and genes/counts cuts
    mito_genes = adata.var_names.str.contains('^MT-|^mt[-:]')
    # for each cell compute fraction of counts in mito genes vs. all genes
    adata.obs['percent_mito'] = np.ravel(np.sum(adata[:, mito_genes].X, axis=1)) / np.ravel(np.sum(adata.X, axis=1))
    #
    # Write filtering value into adata.uns
    #
    if args.cell_filter_strategy == "fixedthresholds":
        adata.uns['sc'] = {
            'scanpy': {
                'filter': {
                    'cellFilterStrategy': args.cell_filter_strategy,
                    'cellFilterMinNCounts': args.min_n_counts,
                    'cellFilterMaxNCounts': args.max_n_counts,
                    'cellFilterMinNGenes': args.min_n_genes,
                    'cellFilterMaxNGenes': args.max_n_genes,
                    'cellFilterMaxPercentMito': args.max_percent_mito,
                    'geneFilterMinNCells': args.min_number_cells,
                }
            }
        }
    elif args.cell_filter_strategy == "adaptivethresholds":
        ncounts_lower, ncounts_upper = get_qc_metric_adaptive_thresholds_by_mad(
            adata=adata,
            qc_metric_name="n_counts"
        )
        ngenes_lower, ngenes_upper = get_qc_metric_adaptive_thresholds_by_mad(
            adata=adata,
            qc_metric_name="n_genes"
        )
        adata.uns['sc'] = {
            'scanpy': {
                'filter': {
                    'cellFilterStrategy': args.cell_filter_strategy,
                    'cellFilterMinNCounts': np.exp(ncounts_lower),
                    'cellFilterMaxNCounts': np.exp(ncounts_upper),
                    'cellFilterMinNGenes': np.exp(ngenes_lower),
                    'cellFilterMaxNGenes': np.exp(ngenes_upper),
                    'cellFilterMaxPercentMito': args.max_percent_mito,
                    'geneFilterMinNCells': args.min_number_cells,
                }
            }
        }
    else:
        raise Exception("VSN ERROR: Not a valid cell filter strategy.")
    return adata


def get_qc_metric_adaptive_thresholds_by_mad(adata, qc_metric_name: str):
    lower_n_mads, upper_n_mads = 3, 3
    #
    # Define a lower limit on metric if set
    #
    lowlimit = np.log(vars(args).get(f"min_{qc_metric_name}")) if vars(args).get(f"min_{qc_metric_name}") > 0 else 0
    #
    # Work in log space
    #
    log_qc_metric_name = f"log_{qc_metric_name}"
    adata.obs[log_qc_metric_name] = np.log(adata.obs[qc_metric_name])
    #
    # Compute median and standard deviation of the given QC metric (defined by qc_metric_name) in the log space
    #
    log_mean, log_stdev = adata.obs[log_qc_metric_name].median(), adata.obs[log_qc_metric_name].std()
    lower, upper = log_mean - lower_n_mads * log_stdev, log_mean + upper_n_mads * log_stdev
    lower = np.max(
        [
            lowlimit,
            lower
        ]
    )
    return lower, upper


def cell_filter(adata, args):
    if args.cell_filter_strategy == "fixedthresholds":
        print("Applying cell filter fixed thresholds strategy...")
        #
        # Filter on min/Max number of counts
        #
        if args.min_n_counts > 0:
            print("Filtering cells by min. number of counts threshold...")
            print(f"Before filter: {adata.shape[0]}")
            adata = adata[adata.obs['n_counts'] > args.min_n_counts, :]
            print(f"After filter: {adata.shape[0]}")
        if args.max_n_counts > 0:
            print("Filtering cells by max. number of counts threshold...")
            print(f"Before filter: {adata.shape[0]}")
            adata = adata[adata.obs['n_counts'] < args.max_n_counts, :]
            print(f"After filter: {adata.shape[0]}")
        #
        # Filter on min/Max number of genes
        #
        if args.min_n_genes > 0:
            print("Filtering cells by min. number of genes threshold...")
            print(f"Before filter: {adata.shape[0]}")
            sc.pp.filter_cells(adata, min_genes=args.min_n_genes)
            print(f"After filter: {adata.shape[0]}")
        if args.max_n_genes > 0:
            print("Filtering cells by max. number of genes threshold...")
            print(f"Before filter: {adata.shape[0]}")
            adata = adata[adata.obs['n_genes'] < args.max_n_genes, :]
            print(f"After filter: {adata.shape[0]}")
        #
        # Filter on percentage of mitochondrial genes
        #
        if args.max_percent_mito > 0:
            adata = adata[adata.obs['percent_mito'] < args.max_percent_mito, :]
        return adata
    elif args.cell_filter_strategy == "adaptivethresholds":
        print("Applying cell filter adaptive thresholds strategy...")
        ncounts_lower, ncounts_upper = get_qc_metric_adaptive_thresholds_by_mad(
            adata=adata,
            qc_metric_name="n_counts"
        )
        ngenes_lower, ngenes_upper = get_qc_metric_adaptive_thresholds_by_mad(
            adata=adata,
            qc_metric_name="n_genes"
        )
        #
        # Apply the filters
        #
        # n_counts upper bound
        print(f"[ADAPTIVE THRESHOLD] Filtering cells by number of counts upper threshold ({np.exp(ncounts_upper)})...")
        print(f"Before filter: {adata.shape[0]}")
        adata = adata[adata.obs['log_n_counts'] <= ncounts_upper]
        print(f"After filter: {adata.shape[0]}")
        # n_counts lower bound
        print(f"[ADAPTIVE THRESHOLD] Filtering cells by number of counts lower threshold ({np.exp(ncounts_lower)})...")
        print(f"Before filter: {adata.shape[0]}")
        adata = adata[adata.obs['log_n_counts'] >= ncounts_lower]
        print(f"After filter: {adata.shape[0]}")
        # n_genes upper bound
        print(f"[ADAPTIVE THRESHOLD] Filtering cells by number of genes upper threshold ({np.exp(ngenes_upper)})...")
        print(f"Before filter: {adata.shape[0]}")
        adata = adata[adata.obs['log_n_genes'] <= ngenes_upper]
        print(f"After filter: {adata.shape[0]}")
        # n_genes lower bound
        print(f"[ADAPTIVE THRESHOLD] Filtering cells by number of genes lower threshold ({np.exp(ngenes_lower)})...")
        print(f"Before filter: {adata.shape[0]}")
        adata = adata[adata.obs['log_n_genes'] >= ngenes_lower]
        print(f"After filter: {adata.shape[0]}")
        # percent_mito upper bound
        if args.max_percent_mito > 0:
            print(f"[FIXED THRESHOLD] Filtering cells by percentage of mitochondrial content upper threshold ({args.max_percent_mito})...")
            print(f"Before filter: {adata.shape[0]}")
            adata = adata[adata.obs['percent_mito'] <= args.max_percent_mito]
            print(f"After filter: {adata.shape[0]}")
        return adata
    else:
        raise Exception("VSN ERROR: Not a valid cell filter strategy.")


def gene_filter(adata, args):
    #
    # Filter on min number of cells
    #
    if args.min_number_cells > -1:
        sc.pp.filter_genes(adata, min_cells=args.min_number_cells)
    return adata


def read_data(args):
    # I/O
    # Expects h5ad file
    try:
        adata = sc.read_h5ad(filename=args.input.name)
    except IOError:
        raise Exception("VSN ERROR: Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(args.input)[0]))
    return adata


def output_data(adata, FILE_PATH_OUT_BASENAME):
    # I/O
    adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

################################################################################
################################################################################


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute and filter on QC metrics from anndata object.")
    subparsers = parser.add_subparsers()
    # add_IO_arguments(parser)

    ##################################################
    # cell parser
    parser_cellfilter = subparsers.add_parser(
        'cellfilter',
        help='Apply cell-level filters.'
    )
    add_IO_arguments(parser_cellfilter)
    add_cell_filters(parser_cellfilter)
    parser_cellfilter.set_defaults(func=cell_filter)

    ##################################################
    # gene parser
    parser_genefilter = subparsers.add_parser(
        'genefilter',
        help='Apply gene-level filters.'
    )
    add_IO_arguments(parser_genefilter)
    add_gene_filters(parser_genefilter)
    parser_genefilter.set_defaults(func=gene_filter)

    ##################################################
    # compute parser (uses all parameters together)
    parent_parser = argparse.ArgumentParser(add_help=False)
    add_IO_arguments(parent_parser)
    add_cell_filters(parent_parser)
    add_gene_filters(parent_parser)
    parser_compute = subparsers.add_parser(
        'compute',
        parents=[parent_parser],
        help='Compute cell- and gene-level statistics without applying filters.'
    )
    parser_compute.set_defaults(func=compute_qc_stats)

    args = parser.parse_args()

    if not hasattr(args, 'func'):
        parser.print_help()
    else:
        # input
        adata = read_data(args)
        # filter/compute
        adata = args.func(adata, args)
        # output
        FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]
        output_data(adata, FILE_PATH_OUT_BASENAME)

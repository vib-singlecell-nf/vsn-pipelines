nextflow.preview.dsl=2

// include getBaseName from '../../utils/files.nf'

process SC__SCENIC__CISTARGET {
    cache 'deep'
    container params.container

    input:
    file filteredloom
    file 'adj.tsv'
    file featherDB
    file annotation
    val type

    output:
    file "reg_${type}.csv"

    """
    pyscenic ctx \
        adj.tsv \
        ${featherDB} \
        --annotations_fname ${annotation} \
        --expression_mtx_fname ${filteredloom} \
        --cell_id_attribute ${params.cell_id_attribute} \
        --gene_attribute ${params.gene_attribute} \
        --mode "dask_multiprocessing" \
        --output reg_${type}.csv \
        --num_workers ${params.numWorkers} \
    """
}

/* options to implement:

        // motif enrichment arguments:
        --rank_threshold RANK_THRESHOLD
        --auc_threshold AUC_THRESHOLD
        --nes_threshold NES_THRESHOLD

        // motif annotation arguments:
        --min_orthologous_identity MIN_ORTHOLOGOUS_IDENTITY
        --max_similarity_fdr MAX_SIMILARITY_FDR
        --annotations_fname ANNOTATIONS_FNAME


        // module generation arguments:
        --thresholds THRESHOLDS [THRESHOLDS ...]
        --top_n_targets TOP_N_TARGETS [TOP_N_TARGETS ...]
        --top_n_regulators TOP_N_REGULATORS [TOP_N_REGULATORS ...]
        --min_genes MIN_GENES
        --expression_mtx_fname EXPRESSION_MTX_FNAME
*/

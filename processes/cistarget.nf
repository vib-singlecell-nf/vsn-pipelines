nextflow.preview.dsl=2

// include getBaseName from '../../utils/files.nf'

process SC__SCENIC__CISTARGET {
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}/cistarget/${params.sc.scenic.numRuns > 1 ? "run_" + runId : ""}", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    maxForks params.sc.scenic.maxForks
    
    input:
    val runId
    file filteredloom
    file "${params.sc.scenic.numRuns > 1 ? "run_" + runId +"__adj.tsv" : "adj.tsv"}"
    file featherDB
    file annotation
    val type

    output:
    file "${params.sc.scenic.numRuns > 1 ? "run_" + runId +"__reg_" + type + ".csv" : "reg_" + type + ".csv"}"

    """
    pyscenic ctx \
        ${params.sc.scenic.numRuns > 1 ? "run_" + runId +"__adj.tsv" : "adj.tsv"} \
        ${featherDB} \
        --annotations_fname ${annotation} \
        --expression_mtx_fname ${filteredloom} \
        --cell_id_attribute ${params.sc.scenic.cell_id_attribute} \
        --gene_attribute ${params.sc.scenic.gene_attribute} \
        --mode "dask_multiprocessing" \
        --output ${params.sc.scenic.numRuns > 1 ? "run_" + runId +"__reg_" + type + ".csv" : "reg_" + type + ".csv"} \
        --num_workers ${params.sc.scenic.numWorkers} \
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

nextflow.preview.dsl=2

process CISTARGET {

    // Process will be submitted as job if params.sc.scenic.labels.processExecutor = 'qsub' (default)
    label "${params.sc.scenic.labels ? params.sc.scenic.labels.processExecutor : "local"}"
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}/${sampleId}/cistarget/${params.sc.scenic.numRuns > 1 ? "run_" + runId : ""}", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=${params.sc.scenic.cistarget.pmem} -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    maxForks params.sc.scenic.cistarget.maxForks
    
    input:
    tuple val(sampleId), file(filteredLoom), file("${params.sc.scenic.numRuns > 1 ? sampleId + "__run_" + runId +"__adj.tsv" : sampleId + "__adj.tsv"}"), val(runId)
    file featherDB
    file annotation
    val type

    output:
    tuple val(sampleId), file(filteredLoom), file("${params.sc.scenic.numRuns > 1 ? sampleId + "__run_" + runId +"__reg_" + type + ".csv" : sampleId + "__reg_" + type + ".csv"}"), val(runId)

    """
    pyscenic ctx \
        ${params.sc.scenic.numRuns > 1 ? sampleId + "__run_" + runId +"__adj.tsv" : sampleId + "__adj.tsv"} \
        ${featherDB} \
        --annotations_fname ${annotation} \
        --expression_mtx_fname ${filteredLoom} \
        --cell_id_attribute ${params.sc.scenic.cell_id_attribute} \
        --gene_attribute ${params.sc.scenic.gene_attribute} \
        --mode "dask_multiprocessing" \
        --output ${params.sc.scenic.numRuns > 1 ? sampleId + "__run_" + runId +"__reg_" + type + ".csv" : sampleId + "__reg_" + type + ".csv"} \
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


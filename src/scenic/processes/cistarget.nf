nextflow.preview.dsl=2

if(!params.containsKey("test")) {
	binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
	binDir = ""
}

toolParams = params.sc.scenic
processParams = params.sc.scenic.cistarget

process CISTARGET {

    // Process will be submitted as job if toolParams.labels.processExecutor = 'qsub' (default)
    label "${toolParams.labels ? toolParams.labels.processExecutor : "local"}"
    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/cistarget/${toolParams.numRuns > 1 ? "run_" + runId : ""}", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${toolParams.numWorkers} -l pmem=${processParams.pmem} -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    maxForks processParams.maxForks
    
    input:
    tuple val(sampleId), path(filteredLoom), path("${toolParams.numRuns > 1 ? sampleId + "__run_" + runId +"__adj.tsv" : sampleId + "__adj.tsv"}"), val(runId)
    file featherDB
    file annotation
    val type

    output:
    tuple val(sampleId), path(filteredLoom), path("${toolParams.numRuns > 1 ? sampleId + "__run_" + runId +"__reg_" + type + ".csv" : sampleId + "__reg_" + type + ".csv"}"), val(runId)

    script:
    """
    pyscenic ctx \
        ${toolParams.numRuns > 1 ? sampleId + "__run_" + runId +"__adj.tsv" : sampleId + "__adj.tsv"} \
        ${featherDB} \
        --annotations_fname ${annotation} \
        --expression_mtx_fname ${filteredLoom} \
        --cell_id_attribute ${toolParams.cell_id_attribute} \
        --gene_attribute ${toolParams.gene_attribute} \
        --mode "dask_multiprocessing" \
        --output ${toolParams.numRuns > 1 ? sampleId + "__run_" + runId +"__reg_" + type + ".csv" : sampleId + "__reg_" + type + ".csv"} \
        --num_workers ${toolParams.numWorkers} \
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


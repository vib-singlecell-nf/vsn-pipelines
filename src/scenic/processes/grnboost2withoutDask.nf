nextflow.preview.dsl=2

if(!params.containsKey("test")) {
	binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
	binDir = ""
}

toolParams = params.sc.scenic
processParams = params.sc.scenic.grn

process GRNBOOST2_WITHOUT_DASK {

    // Process will be submitted as job if toolParams.labels.processExecutor = 'qsub' (default)
    label "${toolParams.labels ? toolParams.labels.processExecutor : "local"}"
    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/grnboost2withoutDask/${toolParams.numRuns > 1 ? "run_" + runId : ""}", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${toolParams.numWorkers} -l pmem=${processParams.pmem} -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    maxForks processParams.maxForks

    input:
        tuple val(sampleId), path(filteredLoom), val(runId)
        file tfs

    output:
        tuple val(sampleId), path(filteredLoom), path("${toolParams.numRuns > 1 ? sampleId + "__run_" + runId +"__adj.tsv" : sampleId + "__adj.tsv"}"), val(runId)

    script:
        """
        ${binDir}grnboost2_without_dask.py \
            $filteredLoom \
            $tfs \
            --output ${toolParams.numRuns > 1 ? sampleId + "__run_" + runId +"__adj.tsv" : sampleId + "__adj.tsv"} \
            --num_workers ${toolParams.numWorkers} \
            --cell_id_attribute ${toolParams.cell_id_attribute} \
            --gene_attribute ${toolParams.gene_attribute}
        """

}

/* options to implement:
        --seed ${params.grn.seed} \

flag parameters not yet implemented:
        --transpose
*/

nextflow.preview.dsl=2

if(!params.containsKey("test")) {
	binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
	binDir = ""
}

toolParams = params.sc.scenic
processParams = params.sc.scenic.aucell

process AUCELL {

    // Process will be submitted as job if toolParams.labels.processExecutor = 'qsub' (degfault)
    label "${toolParams.labels ? toolParams.labels.processExecutor : "local"}"
    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/aucell/${"numRuns" in toolParams && toolParams.numRuns > 1 ? "run_" + runId : ""}", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${toolParams.numWorkers} -l pmem=${processParams.pmem} -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    maxForks processParams.maxForks

    input:
        tuple val(sampleId), path(filteredLoom), path(regulons), val(runId)
        val type

    output:
        tuple val(sampleId), path(filteredLoom), path("${outputFileName}"), val(runId)

    script:
        outputFileName = "numRuns" in toolParams && toolParams.numRuns > 1 ? sampleId + "__run_" + runId +"__auc_" + type + ".loom": sampleId + "__auc_" + type + ".loom"
        """
        pyscenic aucell \
            $filteredLoom \
            $regulons \
            -o "${outputFileName}" \
            --cell_id_attribute ${toolParams.cell_id_attribute} \
            --gene_attribute ${toolParams.gene_attribute} \
            --num_workers ${processParams.numWorkers}
        """

}

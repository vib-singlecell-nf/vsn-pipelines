nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scenic/bin/" : ""

toolParams = params.sc.scenic
processParams = params.sc.scenic.grn

process ARBORETO_WITH_MULTIPROCESSING {

    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/arboreto_with_multiprocessing/${"numRuns" in toolParams && toolParams.numRuns > 1 ? "run_" + runId : ""}", mode: 'link', overwrite: true
    label 'compute_resources__scenic_grn'

    input:
        tuple \
            val(sampleId), \
            path(filteredLoom), \
            val(runId)
        file tfs

    output:
        tuple val(sampleId), \
        file(filteredLoom), \
        file("${outputFileName}"), \
        val(runId)

    script:
        if(toolParams.numRuns > 2 && task.maxForks > 1 && task.executor == "local")
            throw new Exception("Running multi-runs SCENIC is quite computationally extensive. Please submit it as a job instead.")
        outputFileName = "numRuns" in toolParams && toolParams.numRuns > 1 ? sampleId + "__run_" + runId +"__adj.tsv" : sampleId + "__adj.tsv"
        seed = "numRuns" in toolParams && toolParams.numRuns > 1 ? (params.global.seed + runId) : params.global.seed
        """
        arboreto_with_multiprocessing.py \
            $filteredLoom \
            $tfs \
            --output ${outputFileName} \
            --num_workers ${task.cpus} \
            --cell_id_attribute ${toolParams.cell_id_attribute} \
            --gene_attribute ${toolParams.gene_attribute} \
            --method ${processParams.algorithm} \
            --seed ${seed}
        """

}

/* options to implement:
flag parameters not yet implemented:
        --transpose
*/

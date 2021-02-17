nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scenic/bin/" : ""

def toolParams = params.tools.scenic
def processParams = toolParams.grn

process ARBORETO_WITH_MULTIPROCESSING {

    cache 'deep'
    container toolParams.container
    //publishDir "${toolParams.scenicoutdir}/${sampleId}/arboreto_with_multiprocessing/${"numRuns" in toolParams && toolParams.numRuns > 1 ? "run_" + runId : ""}", mode: 'link', overwrite: true
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
            println("Running multi-runs SCENIC is quite computationally extensive. Consider submitting this as a job, or limit the number of parallel processes with 'maxForks'.")
        outputFileName = "numRuns" in toolParams && toolParams.numRuns > 1 ? sampleId + "__run_" + runId +"__adj.tsv.gz" : sampleId + "__adj.tsv.gz"
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

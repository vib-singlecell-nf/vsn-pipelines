nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scenic/bin/" : ""

def toolParams = params.tools.scenic
def processParams = toolParams.grn

process ADD_PEARSON_CORRELATION {

    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/arboreto_with_multiprocessing/${"numRuns" in toolParams && toolParams.numRuns > 1 ? "run_" + runId : ""}", mode: 'link', overwrite: true, pattern: '*__adj.tsv'
    label 'compute_resources__scenic_aucell'

    input:
        tuple \
            val(sampleId), \
            path(filteredLoom), \
            path(adj), \
            val(runId)

    output:
        tuple val(sampleId), \
            file(filteredLoom), \
            file("${outputFileName}"), \
            val(runId)

    script:
        if(toolParams.numRuns > 2 && task.maxForks > 1 && task.executor == "local")
            println("Running multi-runs SCENIC is quite computationally extensive. Consider submitting this as a job, or limit the number of parallel processes with 'maxForks'.")
        outputFileName = "numRuns" in toolParams && toolParams.numRuns > 1 ? sampleId + "__run_" + runId +"__adj.tsv" : sampleId + "__adj.tsv"
        seed = "numRuns" in toolParams && toolParams.numRuns > 1 ? (params.global.seed + runId) : params.global.seed
        """
        pyscenic add_cor \
            ${adj} \
            $filteredLoom \
            --output ${outputFileName} \
            --cell_id_attribute ${toolParams.cell_id_attribute} \
            --gene_attribute ${toolParams.gene_attribute} \
        """

}



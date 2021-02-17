nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scenic/bin/" : ""

def toolParams = params.tools.scenic
def processParams = toolParams.aucell

process AUCELL {

    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}/aucell/${"numRuns" in toolParams && toolParams.numRuns > 1 ? "run_" + runId : ""}", mode: 'link', overwrite: true
    label 'compute_resources__scenic_aucell'

    input:
        tuple \
           val(sampleId), \
           path(filteredLoom), \
           path(regulons), \
           val(runId)
        val type

    output:
        tuple \
           val(sampleId), \
           path(filteredLoom), \
           path("${outputFileName}"), \
           val(runId)

    script:
        if(toolParams.numRuns > 2 && task.maxForks > 1 && task.executor == "local")
            println("Running multi-runs SCENIC is quite computationally extensive. Consider submitting this as a job, or limit the number of parallel processes with 'maxForks'.")

        outputFileName = "numRuns" in toolParams && toolParams.numRuns > 1 ? 
            sampleId + "__run_" + runId +"__auc_" + type + ".loom": 
            sampleId + "__auc_" + type + ".loom"
        seed = "numRuns" in toolParams && toolParams.numRuns > 1 ? 
            (params.global.seed + runId) : 
            params.global.seed
        """
        export MKL_NUM_THREADS=1
		export NUMEXPR_NUM_THREADS=1
		export OMP_NUM_THREADS=1
        pyscenic aucell \
            $filteredLoom \
            $regulons \
            -o "${outputFileName}" \
            --cell_id_attribute ${toolParams.cell_id_attribute} \
            --gene_attribute ${toolParams.gene_attribute} \
            --seed ${seed} \
            --num_workers ${task.cpus}
        """

}


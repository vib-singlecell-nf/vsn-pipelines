nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

toolParams = params.sc.scenic

process AGGR_MULTI_RUNS_REGULONS {

    // Process will be submitted as job if toolParams.labels.processExecutor = 'qsub' (default)
    label "${toolParams.labels ? toolParams.labels.processExecutor : "local"}"
    cache 'deep'
    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${toolParams.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    tuple val(sampleId), path(f)
    val type

    output:
    tuple val(sampleId), path("multi_runs_regulons_${type}")

    """
    ${binDir}aggregate_multi_runs_regulons.py \
        ${f} \
        --output "multi_runs_regulons_${type}" \
    """

}

/* options to implement:
*/


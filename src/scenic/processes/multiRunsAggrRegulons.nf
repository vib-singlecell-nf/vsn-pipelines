nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process SC__SCENIC__MULTI_RUNS_AGGR_REGULONS {
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}", mode: 'copy'
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    file f
    val type

    output:
    file "multi_runs_regulons_${type}"

    """
    ${binDir}aggregate_SCENIC_multi_runs_regulons.py \
        ${f} \
        --regulon-type ${type} \
        --num-runs ${params.sc.scenic.numRuns} \
    """
}

/* options to implement:
*/


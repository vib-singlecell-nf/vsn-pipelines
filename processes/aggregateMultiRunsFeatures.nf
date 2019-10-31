nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process SC__SCENIC__AGGR_MULTI_RUNS_FEATURES {
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}/multi_runs_cistarget/", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    file f
    val type

    output:
    file "multi_runs_features_${type}.csv"

    """
    ${binDir}aggregate_SCENIC_multi_runs_features.py \
        ${f} \
        --output "multi_runs_features_${type}.csv"
    """
}

/* options to implement:
*/


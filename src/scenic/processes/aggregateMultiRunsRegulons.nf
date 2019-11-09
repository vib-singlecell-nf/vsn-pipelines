nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process SC__SCENIC__AGGR_MULTI_RUNS_REGULONS {
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    file f
    val type

    output:
    file "multi_runs_regulons_${type}"

    """
    ${binDir}aggregate_SCENIC_multi_runs_regulons.py \
        ${f} \
        --output "multi_runs_regulons_${type}" \
    """
}

/* options to implement:
*/


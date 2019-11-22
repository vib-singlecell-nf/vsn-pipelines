nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process PUBLISH_LOOM {
    
    clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
    publishDir "${params.sc.scenic.scenicoutdir}/${sampleId}", mode: 'link', overwrite: true, saveAs: { filename -> params.sc.scenic.scenicScopeOutputLoom }

    input:
    tuple val(sampleId), file(f)

    output:
    tuple val(sampleId), file(f)

    """
    """
}


process VISUALIZE {

    container params.sc.scenic.container
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
    tuple val(sampleId), file(inputLoom)

    output:
    tuple val(sampleId), file("scenic_visualize.loom")

    """
    ${binDir}add_visualization.py \
        --loom_input ${inputLoom} \
        --loom_output scenic_visualize.loom \
        --num_workers ${params.sc.scenic.numWorkers}
    """

}


process MERGE_MOTIF_TRACK_LOOMS {

    container params.sc.scenic.container
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    publishDir "${params.sc.scenic.scenicoutdir}/${sampleId}", mode: 'link', overwrite: true

    input:
    tuple val(sampleId), file(motifLoom), file(trackLoom)

    output:
    tuple val(sampleId), file(params.sc.scenic.scenicOutputLoom)

    """
    ${binDir}merge_motif_track_loom.py \
        --loom_motif ${motifLoom} \
        --loom_track ${trackLoom} \
        --loom_output ${params.sc.scenic.scenicOutputLoom}
    """

}

/* options to implement:
*/

process APPEND_SCENIC_LOOM {

    container params.sc.scenic.container
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    publishDir "${params.global.outdir}/loom", mode: 'link', overwrite: true

    input:
    tuple val(sampleId), file(scopeLoom), file(scenicLoom)

    output:
    tuple val(sampleId), file("${sampleId}.${params.sc.scenic.scenicScopeOutputLoom}")

    """
    ${binDir}append_results_to_existing_loom.py \
        --loom_scope ${scopeLoom} \
        --loom_scenic ${scenicLoom} \
        --loom_output ${sampleId}.${params.sc.scenic.scenicScopeOutputLoom}
    """

}

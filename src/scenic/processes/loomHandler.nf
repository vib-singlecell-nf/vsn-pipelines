nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
    binDir = ""
}

toolParams = params.sc.scenic

process PUBLISH_LOOM {

    clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
    publishDir "${toolParams.scenicoutdir}/${sampleId}", mode: 'link', overwrite: true, saveAs: { filename -> toolParams.scenicScopeOutputLoom }

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path(f)

    script:
        """
        """

}


process VISUALIZE {

    container toolParams.container
    clusterOptions "-l nodes=1:ppn=${toolParams.numWorkers} -l pmem=2gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple val(sampleId), path(inputLoom)

    output:
        tuple val(sampleId), path("scenic_visualize.loom")

    script:
        """
        ${binDir}add_visualization.py \
            --loom_input ${inputLoom} \
            --loom_output scenic_visualize.loom \
            --num_workers ${toolParams.numWorkers}
        """

}


process MERGE_MOTIF_TRACK_LOOMS {

    container toolParams.container
    clusterOptions "-l nodes=1:ppn=${toolParams.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    publishDir "${toolParams.scenicoutdir}/${sampleId}", mode: 'link', overwrite: true

    input:
        tuple val(sampleId), path(motifLoom), path(trackLoom)

    output:
        tuple val(sampleId), path(toolParams.scenicOutputLoom)

    script:
        toolParams = params.sc.scenic
        """
        ${binDir}merge_motif_track_loom.py \
            --loom_motif ${motifLoom} \
            --loom_track ${trackLoom} \
            --loom_output ${toolParams.scenicOutputLoom}
        """

}

/* options to implement:
*/

process APPEND_SCENIC_LOOM {

    container toolParams.container
    clusterOptions "-l nodes=1:ppn=${toolParams.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    publishDir "${params.global.outdir}/loom", mode: 'link', overwrite: true

    input:
        tuple val(sampleId), path(scopeLoom), path(scenicLoom)

    output:
        tuple val(sampleId), path("${sampleId}.${toolParams.scenicScopeOutputLoom}")

    script:
        toolParams = params.sc.scenic
        """
        ${binDir}append_results_to_existing_loom.py \
            --loom_scope ${scopeLoom} \
            --loom_scenic ${scenicLoom} \
            --loom_output ${sampleId}.${toolParams.scenicScopeOutputLoom}
        """

}

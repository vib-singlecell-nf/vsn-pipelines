nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scenic/bin/" : ""

toolParams = params.tools.scenic

process PUBLISH_LOOM {

    publishDir "${toolParams.scenicoutdir}/${sampleId}", mode: 'link', overwrite: true, saveAs: { filename -> toolParams.scenicScopeOutputLoom }
    label 'compute_resources__minimal'

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
    label 'compute_resources__scenic'

    input:
        tuple val(sampleId), path(inputLoom)

    output:
        tuple val(sampleId), path("scenic_visualize.loom")

    script:
        """
        ${binDir}add_visualization.py \
            --loom_input ${inputLoom} \
            --loom_output scenic_visualize.loom \
            --num_workers ${task.cpus}
        """

}


process MERGE_MOTIF_TRACK_LOOMS {

    container toolParams.container
    publishDir "${toolParams.scenicoutdir}/${sampleId}", mode: 'link', overwrite: true
    label 'compute_resources__scenic'

    input:
        tuple val(sampleId), path(motifLoom), path(trackLoom)

    output:
        tuple val(sampleId), path(toolParams.scenicOutputLoom)

    script:
        toolParams = params.tools.scenic
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
    publishDir "${params.global.outdir}/loom", mode: 'link', overwrite: true
    label 'compute_resources__scenic'

    input:
        tuple val(sampleId), path(scopeLoom), path(scenicLoom)

    output:
        tuple val(sampleId), path("${sampleId}.${toolParams.scenicScopeOutputLoom}")

    script:
        toolParams = params.tools.scenic
        """
        ${binDir}append_results_to_existing_loom.py \
            --loom_scope ${scopeLoom} \
            --loom_scenic ${scenicLoom} \
            --loom_output ${sampleId}.${toolParams.scenicScopeOutputLoom}
        """

}

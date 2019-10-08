nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process SC__SCENIC__MERGESCENICLOOMS {
    cache 'deep'
    container params.sc.scenic.container
    // publishDir "${params.sc.scenic.scenicoutdir}", mode: 'copy'

    input:
    file motifloom
    file trackloom

    output:
    file params.sc.scenic.scenicOutputLoom

    """
    ${binDir}merge_SCENIC_motif_track_loom.py \
        --loom_motif ${motifloom} \
        --loom_track ${trackloom} \
        --loom_output ${params.sc.scenic.scenicOutputLoom} \
        --num_workers ${params.sc.scenic.numWorkers} \
    """
}

/* options to implement:
*/

process SC__SCENIC__APPENDSCENICLOOM {
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}", mode: 'copy'

    input:
    file scopeloom
    file scenicloom

    output:
    file params.sc.scenic.scenicScopeOutputLoom

    """
    ${binDir}append_SCENIC_results_to_existing_loom.py \
        --loom_scope ${scopeloom} \
        --loom_scenic ${scenicloom} \
        --loom_output ${params.sc.scenic.scenicScopeOutputLoom} \
    """
}


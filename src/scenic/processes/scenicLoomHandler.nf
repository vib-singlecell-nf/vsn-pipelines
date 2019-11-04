nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process SC__SCENIC__PUBLISH_LOOM {
    
    publishDir "${params.sc.scenic.scenicoutdir}", mode: 'link', overwrite: true, saveAs: { filename -> params.sc.scenic.scenicOutputLoom }

    input:
    file f

    output:
    file f

    """
    """
}


process SC__SCENIC__VISUALIZE {
    cache 'deep'
    container params.sc.scenic.container
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
    file input_loom

    output:
    file "scenic_visualize.loom"

    """
    ${binDir}add_visualization.py \
        --loom_input ${input_loom} \
        --loom_output scenic_visualize.loom \
        --num_workers ${params.sc.scenic.numWorkers} \
    """
}


process SC__SCENIC__MERGE_MOTIF_TRACK_LOOMS {
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}", mode: 'copy'
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"

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


nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/scenic/bin/"
} else {
  binDir = ""
}

process PUBLISH_LOOM {
    
    clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
    publishDir "${params.sc.scenic.scenicoutdir}", mode: 'link', overwrite: true, saveAs: { filename -> params.sc.scenic.scenicScopeOutputLoom }

    input:
    file f

    output:
    file f

    """
    """
}


process VISUALIZE {
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


process MERGE_MOTIF_TRACK_LOOMS {
    cache 'deep'
    container params.sc.scenic.container
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    publishDir "${params.sc.scenic.scenicoutdir}", mode: 'link', overwrite: true

    input:
    file motifloom
    file trackloom

    output:
    file params.sc.scenic.scenicOutputLoom

    """
    ${binDir}merge_motif_track_loom.py \
        --loom_motif ${motifloom} \
        --loom_track ${trackloom} \
        --loom_output ${params.sc.scenic.scenicOutputLoom} \
    """
}

/* options to implement:
*/

process APPEND_SCENIC_LOOM {
    cache 'deep'
    container params.sc.scenic.container
    clusterOptions "-l nodes=1:ppn=${params.sc.scenic.numWorkers} -l pmem=2gb -l walltime=24:00:00 -A ${params.global.qsubaccount}"
    publishDir "${params.sc.scenic.scenicoutdir}", mode: 'link', overwrite: true

    input:
    file scopeloom
    file scenicloom

    output:
    file params.sc.scenic.scenicScopeOutputLoom

    """
    ${binDir}append_results_to_existing_loom.py \
        --loom_scope ${scopeloom} \
        --loom_scenic ${scenicloom} \
        --loom_output ${params.sc.scenic.scenicScopeOutputLoom} \
    """
}


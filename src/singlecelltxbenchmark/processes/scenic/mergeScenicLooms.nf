nextflow.preview.dsl=2

process SC__SCENIC__MERGESCENICLOOMS {
    cache 'deep'
    container params.container
    publishDir "${params.scenicoutdir}", mode: 'copy'

    input:
    file motifloom
    file trackloom

    output:
    file params.scenicoutputloom

    """
    merge_SCENIC_motif_track_loom.py \
        --loom_motif ${motifloom} \
        --loom_track ${trackloom} \
        --loom_output ${params.scenicoutputloom} \
        --num_workers ${params.numWorkers} \
    """
}

/* options to implement:
*/


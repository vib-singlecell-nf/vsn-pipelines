nextflow.preview.dsl=2

process SC__SCENIC__MERGESCENICLOOMS {
    cache 'deep'
    container params.sc.scenic.container
    publishDir "${params.sc.scenic.scenicoutdir}", mode: 'copy'

    input:
    file motifloom
    file trackloom

    output:
    file params.sc.scenic.scenicoutputloom

    """
    merge_SCENIC_motif_track_loom.py \
        --loom_motif ${motifloom} \
        --loom_track ${trackloom} \
        --loom_output ${params.sc.scenic.scenicoutputloom} \
        --num_workers ${params.sc.scenic.numWorkers} \
    """
}

/* options to implement:
*/


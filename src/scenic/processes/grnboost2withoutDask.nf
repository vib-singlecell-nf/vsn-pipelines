nextflow.preview.dsl=2

// include getBaseName from '../../utils/files.nf'


process SC__SCENIC__GRNBOOST2WITHOUTDASK {
    cache 'deep'
    container params.sc.scenic.container

    input:
    file filteredloom
    file tfs

    output:
    file 'adj.tsv'

    """
    grnboost2_without_dask.py \
        $filteredloom \
        $tfs \
        --output adj.tsv \
        --num_workers ${params.sc.scenic.numWorkers} \
        --cell_id_attribute ${params.sc.scenic.cell_id_attribute} \
        --gene_attribute ${params.sc.scenic.gene_attribute}
    """
}

/* options to implement:
        --seed ${params.grn.seed} \

flag parameters not yet implemented:
        --transpose
*/

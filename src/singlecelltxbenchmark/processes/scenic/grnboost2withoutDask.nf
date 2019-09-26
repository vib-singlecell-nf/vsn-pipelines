nextflow.preview.dsl=2

// include getBaseName from '../../utils/files.nf'

process SC__SCENIC__GRNBOOST2WITHOUTDASK {
    cache 'deep'
    container params.container

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
        --num_workers ${params.numWorkers} \
        --cell_id_attribute ${params.cell_id_attribute} \
        --gene_attribute ${params.gene_attribute}
    """
}

/* options to implement:
        --seed ${params.grn.seed} \

flag parameters not yet implemented:
        --transpose
*/

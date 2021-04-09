nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SEURAT__MARKER_GENES;
} from '../processes/marker_genes.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow DIFFERENTIAL_GENE_EXPRESSION {

    take:
        data
    
    main:
        SC__SEURAT__MARKER_GENES( data )

    emit:
        SC__SEURAT__MARKER_GENES.out
}
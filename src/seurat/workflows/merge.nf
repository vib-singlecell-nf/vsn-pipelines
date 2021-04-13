nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SEURAT__MERGE;
} from '../processes/merge.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow MERGE {

    take:
        data
    
    main:
        SC__SEURAT__MERGE( data )

    emit:
        SC__SEURAT__MERGE.out
}

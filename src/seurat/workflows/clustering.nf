nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SEURAT__CLUSTERING;
} from '../processes/cluster.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow CLUSTERING {

    take:
        data
    
    main:
        SC__SEURAT__CLUSTERING( data )

    emit:
        SC__SEURAT__CLUSTERING.out
}
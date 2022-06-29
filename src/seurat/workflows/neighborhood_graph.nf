nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SEURAT__NEIGHBORHOOD_GRAPH;
} from '../processes/neighborhood_graph.nf' params(params)

workflow NEIGHBORHOOD_GRAPH {

    take:
        data
    
    main:
        SC__SEURAT__NEIGHBORHOOD_GRAPH( data )
    
    emit:
        SC__SEURAT__NEIGHBORHOOD_GRAPH.out
}
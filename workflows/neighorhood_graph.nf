nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include SC__SCANPY__NEIGHBORHOOD_GRAPH from './../processes/neighborhood_graph' params(params)

//////////////////////////////////////////////////////

workflow NEIGHBORHOOD_GRAPH {

    take:
        data

    main:
        SC__SCANPY__NEIGHBORHOOD_GRAPH( data )

    emit:
        SC__SCANPY__NEIGHBORHOOD_GRAPH.out

}

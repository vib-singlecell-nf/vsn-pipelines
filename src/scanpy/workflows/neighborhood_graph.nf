nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include {
    SC__SCANPY__NEIGHBORHOOD_GRAPH;
} from './../processes/neighborhood_graph.nf' params(params)

//////////////////////////////////////////////////////

workflow NEIGHBORHOOD_GRAPH {

    take:
        data

    main:
        out = SC__SCANPY__NEIGHBORHOOD_GRAPH( data )

    emit:
        out

}

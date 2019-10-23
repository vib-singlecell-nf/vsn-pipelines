nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include SC__SCANPY__CLUSTERING from '../processes/cluster.nf' params(params.sc.scanpy.clustering + params.global + params)
include SC__SCANPY__MARKER_GENES from '../processes/marker_genes.nf' params(params.sc.scanpy.marker_genes + params.global + params)

//////////////////////////////////////////////////////

workflow CLUSTER_IDENTIFICATION {
    get:
        data
    main:
        SC__SCANPY__CLUSTERING( data )
        SC__SCANPY__MARKER_GENES( SC__SCANPY__CLUSTERING.out )
    emit:
        SC__SCANPY__MARKER_GENES.out
}

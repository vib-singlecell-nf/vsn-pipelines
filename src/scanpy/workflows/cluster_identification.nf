nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include SC__SCANPY__CLUSTERING from '../processes/cluster.nf' params(params.sc.scanpy.clustering + params.global + params)
include SC__SCANPY__MARKER_GENES from '../processes/marker_genes.nf' params(params.sc.scanpy.marker_genes + params.global + params)

// reporting:
include GENERATE_REPORT from './create_report.nf' params(params.sc.scanpy.feature_scaling + params)

//////////////////////////////////////////////////////

workflow CLUSTER_IDENTIFICATION {
    get:
        data
    main:
        SC__SCANPY__CLUSTERING( data )
        GENERATE_REPORT(
            SC__SCANPY__CLUSTERING.out,
            file(params.sc.scanpy.clustering.report_ipynb),
            "SC_clustering_report"
        )
        SC__SCANPY__MARKER_GENES( SC__SCANPY__CLUSTERING.out )
    emit:
        SC__SCANPY__MARKER_GENES.out
}


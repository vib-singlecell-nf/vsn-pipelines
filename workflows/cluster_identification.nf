nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include SC__SCANPY__CLUSTERING from '../processes/cluster.nf' params(params)
include SC__SCANPY__MARKER_GENES from '../processes/marker_genes.nf' params(params)

// reporting:
include GENERATE_REPORT from './create_report.nf' params(params)

//////////////////////////////////////////////////////

workflow CLUSTER_IDENTIFICATION {

    get:
        normalizedTransformedData
        data

    main:
        SC__SCANPY__CLUSTERING( data )
        report = GENERATE_REPORT(
            SC__SCANPY__CLUSTERING.out,
            file(workflow.projectDir + params.sc.scanpy.clustering.report_ipynb),
            "SC_clustering_report"
        )
        marker_genes = SC__SCANPY__MARKER_GENES(
            normalizedTransformedData.join(SC__SCANPY__CLUSTERING.out)
        )

    emit:
        marker_genes
        report

}

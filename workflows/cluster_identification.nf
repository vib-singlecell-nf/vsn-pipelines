nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include '../../utils/processes/utils.nf'

// scanpy:
include '../processes/cluster' params(params)
include '../processes/marker_genes.nf' params(params)

// reporting:
include GENERATE_REPORT from './create_report.nf' params(params)

workflow CLUSTER_IDENTIFICATION {

    take:
        normalizedTransformedData
        data
        tag

    main:
        // To run multiple clustering, we need at least 1 argument that is a list
        def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.sc.scanpy.clustering) )
        if(clusteringParams.isBenchmarkMode()) {
            // Run
            out = SC__SCANPY__BENCHMARK_CLUSTERING(
                data.map{ 
                    // Remove the runtimeParams
                    it -> tuple(it[0], it[1], it[2])
                }.combine(
                    // Add the runtimeParams
                    clusteringParams.$(tag)
                )
            )
        } else {
            // Run
            out = SC__SCANPY__CLUSTERING( data.map { it -> tuple(it[0], it[1]) } )
        }

        // Generate the report
        report = GENERATE_REPORT(
            "CLUSTERING",
            out,
            file(workflow.projectDir + params.sc.scanpy.clustering.report_ipynb),
            clusteringParams.isBenchmarkMode()
        )

        // Find marker genes for each of clustering
        if(clusteringParams.isBenchmarkMode()) {
            marker_genes = SC__SCANPY__BENCHMARK_MARKER_GENES(
                normalizedTransformedData.combine(out, by: 0)
            )
        } else {
            marker_genes = SC__SCANPY__MARKER_GENES(
                normalizedTransformedData.join(out)
            )
        }

    emit:
        marker_genes
        report

}

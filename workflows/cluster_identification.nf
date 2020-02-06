nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include SC__SCANPY__MULTI_CLUSTERING from '../processes/cluster.nf' params(params)
include SC__SCANPY__CLUSTERING from '../processes/cluster.nf' params(params)
include SC__SCANPY__MARKER_GENES from '../processes/marker_genes.nf' params(params)

// reporting:
include GENERATE_REPORT from './create_report.nf' params(params)

//////////////////////////////////////////////////////

workflow CLUSTER_IDENTIFICATION {

    take:
        normalizedTransformedData
        data

    main:
        def isMultiArgs = false
        // Properly define the arguments and handle if not present in config
        def scp = params.sc.scanpy.clustering
        def method = !scp.containsKey("clusteringMethod") ? "NULL" : scp.clusteringMethod
        def resolution = !scp.containsKey("resolution") ? "NULL" : scp.resolution

        // To run multiple clustering, we need at least 1 argument that is a list
        if(method instanceof List
            || resolution instanceof List) {
            // Set 
            isMultiArgs = true
            // Prepare arguments stream
            $method = Channel.from(method)
            $resolution = Channel.from(resolution)
            $args = $method
                .combine($resolution)
                .view()
            // Run
            out = SC__SCANPY__MULTI_CLUSTERING( 
                data.combine( $args )
            )
        } else {
            // Run
            out = SC__SCANPY__CLUSTERING( data )
        }
        // Generate the report
        report = GENERATE_REPORT(
            out,
            file(workflow.projectDir + params.sc.scanpy.clustering.report_ipynb),
            "SC_Scanpy_clustering_report",
            isMultiArgs
        )
        // Find marker genes for each of clustering
    //     marker_genes = SC__SCANPY__MARKER_GENES(
    //         normalizedTransformedData.join(out)
    //     )

    // emit:
    //     marker_genes
    //     report

}

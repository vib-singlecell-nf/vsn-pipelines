nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    clean;
} from '../../utils/processes/utils.nf' params(params)

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SCANPY__CLUSTERING;
    SC__SCANPY__CLUSTERING_PREFLIGHT_CHECKS;
    SC__SCANPY__CLUSTERING_PARAMS;
    SC__SCANPY__PARAM_EXPLORE_CLUSTERING;
} from '../processes/cluster' params(params)
include {
    SC__SCANPY__PARAM_EXPLORE_MARKER_GENES;
    SC__SCANPY__MARKER_GENES;
} from '../processes/marker_genes.nf' params(params)
// reporting:
include {
    GENERATE_REPORT
} from './create_report.nf' params(params)

workflow CLUSTER_IDENTIFICATION {

    take:
        normalizedTransformedData
        data
        tag

    main:
        // To run multiple clustering, we need at least 1 argument that is a list
        def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.tools.scanpy.clustering) )
        // Run sanity checks
        if(params.tools.scanpy.clustering?.preflight_checks) {
            $data = SC__SCANPY__CLUSTERING_PREFLIGHT_CHECKS( data.map { it -> tuple(it[0], it[1]) } )
        } else {
            $data = data
        }

        if(clusteringParams.isParameterExplorationModeOn()) {
            // Run
            out = SC__SCANPY__PARAM_EXPLORE_CLUSTERING(
                $data.map{ 
                    // Remove the runtimeParams
                    it -> tuple(it[0], it[1], it[2])
                }.combine(
                    // Add the runtimeParams
                    clusteringParams.$(tag)
                )
            )
        } else {
            // Run
            out = SC__SCANPY__CLUSTERING( $data.map { it -> tuple(it[0], it[1]) } )
        }

        // Generate the report
        report = GENERATE_REPORT(
            "CLUSTERING",
            out,
            file(workflow.projectDir + params.tools.scanpy.clustering.report_ipynb),
            clusteringParams.isParameterExplorationModeOn()
        )

        // Find marker genes for each of clustering
        if(clusteringParams.isParameterExplorationModeOn()) {
            marker_genes = SC__SCANPY__PARAM_EXPLORE_MARKER_GENES(
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

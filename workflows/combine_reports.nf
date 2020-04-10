nextflow.preview.dsl=2

include '../../utils/processes/utils.nf' params(params)
include '../processes/cluster.nf' params(params)

workflow COMBINE_REPORTS {

    take:
        samples
        workflow_config_report
        qc_filter_report
        hvg_selection_report
        dim_reduction_report
        cluster_report
    
    main:
        def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.sc.scanpy.clustering) )
        ipynbs = qc_filter_report.map {
            it -> tuple(it[0], it[1])
        }.mix(
            samples.combine(workflow_config_report),
            hvg_selection_report.map {
                it -> tuple(it[0], it[1])
            },
            dim_reduction_report.map {
                it -> tuple(it[0], it[1])
            }
        ).groupTuple().map {
            it -> tuple(it[0], *it[1])
        }.join(
            cluster_report,
            by: 0
        )

        if(!clusteringParams.isParameterExplorationModeOn()) {
            ipynbs = ipynbs.map {
                it -> tuple(it[0], it[1..it.size()-2], null)
            }
        } else {
            ipynbs = ipynbs.map {
                it -> tuple(it[0], it[1..it.size()-2], it[it.size()-1])
            }
        }
    emit:
        out = ipynbs
}
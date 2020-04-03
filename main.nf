nextflow.preview.dsl=2

import static groovy.json.JsonOutput.*

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include '../utils/workflows/utils.nf' params(params)
INIT()
include '../utils/processes/utils.nf' params(params)

include QC_FILTER from './workflows/qc_filter.nf' params(params)
include NORMALIZE_TRANSFORM from './workflows/normalize_transform.nf' params(params)
include HVG_SELECTION from './workflows/hvg_selection.nf' params(params)
include SC__SCANPY__REGRESS_OUT from './processes/regress_out.nf' params(params)
include NEIGHBORHOOD_GRAPH from './workflows/neighborhood_graph.nf' params(params)
include DIM_REDUCTION_PCA from './workflows/dim_reduction_pca.nf' params(params)
include './workflows/dim_reduction.nf' params(params)
include './processes/cluster.nf' params(params)
include CLUSTER_IDENTIFICATION from './workflows/cluster_identification.nf' params(params)

// reporting:
include UTILS__GENERATE_WORKFLOW_CONFIG_REPORT from '../utils/processes/reports.nf' params(params)
include SC__SCANPY__MERGE_REPORTS from './processes/reports.nf' params(params)
include SC__SCANPY__REPORT_TO_HTML from './processes/reports.nf' params(params)
include COMBINE_REPORTS from './workflows/combine_reports.nf' params(params)

workflow single_sample {

    take:
        data

    main:
        // Process the data
        out = params.sc.scanpy.containsKey("filter") ? QC_FILTER( data ).filtered : data
        out = params.sc.scanpy.containsKey("data_transformation") && 
            params.sc.scanpy.containsKey("normalization") ? NORMALIZE_TRANSFORM( out ) : out
        out = HVG_SELECTION( out )
        out = params.sc.scanpy.containsKey("regress_out") ? SC__SCANPY__REGRESS_OUT( out.scaled ) : out.scaled
        DIM_REDUCTION_PCA( out )
        NEIGHBORHOOD_GRAPH( DIM_REDUCTION_PCA.out )
        DIM_REDUCTION_TSNE_UMAP( NEIGHBORHOOD_GRAPH.out )
        CLUSTER_IDENTIFICATION(
            NORMALIZE_TRANSFORM.out,
            DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
            "No Batch Effect Correction"
        )

        // Reporting
        samples = data.map { it -> it[0] }.view()
        UTILS__GENERATE_WORKFLOW_CONFIG_REPORT(
            file(workflow.projectDir + params.utils.workflow_configuration.report_ipynb)
        )

        ipynbs = COMBINE_REPORTS(
            samples,
            UTILS__GENERATE_WORKFLOW_CONFIG_REPORT.out,
            QC_FILTER.out.report,
            HVG_SELECTION.out.report,
            DIM_REDUCTION_TSNE_UMAP.out.report,
            CLUSTER_IDENTIFICATION.out.report
        )

        def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.sc.scanpy.clustering) )
        merged_report = SC__SCANPY__MERGE_REPORTS(
            ipynbs,
            "merged_report",
            clusteringParams.isParameterExplorationModeOn()
        )
        SC__SCANPY__REPORT_TO_HTML(SC__SCANPY__MERGE_REPORTS.out)

    emit:
        filtered_data = params.sc.scanpy.containsKey("filter") ? QC_FILTER.out.filtered : Channel.empty()
        hvg_data = HVG_SELECTION.out.hvg
        dr_pca_data = DIM_REDUCTION_PCA.out.view()
        final_processed_data = CLUSTER_IDENTIFICATION.out.marker_genes
        reports = ipynbs
        merged_report

}

workflow single_sample_scrublet {

    take:
        data

    main:
        single_sample( data )
        

}


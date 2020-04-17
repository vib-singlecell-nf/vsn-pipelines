nextflow.preview.dsl=2

import static groovy.json.JsonOutput.*

//////////////////////////////////////////////////////
//  Import sub-workflows

include '../../utils/processes/utils.nf' params(params)
include SC__H5AD_TO_FILTERED_LOOM from '../../utils/processes/h5adToLoom.nf' params(params)
include FILE_CONVERTER from '../../utils/workflows/fileConverter.nf' params(params)

include QC_FILTER from './qc_filter.nf' params(params)
include NORMALIZE_TRANSFORM from './normalize_transform.nf' params(params)
include HVG_SELECTION from './hvg_selection.nf' params(params)
include SC__SCANPY__REGRESS_OUT from '../processes/regress_out.nf' params(params)
include NEIGHBORHOOD_GRAPH from './neighborhood_graph.nf' params(params)
include DIM_REDUCTION_PCA from './dim_reduction_pca.nf' params(params)
include './dim_reduction.nf' params(params)
include '../processes/cluster.nf' params(params)
include CLUSTER_IDENTIFICATION from './cluster_identification.nf' params(params)

// reporting:
include UTILS__GENERATE_WORKFLOW_CONFIG_REPORT from '../../utils/processes/reports.nf' params(params)
include SC__SCANPY__MERGE_REPORTS from '../processes/reports.nf' params(params)
include SC__SCANPY__REPORT_TO_HTML from '../processes/reports.nf' params(params)
include COMBINE_REPORTS from './combine_reports.nf' params(params)

workflow SINGLE_SAMPLE {

    take:
        // Expects (sampleId, h5ad)
        data

    main:
        // Process the data
        filtered = params.sc.scanpy.containsKey("filter") 
            ? QC_FILTER( data ).filtered : data
        transformed_normalized = params.sc.scanpy.containsKey("data_transformation") \
            && params.sc.scanpy.containsKey("normalization") 
            ? NORMALIZE_TRANSFORM( filtered ) : filtered
        hvg_selection = HVG_SELECTION( transformed_normalized )
        out = params.sc.scanpy.containsKey("regress_out") 
            ? SC__SCANPY__REGRESS_OUT( hvg_selection.scaled ) : hvg_selection.scaled
        DIM_REDUCTION_PCA( out )
        NEIGHBORHOOD_GRAPH( DIM_REDUCTION_PCA.out )
        DIM_REDUCTION_TSNE_UMAP( NEIGHBORHOOD_GRAPH.out )
        CLUSTER_IDENTIFICATION(
            transformed_normalized,
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
            params.sc.scanpy.containsKey("filter") ? QC_FILTER.out.report : Channel.empty(),
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

        // Conversion
        // Convert h5ad to X (here we choose: loom format)
        if(params.sc.scanpy.containsKey("filter")) {
            filtered_loom = SC__H5AD_TO_FILTERED_LOOM( QC_FILTER.out.filtered )
            // In parameter exploration mode, this automatically merge all the results into the resulting loom
            final_processed_scope_loom = FILE_CONVERTER(
                CLUSTER_IDENTIFICATION.out.marker_genes.groupTuple(),
                'loom',
                QC_FILTER.out.filtered
            )
        } else {
            // We expect here the data channel to contain a filtered h5ad
            filtered_loom = SC__H5AD_TO_FILTERED_LOOM( data )
            final_processed_scope_loom = FILE_CONVERTER(
                CLUSTER_IDENTIFICATION.out.marker_genes.groupTuple(),
                'loom',
                data
            )
        }

        // Publishing
        final_published_data = SC__PUBLISH_H5AD(
            CLUSTER_IDENTIFICATION.out.marker_genes.map { 
                it -> tuple(it[0], it[1], null)
            },
            params.global.project_name+".single_sample.output"
        )

    emit:
        filtered_data = params.sc.scanpy.containsKey("filter") ? QC_FILTER.out.filtered : Channel.empty()
        filtered_loom
        hvg_data = HVG_SELECTION.out.hvg
        dr_pca_data = DIM_REDUCTION_PCA.out
        final_processed_data = CLUSTER_IDENTIFICATION.out.marker_genes
        final_published_data
        final_processed_scope_loom
        reports = ipynbs
        merged_report

}

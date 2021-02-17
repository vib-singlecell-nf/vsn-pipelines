nextflow.enable.dsl=2

import static groovy.json.JsonOutput.*

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    clean;
} from '../../utils/processes/utils.nf' params(params)
include {
    UTILS__GENERATE_WORKFLOW_CONFIG_REPORT;
} from '../../utils/processes/reports.nf' params(params)
include {
    PUBLISH;
} from '../../utils/workflows/utils.nf' params(params)
include {
    FINALIZE;
} from '../../utils/workflows/finalize.nf' params(params)
include {
    SC__H5AD_TO_FILTERED_LOOM;
} from '../../utils/processes/h5adToLoom.nf' params(params)
include {
    FILE_CONVERTER;
} from '../../utils/workflows/fileConverter.nf' params(params)
include {
    FILTER_AND_ANNOTATE_AND_CLEAN;
} from '../../utils/workflows/filterAnnotateClean.nf' params(params)

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    QC_FILTER;
} from './qc_filter.nf' params(params)
include {
    NORMALIZE_TRANSFORM;
} from './normalize_transform.nf' params(params)
include {
    HVG_SELECTION;
} from './hvg_selection.nf' params(params)
include {
    NEIGHBORHOOD_GRAPH;
} from './neighborhood_graph.nf' params(params)
include {
    DIM_REDUCTION_PCA;
} from './dim_reduction_pca.nf' params(params)
include {
    DIM_REDUCTION_TSNE_UMAP;
} from './dim_reduction.nf' params(params)
include {
    SC__SCANPY__CLUSTERING_PARAMS;
} from '../processes/cluster.nf' params(params)
include {
    CLUSTER_IDENTIFICATION;
} from './cluster_identification.nf' params(params)

// reporting:
include {
    SC__SCANPY__MERGE_REPORTS;
    SC__SCANPY__REPORT_TO_HTML;
} from '../processes/reports.nf' params(params)
include {
    COMBINE_REPORTS;
} from './combine_reports.nf' params(params)


workflow SINGLE_SAMPLE {

    take:
        // Expects (sampleId, h5ad)
        data

    main:
        // Prefilter the data
        out = FILTER_AND_ANNOTATE_AND_CLEAN( data )

        def scanpyParams = params.getToolParams("scanpy")

        filtered = scanpyParams?.filter ? QC_FILTER( out ).filtered : out
        transformed_normalized = scanpyParams?.data_transformation && scanpyParams?.normalization
            ? NORMALIZE_TRANSFORM( filtered ) : filtered
        out = HVG_SELECTION( transformed_normalized )
        DIM_REDUCTION_PCA( out.scaled )
        NEIGHBORHOOD_GRAPH( DIM_REDUCTION_PCA.out )
        DIM_REDUCTION_TSNE_UMAP( NEIGHBORHOOD_GRAPH.out )
        CLUSTER_IDENTIFICATION(
            transformed_normalized,
            DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
            "No Batch Effect Correction"
        )


        // Reporting
        samples = data.map { it -> it[0] }
        UTILS__GENERATE_WORKFLOW_CONFIG_REPORT(
            file(workflow.projectDir + params.getUtilsParams("workflow_configuration").report_ipynb)
        )

        ipynbs = COMBINE_REPORTS(
            samples,
            UTILS__GENERATE_WORKFLOW_CONFIG_REPORT.out,
            scanpyParams?.filter ? QC_FILTER.out.report : Channel.empty(),
            HVG_SELECTION.out.report,
            DIM_REDUCTION_TSNE_UMAP.out.report,
            CLUSTER_IDENTIFICATION.out.report
        )

        def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(scanpyParams.clustering) )
        merged_report = SC__SCANPY__MERGE_REPORTS(
            ipynbs,
            "merged_report",
            clusteringParams.isParameterExplorationModeOn()
        )
        SC__SCANPY__REPORT_TO_HTML(SC__SCANPY__MERGE_REPORTS.out)



        // Finalize
        FINALIZE(
            scanpyParams?.filter ? QC_FILTER.out.filtered : data,
            CLUSTER_IDENTIFICATION.out.marker_genes,
            'SINGLE_SAMPLE.final_output'
        )

        marker_genes = CLUSTER_IDENTIFICATION.out.marker_genes.map { 
            it -> tuple(it[0], it[1], null)
        }

        // Publishing
        PUBLISH(
            marker_genes,
            params.global.project_name+".single_sample.final_output",
            "h5ad",
            null,
            clusteringParams.isParameterExplorationModeOn()
        )

    emit:
        filtered_data = scanpyParams?.filter ? QC_FILTER.out.filtered : Channel.empty()
        filtered_loom = FINALIZE.out.filteredloom
        hvg_data = HVG_SELECTION.out.hvg
        dr_pca_data = DIM_REDUCTION_PCA.out
        final_processed_data = marker_genes
        final_processed_scanpy_h5ad = FINALIZE.out.scanpyh5ad
        final_processed_scope_loom = FINALIZE.out.scopeloom
        reports = ipynbs
        merged_report

}

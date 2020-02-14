nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include '../src/utils/processes/utils.nf' params(params)

include star as STAR from '../workflows/star.nf' params(params)
include QC_FILTER from '../src/scanpy/workflows/qc_filter.nf' params(params)
include NORMALIZE_TRANSFORM from '../src/scanpy/workflows/normalize_transform.nf' params(params)
include HVG_SELECTION from '../src/scanpy/workflows/hvg_selection.nf' params(params)
include NEIGHBORHOOD_GRAPH from '../src/scanpy/workflows/neighborhood_graph.nf' params(params)
include DIM_REDUCTION_PCA from '../src/scanpy/workflows/dim_reduction_pca.nf' params(params)
include '../src/scanpy/workflows/dim_reduction.nf' params(params)
include '../src/scanpy/processes/cluster.nf' params(params)
include CLUSTER_IDENTIFICATION from '../src/scanpy/workflows/cluster_identification.nf' params(params)
include FILE_CONVERTER from '../src/utils/workflows/fileConverter.nf' params(params)

// data channel to start from 10x data:
include getChannel as getTenXChannel from '../src/channels/tenx.nf' params(params)

// reporting:
include UTILS__GENERATE_WORKFLOW_CONFIG_REPORT from '../src/utils/processes/reports.nf' params(params)
include SC__SCANPY__MERGE_REPORTS from '../src/scanpy/processes/reports.nf' params(params)
include SC__SCANPY__REPORT_TO_HTML from '../src/scanpy/processes/reports.nf' params(params)

workflow single_sample_star {
    
    // run the pipeline
    data = STAR()
    samples = data.map { it -> it[0] }
    UTILS__GENERATE_WORKFLOW_CONFIG_REPORT(
        file(workflow.projectDir + params.utils.workflow_configuration.report_ipynb)
    )
    QC_FILTER( data )
    NORMALIZE_TRANSFORM( QC_FILTER.out.filtered )
    HVG_SELECTION( NORMALIZE_TRANSFORM.out )
    DIM_REDUCTION_PCA( HVG_SELECTION.out.scaled )
    NEIGHBORHOOD_GRAPH( DIM_REDUCTION_PCA.out )
    DIM_REDUCTION_TSNE_UMAP( NEIGHBORHOOD_GRAPH.out )
    CLUSTER_IDENTIFICATION(
        NORMALIZE_TRANSFORM.out,
        DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
        "No Batch Effect Correction"
    )

    // Conversion
    // Convert h5ad to X (here we choose: loom format)
    filteredloom = SC__H5AD_TO_FILTERED_LOOM( QC_FILTER.out.filtered )
    scopeloom = FILE_CONVERTER(
        CLUSTER_IDENTIFICATION.out.marker_genes,
        'loom',
        QC_FILTER.out.filtered
    )

    // Publishing
    SC__PUBLISH_H5AD( 
        CLUSTER_IDENTIFICATION.out.marker_genes,
        "single_sample.output"
    )

    // Reporting:
    def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.sc.scanpy.clustering) )
    SC__SCANPY__MERGE_REPORTS(
        QC_FILTER.out.report.mix(
            samples.combine(UTILS__GENERATE_WORKFLOW_CONFIG_REPORT.out),
            HVG_SELECTION.out.report,
            DIM_REDUCTION_TSNE_UMAP.out.report,
            CLUSTER_IDENTIFICATION.out.report
        ).groupTuple(),
        "merged_report",
        clusteringParams.isBenchmarkMode()
    )
    SC__SCANPY__REPORT_TO_HTML(SC__SCANPY__MERGE_REPORTS.out)

    emit:
        filteredloom
        scopeloom

}

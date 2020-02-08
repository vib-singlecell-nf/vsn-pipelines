nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include QC_FILTER from '../src/scanpy/workflows/qc_filter.nf' params(params)
include SC__FILE_CONCATENATOR from '../src/utils/processes/utils.nf' params(params.sc.file_concatenator + params.global + params)
include NORMALIZE_TRANSFORM from '../src/scanpy/workflows/normalize_transform.nf' params(params + params.global)
include HVG_SELECTION from '../src/scanpy/workflows/hvg_selection.nf' params(params + params.global)
include DIM_REDUCTION from '../src/scanpy/workflows/dim_reduction.nf' params(params + params.global)
// CLUSTER_IDENTIFICATION
include '../src/scanpy/workflows/cluster_identification.nf' params(params + params.global) // Don't only import a specific process (the function needs also to be imported)
include SC__H5AD_TO_FILTERED_LOOM from '../src/utils/processes/h5adToLoom.nf' params(params + params.global)
include FILE_CONVERTER from '../src/utils/workflows/fileConverter.nf' params(params)
include BEC_BBKNN from '../src/scanpy/workflows/bec_bbknn.nf' params(params)

// data channel to start from 10x data:
include getChannel as getTenXChannel from '../src/channels/tenx.nf' params(params)

// reporting:
include UTILS__GENERATE_WORKFLOW_CONFIG_REPORT from '../src/utils/processes/reports.nf' params(params)
include SC__SCANPY__MERGE_REPORTS from '../src/scanpy/processes/reports.nf' params(params + params.global)
include SC__SCANPY__REPORT_TO_HTML from '../src/scanpy/processes/reports.nf' params(params + params.global)


workflow bbknn_base {

    take:
        data

    main:
        // run the pipeline
        QC_FILTER( data ) // Remove concat 
        SC__FILE_CONCATENATOR( QC_FILTER.out.filtered.map{it -> it[1]}.collect() )
        NORMALIZE_TRANSFORM( SC__FILE_CONCATENATOR.out )
        HVG_SELECTION( NORMALIZE_TRANSFORM.out )
        DIM_REDUCTION( HVG_SELECTION.out.scaled )

        //// Perform the clustering step w/o batch effect correction (for comparison matter)
        clusterIdentificationPreBatchEffectCorrection = CLUSTER_IDENTIFICATION( 
            NORMALIZE_TRANSFORM.out,
            DIM_REDUCTION.out.dimred_pca_tsne_umap,
            "Pre Batch Effect Correction"
        )

        //// Perform the batch effect correction
        BEC_BBKNN(
            NORMALIZE_TRANSFORM.out,
            //// include only PCA and t-SNE pre-merge dim reductions. Omit UMAP for clarity since it will have to be overwritten by BEC_BBKNN
            DIM_REDUCTION.out.dimred_pca_tsne,
            clusterIdentificationPreBatchEffectCorrection.marker_genes
        )

        // // conversion
        //// convert h5ad to X (here we choose: loom format)
        filteredloom = SC__H5AD_TO_FILTERED_LOOM( SC__FILE_CONCATENATOR.out )
        scopeloom = FILE_CONVERTER(
            BEC_BBKNN.out.data.groupTuple(),
            'loom',
            SC__FILE_CONCATENATOR.out
        )

        project = BEC_BBKNN.out.data.map { it -> it[0] }
        UTILS__GENERATE_WORKFLOW_CONFIG_REPORT(
            file(workflow.projectDir + params.utils.workflow_configuration.report_ipynb)
        )

        // collect the reports:
        ipynbs = project.combine(
            UTILS__GENERATE_WORKFLOW_CONFIG_REPORT.out
        ).join(
            HVG_SELECTION.out.report
        ).join(
            BEC_BBKNN.out.cluster_report
        ).combine(
            BEC_BBKNN.out.bbknn_report,
            by: 0
        ).map{ 
            tuple( it[0], it.drop(1) ) 
        }.view()
        // reporting:
        SC__SCANPY__MERGE_REPORTS(
            ipynbs,
            "merged_report",
            true
        )
        SC__SCANPY__REPORT_TO_HTML(SC__SCANPY__MERGE_REPORTS.out)

    emit:
        filteredloom
        scopeloom

}

workflow bbknn_standalone {

    main:
        data = getTenXChannel( params.data.tenx.cellranger_outs_dir_path ).view()
        bbknn_base( data )

    emit:
        filteredloom = bbknn_base.out.filteredloom
        scopeloom = bbknn_base.out.scopeloom

}

workflow bbknn {

    take:
        data

    main:
        bbknn_base( data )

    emit:
        filteredloom = bbknn_base.out.filteredloom
        scopeloom = bbknn_base.out.scopeloom

}

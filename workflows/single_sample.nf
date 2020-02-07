nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

//include CELLRANGER from '../src/cellranger/main.nf' params(params)
include QC_FILTER from '../src/scanpy/workflows/qc_filter.nf' params(params)
// include SC__FILE_CONCATENATOR from '../src/utils/processes/utils.nf' params(params)
include NORMALIZE_TRANSFORM from '../src/scanpy/workflows/normalize_transform.nf' params(params)
include HVG_SELECTION from '../src/scanpy/workflows/hvg_selection.nf' params(params)
include DIM_REDUCTION from '../src/scanpy/workflows/dim_reduction.nf' params(params)
include CLUSTER_IDENTIFICATION from '../src/scanpy/workflows/cluster_identification.nf' params(params)
include SC__H5AD_TO_FILTERED_LOOM from '../src/utils/processes/h5adToLoom.nf' params(params)
include FILE_CONVERTER from '../src/utils/workflows/fileConverter.nf' params(params)
include SC__PUBLISH_H5AD from '../src/utils/processes/utils.nf' params(params)

// data channel to start from 10x data:
include getChannel as getTenXChannel from '../src/channels/tenx.nf' params(params)

// reporting:
include UTILS__GENERATE_WORKFLOW_CONFIG_REPORT from '../src/utils/processes/reports.nf' params(params)
include SC__SCANPY__MERGE_REPORTS from '../src/scanpy/processes/reports.nf' params(params)
include SC__SCANPY__REPORT_TO_HTML from '../src/scanpy/processes/reports.nf' params(params)

workflow single_sample_base {

    take:
        data

    main:
        // run the pipeline
        samples = data.map { it -> it[0] }
        UTILS__GENERATE_WORKFLOW_CONFIG_REPORT(
            file(workflow.projectDir + params.utils.workflow_configuration.report_ipynb)
        )
        QC_FILTER( data )
        NORMALIZE_TRANSFORM( QC_FILTER.out.filtered )
        HVG_SELECTION( NORMALIZE_TRANSFORM.out )
        DIM_REDUCTION( HVG_SELECTION.out.scaled )
        CLUSTER_IDENTIFICATION(
            NORMALIZE_TRANSFORM.out,
            DIM_REDUCTION.out.dimred_pca_tsne_umap
        )

        // conversion
        //// convert h5ad to X (here we choose: loom format)
        filteredloom = SC__H5AD_TO_FILTERED_LOOM( QC_FILTER.out.filtered )

        // In benchmark mode, this automatically merge all the results into the resulting loom
        scopeloom = FILE_CONVERTER(
            CLUSTER_IDENTIFICATION.out.marker_genes.groupTuple(),
            'loom',
            QC_FILTER.out.filtered
        )

        // publishing
        SC__PUBLISH_H5AD( 
            CLUSTER_IDENTIFICATION.out.marker_genes,
            params.global.project_name+".single_sample.output"
        )

        // reporting:
        SC__SCANPY__MERGE_REPORTS(
            QC_FILTER.out.report.mix(
                samples.combine(UTILS__GENERATE_WORKFLOW_CONFIG_REPORT.out),
                HVG_SELECTION.out.report,
                DIM_REDUCTION.out.report,
                CLUSTER_IDENTIFICATION.out.report
            ).groupTuple(),
            "merged_report"
        )
        SC__SCANPY__REPORT_TO_HTML(SC__SCANPY__MERGE_REPORTS.out)

    emit:
        filteredloom
        scopeloom

}

workflow single_sample_standalone {

    main:
        data = getTenXChannel( params.data.tenx.cellranger_outs_dir_path ).view()
        single_sample_base( data )

    emit:
        filteredloom = single_sample_base.out.filteredloom
        scopeloom = single_sample_base.out.scopeloom

}

workflow single_sample {

    take:
        data

    main:
        single_sample_base( data )

    emit:
        single_sample_base.out.filteredloom
        single_sample_base.out.scopeloom

}

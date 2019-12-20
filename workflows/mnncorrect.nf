nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include QC_FILTER from '../src/scanpy/workflows/qc_filter.nf' params(params)
include SC__FILE_CONCATENATOR from '../src/utils/processes/utils.nf' params(params)
include NORMALIZE_TRANSFORM from '../src/scanpy/workflows/normalize_transform.nf' params(params)
include HVG_SELECTION from '../src/scanpy/workflows/hvg_selection.nf' params(params)
include SC__SCANPY__ADJUSTMENT from '../src/scanpy/processes/adjust.nf' params(params)
include BEC_MNN_CORRECT from '../src/scanpy/workflows/bec_mnn_correct.nf' params(params)
include SC__H5AD_TO_FILTERED_LOOM from '../src/utils/processes/h5adToLoom.nf' params(params)
include FILE_CONVERTER from '../src/utils/workflows/fileConverter.nf' params(params)

// data channel to start from 10x data:
include getChannel as getTenXChannel from '../src/channels/tenx.nf' params(params)

// reporting:
include UTILS__GENERATE_WORKFLOW_CONFIG_REPORT from '../src/utils/processes/reports.nf' params(params)
include SC__SCANPY__MERGE_REPORTS from '../src/scanpy/processes/reports.nf' params(params + params.global)
include SC__SCANPY__REPORT_TO_HTML from '../src/scanpy/processes/reports.nf' params(params + params.global)

workflow mnncorrect {

    // run the pipeline
    data = getTenXChannel( params.data.tenx.cellranger_outs_dir_path ).view()
    QC_FILTER( data ) // Remove concat
    SC__FILE_CONCATENATOR( QC_FILTER.out.filtered.map{it -> it[1]}.collect() )
    NORMALIZE_TRANSFORM( SC__FILE_CONCATENATOR.out )
    HVG_SELECTION( NORMALIZE_TRANSFORM.out )
    BEC_MNN_CORRECT(
        NORMALIZE_TRANSFORM.out,
        HVG_SELECTION.out.scaled
    )
    // SC__SCANPY__ADJUSTMENT( HVG_SELECTION.out.scaled )

    // conversion
    //// convert h5ad to X (here we choose: loom format)
    filteredloom = SC__H5AD_TO_FILTERED_LOOM( SC__FILE_CONCATENATOR.out )
    scopeloom = FILE_CONVERTER(
        BEC_MNN_CORRECT.out.data,
        'loom',
        SC__FILE_CONCATENATOR.out,
    )

    project = BEC_MNN_CORRECT.out.data.map { it -> it[0] }
    UTILS__GENERATE_WORKFLOW_CONFIG_REPORT(
        file(workflow.projectDir + params.utils.workflow_configuration.report_ipynb)
    )

    // collect the reports:
    ipynbs = project.combine(UTILS__GENERATE_WORKFLOW_CONFIG_REPORT.out)
        .join(HVG_SELECTION.out.report)
        .join(BEC_MNN_CORRECT.out.cluster_report)
        .join(BEC_MNN_CORRECT.out.bbknn_report)
        .map{ tuple( it[0], it.drop(1) ) }
    // reporting:
    SC__SCANPY__MERGE_REPORTS(
        ipynbs,
        "merged_report"
    )
    SC__SCANPY__REPORT_TO_HTML(SC__SCANPY__MERGE_REPORTS.out)

    emit:
        filteredloom
        scopeloom
}

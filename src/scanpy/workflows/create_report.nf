//
// Version: 
// Test: 
// Command: 
//
/*
 * QC workflow 
 * Source:
 * 
 * Steps considered: 
 * - filter (cell, gene) + qc report
 */ 

nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include '../processes/reports.nf' params(params.sc.scanpy.filter + params.global + params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow GENERATE_QC_REPORT {
    get:
        prefiltered  // anndata
        postfiltered // anndata
        ipynb
        report_title
    main:
        report_notebook = SC__SCANPY__FILTER_QC_REPORT(
            ipynb, prefiltered, postfiltered, report_title )
        SC__SCANPY__REPORT_TO_HTML( report_notebook, report_title )
    emit:
        report_notebook
}

workflow GENERATE_REPORT {
    get:
        data // anndata
        ipynb
        report_title
    main:
        report_notebook = SC__SCANPY__GENERATE_REPORT(
            ipynb, data, report_title )
        SC__SCANPY__REPORT_TO_HTML( report_notebook, report_title )
    emit:
        report_notebook
}


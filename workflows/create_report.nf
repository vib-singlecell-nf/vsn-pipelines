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

include '../processes/reports.nf' params(params.global + params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow GENERATE_DUAL_INPUT_REPORT {
    get:
        data1  // anndata
        data2 // anndata
        ipynb
        report_title
    main:
        report_notebook = SC__SCANPY__GENERATE_DUAL_INPUT_REPORT(
            ipynb,
            data1.join(data2),
            report_title
        )
        SC__SCANPY__REPORT_TO_HTML(report_notebook)
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
            ipynb,
            data,
            report_title
        )
        SC__SCANPY__REPORT_TO_HTML(report_notebook)
    emit:
        report_notebook
}


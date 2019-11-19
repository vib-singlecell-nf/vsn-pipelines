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

// utils:
include SC__FILE_CONVERTER from '../../utils/processes/utils.nf' params(params.sc.file_converter + params.global + params)
include SC__FILE_ANNOTATOR from '../../utils/processes/utils.nf' params(params.sc.file_annotator + params.global + params)

// scanpy:
include '../processes/filter.nf' params(params.sc.scanpy.filter + params.global + params)

// reporting:
include GENERATE_QC_REPORT from './create_report.nf' params(params.sc.scanpy.filter + params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow QC_FILTER {
    get:
        data
    main:
        data = SC__FILE_CONVERTER( data )
        if (params.sc.file_annotator.metaDataFilePath && params.sc.file_annotator.metaDataFilePath != '') {
            data = SC__FILE_ANNOTATOR( SC__FILE_CONVERTER.out, file(params.sc.file_annotator.metaDataFilePath) )
        }
        unfiltered = SC__SCANPY__COMPUTE_QC_STATS( data )
        SC__SCANPY__GENE_FILTER( unfiltered )
        filtered = SC__SCANPY__CELL_FILTER( SC__SCANPY__GENE_FILTER.out )
        report = GENERATE_QC_REPORT( 
            unfiltered, 
            filtered,
            file(workflow.projectDir + params.sc.scanpy.filter.report_ipynb),
            'SC_QC_filtering_report'
        )
    emit:
        filtered
        report
}

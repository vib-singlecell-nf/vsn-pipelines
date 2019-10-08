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
include SC__FILE_CONVERTER from '../utils/processes/utils.nf' params(params.sc.file_converter + params.global + params)
include SC__FILE_ANNOTATOR from '../utils/processes/utils.nf' params(params.sc.file_annotator + params.global + params)

// scanpy:
include './processes/filter.nf' params(params.sc.scanpy.filter + params.global + params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow QC_FILTER {
    get:
        cellrangercounts
    main:
        SC__FILE_CONVERTER( cellrangercounts )
        SC__FILE_ANNOTATOR( SC__FILE_CONVERTER.out, file(params.sc.file_annotator.metaDataFilePath) )
        SC__SCANPY__GENE_FILTER( SC__FILE_ANNOTATOR.out )
        filtered = SC__SCANPY__CELL_FILTER( SC__SCANPY__GENE_FILTER.out )
        report = SC__SCANPY__PREPARE_FILTER_QC_REPORT()
        SC__SCANPY__FILTER_QC_REPORT( report, filtered )
    emit:
        filtered
}


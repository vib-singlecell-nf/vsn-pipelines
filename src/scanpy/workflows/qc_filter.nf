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
include SC__FILE_CONVERTER from '../../utils/processes/utils.nf' params(params)
include FILE_CONVERTER from '../../utils/workflows/fileConverter.nf' params(params)
include SC__ANNOTATE_BY_SAMPLE_METADATA from '../../utils/processes/h5adAnnotate.nf' params(params)
include SC__ANNOTATE_BY_CELL_METADATA from '../../utils/processes/h5adAnnotate' params(params)
include FILTER_BY_CELL_METADATA from '../../utils/workflows/filterByCellMetadata.nf' params(params)

// scanpy:
include '../processes/filter.nf' params(params)

// reporting:
include GENERATE_DUAL_INPUT_REPORT from './create_report.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow QC_FILTER {

    get:
        data

    main:
        data = FILE_CONVERTER(
            data,
            'h5ad',
            Channel.empty()
        )
        if(params.sc.cell_filter) {
            data = FILTER_BY_CELL_METADATA( data )
        }
        if(params.sc.cell_annotate) {
            data = SC__ANNOTATE_BY_CELL_METADATA( data )
        }
        if (params.sc.sample_annotate
            && params.sc.sample_annotate.metaDataFilePath
            && params.sc.sample_annotate.metaDataFilePath != ''
        ) {
            data = SC__ANNOTATE_BY_SAMPLE_METADATA( data )
        }
        unfiltered = SC__SCANPY__COMPUTE_QC_STATS( data )
        SC__SCANPY__GENE_FILTER( unfiltered )
        filtered = SC__SCANPY__CELL_FILTER( SC__SCANPY__GENE_FILTER.out )
        report = GENERATE_DUAL_INPUT_REPORT(
            unfiltered, 
            filtered,
            file(workflow.projectDir + params.sc.scanpy.filter.report_ipynb),
            'SC_QC_filtering_report'
        )

    emit:
        filtered
        report

}

nextflow.preview.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    FILTER_AND_ANNOTATE_AND_CLEAN
} from '../../utils/workflows/filterAnnotateClean.nf' params(params)

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SCANPY__COMPUTE_QC_STATS;
    SC__SCANPY__CELL_FILTER;
    SC__SCANPY__GENE_FILTER
} from '../processes/filter.nf' params(params)

// reporting:
include {
    GENERATE_DUAL_INPUT_REPORT
} from './create_report.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow QC_FILTER {

    take:
        data

    main:
        FILTER_AND_ANNOTATE_AND_CLEAN( data )
        unfiltered = SC__SCANPY__COMPUTE_QC_STATS( FILTER_AND_ANNOTATE_AND_CLEAN.out )
        SC__SCANPY__CELL_FILTER( unfiltered )
        filtered = SC__SCANPY__GENE_FILTER( SC__SCANPY__CELL_FILTER.out )
        report = GENERATE_DUAL_INPUT_REPORT(
            unfiltered.join(filtered).map { 
                it -> tuple(*it[0..(it.size()-1)], null)
            },
            file(workflow.projectDir + params.sc.scanpy.filter.report_ipynb),
            'SC_QC_filtering_report',
            false
        )

    emit:
        filtered
        report

}

nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:

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
        filtered = data | \
            SC__SCANPY__COMPUTE_QC_STATS | \
            SC__SCANPY__CELL_FILTER | \
            SC__SCANPY__GENE_FILTER
        
        report = !params.tools.scanpy.filter?.report_ipynb ? Channel.empty() :
            GENERATE_DUAL_INPUT_REPORT(
                SC__SCANPY__COMPUTE_QC_STATS.out.join(filtered).map { 
                    it -> tuple(*it[0..(it.size()-1)], null)
                },
                file(workflow.projectDir + params.tools.scanpy.filter.report_ipynb),
                'SC_QC_filtering_report',
                false
            )

    emit:
        filtered
        report

}

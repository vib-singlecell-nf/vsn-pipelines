nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SEURAT__CLUSTERING;
} from '../processes/cluster.nf' params(params)
include {
    GENERATE_REPORT;
} from './create_report.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow CLUSTERING {

    take:
        data
    
    main:
        clustered = SC__SEURAT__CLUSTERING( data )
        report = GENERATE_REPORT(
            "CLUSTERING",
            clustered,
            file(workflow.projectDir + params.tools.seurat.clustering.report_rmd)
        )

    emit:
        clustered
        report
}
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
        SC__SEURAT__CLUSTERING( data )
        report = GENERATE_REPORT(
            "CLUSTERING",
            SC__SEURAT__CLUSTERING.out,
            file(workflow.projectDir + params.tools.seurat.clustering.report_rmd)
        )

    emit:
        SC__SEURAT__CLUSTERING.out
}
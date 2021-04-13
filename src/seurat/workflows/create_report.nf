nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include {
    SC__SEURAT__GENERATE_REPORT;
    SC__SEURAT__REPORT_TO_HTML
} from '../processes/reports.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow GENERATE_REPORT {

    take:
        pipelineStep
        data
        rmd

    main:
        def reportTitle = "SC_Seurat_" + pipelineStep.toLowerCase() + "_report"

        SC__SEURAT__GENERATE_REPORT(
            rmd,
            data,
            reportTitle
        )
        // FIXME: the only reason we pass data to the html converter is to make sure the rmd can find the rds file it needs to load
        // This can probably be improved!
        SC__SEURAT__REPORT_TO_HTML(
            SC__SEURAT__GENERATE_REPORT.out,
            data,
            reportTitle
        )
        
    emit:
        SC__SEURAT__GENERATE_REPORT.out

}

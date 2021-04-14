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
        ) | SC__SEURAT__REPORT_TO_HTML
        
    emit:
        SC__SEURAT__GENERATE_REPORT.out

}

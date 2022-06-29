nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include {
    SC__SEURAT__GENERATE_REPORT;
    SC__SEURAT__GENERATE_DUAL_INPUT_REPORT;
    SC__SEURAT__REPORT_TO_HTML;
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

workflow GENERATE_DUAL_INPUT_REPORT {
    
    take:
        pipelineStep
        data1
        data2
        rmd
    
    main:
        def reportTitle = "SC_Seurat_" + pipelineStep.toLowerCase() + "_report"

        SC__SEURAT__GENERATE_DUAL_INPUT_REPORT(
            rmd,
            data1,
            data2,
            reportTitle
        ) | SC__SEURAT__REPORT_TO_HTML
        
    emit:
        SC__SEURAT__GENERATE_DUAL_INPUT_REPORT.out
}

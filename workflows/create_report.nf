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

include '../processes/reports.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow GENERATE_DUAL_INPUT_REPORT {

    take:
        data  // Expects (sampleIdd, anndata1, anndata2)
        ipynb
        reportTitle
        isBenchmarkMode

    main:
        report_notebook = SC__SCANPY__GENERATE_DUAL_INPUT_REPORT(
            ipynb,
            data,
            reportTitle,
            isBenchmarkMode
        )
        SC__SCANPY__REPORT_TO_HTML(report_notebook)

    emit:
        report_notebook

}

workflow GENERATE_REPORT {

    take:
        pipelineStep
        data // anndata
        ipynb
        isBenchmarkMode

    main:
        def reportTitle = "SC_Scanpy_" + pipelineStep.toLowerCase() + "_report"
        if(isBenchmarkMode) {
            switch(pipelineStep) {
                case "CLUSTERING":
                    report_notebook = SC__SCANPY__BENCHMARK_CLUSTERING_GENERATE_REPORT(
                        ipynb,
                        // expects (sample_id, adata, ...arguments)
                        data,
                        reportTitle
                    )
                    break;
                default: 
                    throw new Exception("Invalid pipeline step")
                    break; 
            }
        } else {
            report_notebook = SC__SCANPY__GENERATE_REPORT(
                ipynb,
                // expects (sample_id, adata)
                data,
                reportTitle
            )
        }
        SC__SCANPY__REPORT_TO_HTML(report_notebook)

    emit:
        report_notebook

}

nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include {
    SC__SCANPY__GENERATE_DUAL_INPUT_REPORT;
    SC__SCANPY__PARAM_EXPLORE_CLUSTERING_GENERATE_REPORT;
    SC__SCANPY__GENERATE_REPORT;
    SC__SCANPY__REPORT_TO_HTML;
} from '../processes/reports.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow GENERATE_DUAL_INPUT_REPORT {

    take:
        data  // Expects (sampleIdd, anndata1, anndata2, stashedParams)
        ipynb
        reportTitle
        isParameterExplorationModeOn

    main:
        report_notebook = SC__SCANPY__GENERATE_DUAL_INPUT_REPORT(
            ipynb,
            data,
            reportTitle,
            isParameterExplorationModeOn
        )
        SC__SCANPY__REPORT_TO_HTML(report_notebook.map {
            it -> tuple(it[0], it[1])
        })

    emit:
        report_notebook

}

workflow GENERATE_REPORT {

    take:
        pipelineStep
        data // anndata
        ipynb
        isParameterExplorationModeOn

    main:
        def reportTitle = "SC_Scanpy_" + pipelineStep.toLowerCase() + "_report"
        if(isParameterExplorationModeOn) {
            switch(pipelineStep) {
                case "CLUSTERING":
                    report_notebook = SC__SCANPY__PARAM_EXPLORE_CLUSTERING_GENERATE_REPORT(
                        ipynb,
                        // expects (sample_id, adata, ...arguments)
                        data,
                        reportTitle
                    ).map {
                        it -> tuple(it[0], it[1], it[2..(it.size()-1)]) // stash params
                    }
                    break;
                default: 
                    throw new Exception("Invalid pipeline step")
                    break;
            }
        } else {
            switch(pipelineStep) {
                default:
                    report_notebook = SC__SCANPY__GENERATE_REPORT(
                        ipynb,
                        // expects (sample_id, adata)
                        data,
                        reportTitle
                    ).map {
                        it -> tuple(it[0], it[1], null)
                    }
                    break;
            }
        }
        SC__SCANPY__REPORT_TO_HTML(report_notebook.map {
            it -> tuple(it[0], it[1])
        })

    emit:
        report_notebook

}

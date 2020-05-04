nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:
// utils
include PUBLISH as PUBLISH_H5AD_HVG from "../../utils/workflows/utils.nf"

// scanpy:
include SC__SCANPY__FIND_HIGHLY_VARIABLE_GENES from '../processes/feature_selection.nf' params(params)
include SC__SCANPY__SUBSET_HIGHLY_VARIABLE_GENES from '../processes/feature_selection.nf' params(params)
include SC__SCANPY__FEATURE_SCALING from '../processes/transform.nf' params(params)

// reporting:
include GENERATE_REPORT from './create_report.nf' params(params)


//////////////////////////////////////////////////////

workflow HVG_SELECTION {

    take:
        data

    main:
        hvg = data \
            | SC__SCANPY__FIND_HIGHLY_VARIABLE_GENES \
            | SC__SCANPY__SUBSET_HIGHLY_VARIABLE_GENES
        PUBLISH_H5AD_HVG(
            SC__SCANPY__FIND_HIGHLY_VARIABLE_GENES.out
            null,
            "scanpy"
        )
        scaled = SC__SCANPY__FEATURE_SCALING( hvg )
        report = GENERATE_REPORT(
            "HVG",
            SC__SCANPY__FIND_HIGHLY_VARIABLE_GENES.out,
            file(workflow.projectDir + params.sc.scanpy.feature_selection.report_ipynb),
            false
        )

    emit:
        hvg
        scaled
        report

}

nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    PUBLISH as PUBLISH_SEURAT_RDS_NORMALIZED;
} from '../../utils/workflows/utils.nf' params(params)

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SEURAT__NORMALIZATION;
    SC__SEURAT__NORMALIZATION_SCT
    SC__SEURAT__SCALING;
} from '../processes/normalize_transform.nf' params(params)
include {
    GENERATE_REPORT;
} from './create_report.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow NORMALIZE {

    take:
        filtered
    
    main:
        SC__SEURAT__NORMALIZATION( filtered )
        PUBLISH_SEURAT_RDS_NORMALIZED(
            SC__SEURAT__NORMALIZATION.out,
            'SEURAT.normalized_output',
            'Rds',
            'seurat',
            false
        )
    emit:
        SC__SEURAT__NORMALIZATION.out
}

workflow NORMALIZE_SCALE_SCT {

    take:
        filtered
    
    main:
        scaled = SC__SEURAT__NORMALIZATION_SCT( filtered )
        PUBLISH_SEURAT_RDS_NORMALIZED(
            scaled,
            'SEURAT.normalized_sct_output',
            'Rds',
            'seurat',
            false
        )

        report = GENERATE_REPORT(
            "HVG",
            scaled,
            file(workflow.projectDir + params.tools.seurat.feature_selection.report_rmd)
        )

    emit:
        scaled
        report

}
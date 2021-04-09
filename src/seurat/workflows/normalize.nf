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
            'rds',
            'seurat',
            false
        )

        SC__SEURAT__SCALING( SC__SEURAT__NORMALIZATION.out )

    emit:
        SC__SEURAT__SCALING.out
}

workflow NORMALIZE_SCALE_SCT {

    take:
        filtered
    
    main:
        SC__SEURAT__NORMALIZATION_SCT( filtered )
        PUBLISH_SEURAT_RDS_NORMALIZED(
            SC__SEURAT__NORMALIZATION_SCT.out,
            'SEURAT.normalized_sct_output',
            'rds',
            'seurat',
            false
        )

    emit:
        SC__SEURAT__NORMALIZATION_SCT.out

}
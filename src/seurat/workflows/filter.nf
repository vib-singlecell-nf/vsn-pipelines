nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    PUBLISH as PUBLISH_SEURAT_RDS_FILTERED;
} from '../../utils/workflows/utils.nf' params(params)

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SEURAT__CELL_FILTER;
    SC__SEURAT__FEATURE_FILTER
} from '../processes/filter.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow FILTER {

    take:
        data
    
    main:
        filtered = data | \
            SC__SEURAT__CELL_FILTER | \
            SC__SEURAT__FEATURE_FILTER

        PUBLISH_SEURAT_RDS_FILTERED(
            filtered,
            'SEURAT.filtered_output',
            'Rds',
            'seurat',
            false
        )

    emit:
        filtered
}
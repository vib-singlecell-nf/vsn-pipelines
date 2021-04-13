nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    PUBLISH as PUBLISH_SEURAT_RDS_MERGED;
} from '../../utils/workflows/utils.nf' params(params)


////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SEURAT__MERGE;
} from '../processes/merge.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow MERGE {

    take:
        data
    
    main:
        SC__SEURAT__MERGE( data )

        PUBLISH_SEURAT_RDS_MERGED(
            SC__SEURAT__MERGE.out,
            'SEURAT.merged_output',
            'Rds',
            'seurat',
            false
        )

    emit:
        SC__SEURAT__MERGE.out
}

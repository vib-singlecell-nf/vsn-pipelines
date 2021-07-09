nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Process imports:
include {
    UPDATE_FEATURE_NOMENCLATURE
} from './updateFeatureNomenclature.nf' params(params)
include {
    FILTER_BY_CELL_METADATA
} from './filterByCellMetadata.nf' params(params)
include {
    STATIC__ANNOTATE_BY_CELL_METADATA
} from './annotateByCellMetadata.nf' params(params)
include {
    hasMetadataFilePath;
    SC__ANNOTATE_BY_SAMPLE_METADATA
} from '../processes/h5adAnnotate.nf' params(params)
include {
    SC__H5AD_BEAUTIFY;
} from '../processes/h5adUpdate.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow FILTER_AND_ANNOTATE_AND_CLEAN {

    take:
        // Expects (sampleId, h5ad) : Channel
        data

    main:
        out = data
        if(params.utils?.update_feature_metadata_index) {
            out = UPDATE_FEATURE_NOMENCLATURE( data )
        }
        // Filter cells based on an indexed cell-based metadata table
        if(params.utils?.cell_filter) {
            out = FILTER_BY_CELL_METADATA( out, 'NULL' )
        }
        // Annotate cells based on an indexed cell-based metadata table
        if(params.utils?.cell_annotate) {
            out = STATIC__ANNOTATE_BY_CELL_METADATA( 
                out,
                null
            )
        }
        // Annotate cells based on an indexed sample-based metadata table
        if(params.utils?.sample_annotate) {
            if (!hasMetadataFilePath(params.utils.sample_annotate)) {
                throw new Exception("The metadataFilePath param is missing in sample_annotate.")
            }
            out = SC__ANNOTATE_BY_SAMPLE_METADATA( out )
        }
        // Clean
        // e.g.: 
        // - h5ad: rename adata.obs values, remove adata.obs columns
        if(params.utils?.file_cleaner) {
            out = SC__H5AD_BEAUTIFY( out )
        }

    emit:
        out

}

/*
 * Conversion workflow 
 * Source:
 * 
 */ 

nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// Imports
include {
    SC__UTILS__EXTRACT_FEATURE_METADATA;
} from './../processes/h5adExtractMetadata' params(params)
include {
    FLYBASER__CONVERT_FBGN_TO_GENE_SYMBOL;
} from './../../flybaser/processes/convertNomenclature' params(params)
include {
    SC__UTILS__UPDATE_FEATURE_METADATA_INDEX;
} from './../processes/h5adUpdateMetadata' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow UPDATE_FEATURE_NOMENCLATURE {

    take:
        // Expects (sampleId, data)
        data

    main:
        SC__UTILS__EXTRACT_FEATURE_METADATA( data )
        FLYBASER__CONVERT_FBGN_TO_GENE_SYMBOL( SC__UTILS__EXTRACT_FEATURE_METADATA.out )
        out = SC__UTILS__UPDATE_FEATURE_METADATA_INDEX( data.join(FLYBASER__CONVERT_FBGN_TO_GENE_SYMBOL.out) )

    emit:
        out

}

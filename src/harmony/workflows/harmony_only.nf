nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include {
    SC__HARMONY__HARMONY_MATRIX;
} from './../processes/runHarmony.nf' params(params)
include {
SC__H5AD_UPDATE_X_PCA;
} from './../../utils/processes/h5adUpdate.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow BEC_HARMONY {

    take:
        data

    main:
        // Run Harmony
        harmony_embeddings = SC__HARMONY__HARMONY_MATRIX( 
            data.map { 
                it -> tuple(it[0], it[1])
            } 
        )
        SC__H5AD_UPDATE_X_PCA( 
            data.map {
                it -> tuple(it[0], it[1]) 
            }.join(harmony_embeddings) 
        )

    emit:
        SC__H5AD_UPDATE_X_PCA.out
}

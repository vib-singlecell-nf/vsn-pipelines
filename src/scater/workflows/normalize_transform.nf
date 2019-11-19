nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include SC__SCATER__DATA_TRANSFORMATION from '../processes/transform.nf' params(params.sc.scater.data_transformation + params.global + params)
include SC__SCATER__NORMALIZATION from '../processes/transform.nf' params(params.sc.scater.normalization + params.global + params)

//////////////////////////////////////////////////////

workflow NORMALIZE_TRANSFORM {
    get:
        filtered
    main:
        SC__SCATER__NORMALIZATION( filtered )
        SC__SCATER__DATA_TRANSFORMATION( SC__SCATER__NORMALIZATION.out )
    emit:
        SC__SCATER__DATA_TRANSFORMATION.out
}

nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include SC__SCANPY__DATA_TRANSFORMATION from '../processes/transform.nf' params(params)
include SC__SCANPY__NORMALIZATION from '../processes/transform.nf' params(params)

//////////////////////////////////////////////////////

workflow NORMALIZE_TRANSFORM {
    get:
        filtered
    main:
        SC__SCANPY__NORMALIZATION( filtered )
        SC__SCANPY__DATA_TRANSFORMATION( SC__SCANPY__NORMALIZATION.out )
    emit:
        SC__SCANPY__DATA_TRANSFORMATION.out
}

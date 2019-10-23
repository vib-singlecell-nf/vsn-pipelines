nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include SC__SCANPY__FEATURE_SELECTION from '../processes/feature_selection.nf' params(params.sc.scanpy.feature_selection + params.global + params)
include SC__SCANPY__FEATURE_SCALING from '../processes/transform.nf' params(params.sc.scanpy.feature_scaling + params.global + params)

//////////////////////////////////////////////////////

workflow HVG_SELECTION {
    get:
        data
    main:
        SC__SCANPY__FEATURE_SELECTION( data )
        SC__SCANPY__FEATURE_SCALING( SC__SCANPY__FEATURE_SELECTION.out )
    emit:
        SC__SCANPY__FEATURE_SCALING.out
}

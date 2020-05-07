nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:
// utils
include PUBLISH as PUBLISH_H5AD_NORMALIZED from "../../utils/workflows/utils.nf" params(params)


// scanpy:
include SC__SCANPY__DATA_TRANSFORMATION from '../processes/transform.nf' params(params)
include SC__SCANPY__NORMALIZATION from '../processes/transform.nf' params(params)

//////////////////////////////////////////////////////

workflow NORMALIZE_TRANSFORM {

    take:
        filtered

    main:
        SC__SCANPY__NORMALIZATION( filtered )
        PUBLISH_H5AD_NORMALIZED(
            SC__SCANPY__NORMALIZATION.out,
            "SCANPY.normalized_output",
            "scanpy",
            false
        )
        SC__SCANPY__DATA_TRANSFORMATION( SC__SCANPY__NORMALIZATION.out )

    emit:
        SC__SCANPY__DATA_TRANSFORMATION.out

}

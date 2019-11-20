nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include QC_FILTER from '../src/scanpy/workflows/qc_filter.nf' params(params)
include SC__FILE_CONCATENATOR from '../src/utils/processes/utils.nf' params(params.sc.file_concatenator + params.global + params)
include NORMALIZE_TRANSFORM from '../src/scanpy/workflows/normalize_transform.nf' params(params + params.global)
include HVG_SELECTION from '../src/scanpy/workflows/hvg_selection.nf' params(params + params.global)
include SC__SCANPY__ADJUSTMENT from '../src/scanpy/processes/adjust.nf' params(params + params.global)
include SC__H5AD_TO_FILTERED_LOOM from '../src/utils/processes/h5adToLoom.nf' params(params + params.global)
include BEC_MNN_CORRECT from '../src/scanpy/workflows/bec_mnn_correct.nf' params(params.sc.file_concatenator + params)

// data channel to start from 10x data:
include getChannel as getTenXChannel from '../src/channels/tenx.nf' params(params)

workflow mnncorrect {

    data = getTenXChannel( params.global.tenx_folder ).view()
    QC_FILTER( data ) // Remove concat
    SC__FILE_CONCATENATOR( QC_FILTER.out.filtered.map{it -> it[1]}.collect() )
    NORMALIZE_TRANSFORM( SC__FILE_CONCATENATOR.out )
    HVG_SELECTION( NORMALIZE_TRANSFORM.out )
    // SC__SCANPY__ADJUSTMENT( HVG_SELECTION.out.scaled )
    filteredloom = SC__H5AD_TO_FILTERED_LOOM( SC__FILE_CONCATENATOR.out )
    scopeloom = BEC_MNN_CORRECT( HVG_SELECTION.out.scaled )

    emit:
        filteredloom
        scopeloom
}


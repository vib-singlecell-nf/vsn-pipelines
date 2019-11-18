//
// Version:
// Test:
// Command: 
// 
//
/*
 * Remote run test
 * Source:
 * 
 * Steps considered: 

 */ 
import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2

// print all parameters:
// println(prettyPrint(toJson( params )))


//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include QC_FILTER from '../src/scanpy/workflows/qc_filter.nf' params(params)
include SC__FILE_CONCATENATOR from '../src/utils/processes/utils.nf' params(params.sc.file_concatenator + params.global + params)
include NORMALIZE_TRANSFORM from '../src/scanpy/workflows/normalize_transform.nf' params(params + params.global)
include HVG_SELECTION from '../src/scanpy/workflows/hvg_selection.nf' params(params + params.global)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__PCA from '../src/scanpy/processes/dim_reduction.nf' params(params.sc.scanpy.dim_reduction.pca + params.global + params)
include SC__H5AD_TO_FILTERED_LOOM from '../src/utils/processes/h5ad_to_loom.nf' params(params + params.global)
include BEC_BBKNN from '../src/scanpy/workflows/bec_bbknn.nf' params(params)

// data channel to start from 10x data:
include getChannel as getTenXChannel from '../src/channels/tenx.nf' params(params)

workflow bbknn {

    data = getTenXChannel( params.global.tenx_folder ).view()
    QC_FILTER( data ) // Remove concat
    SC__FILE_CONCATENATOR( QC_FILTER.out.filtered.map{it -> it[1]}.collect() )
    NORMALIZE_TRANSFORM( SC__FILE_CONCATENATOR.out )
    HVG_SELECTION( NORMALIZE_TRANSFORM.out )
    SC__SCANPY__DIM_REDUCTION__PCA( HVG_SELECTION.out.scaled )
    filteredloom = SC__H5AD_TO_FILTERED_LOOM( SC__FILE_CONCATENATOR.out )
    scopeloom = BEC_BBKNN( SC__SCANPY__DIM_REDUCTION__PCA.out )

    emit:
        filteredloom
        scopeloom
}


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

include CELLRANGER from './src/cellranger/main.nf' params(params)
include QC_FILTER from './src/scanpy/qc_filter.nf' params(params)
include BEC_BBKNN from './src/scanpy/bec_bbknn.nf' params(params)

include SC__H5AD_TO_FILTERED_LOOM from './src/utils/processes/h5ad_to_loom.nf' params(params + params.global)
include SCENIC_append from './src/scenic/main.nf' params(params)

// data channel to start from 10x data:
include getChannel as getTenXChannel from './src/channels/tenx.nf' params(params)


workflow {
    // CELLRANGER()
    // QC_FILTER( CELLRANGER.out )
    // BEC_BBKNN( QC_FILTER.out )

    // to run starting from the 10x output:
    data = getTenXChannel( params.global.tenx_folder ).view()
    QC_FILTER( data )
    scopeloom = BEC_BBKNN( QC_FILTER.out )
    filteredloom = SC__H5AD_TO_FILTERED_LOOM( QC_FILTER.out )
    SCENIC_append( filteredloom, scopeloom )

}



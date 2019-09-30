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
println(prettyPrint(toJson( params )))


//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include CELLRANGER from './src/cellranger/main.nf' params(params)
include BEC_BBKNN from './src/scanpy/bec_bbknn.nf' params(params)
include SCENIC from './src/scenic/main.nf' params(params)

// data channel to start from 10x data:
include getChannel as getTenXChannel from './src/channels/tenx.nf' params(params)


workflow {
    CELLRANGER()
    BEC_BBKNN( CELLRANGER.out )

    // to run BEC_BBKNN starting from the 10x output:
    // data = getTenXChannel( params.global.tenx_folder ).view()
    // BEC_BBKNN( data )

    // SCENIC( file( params.sc.scenic.filteredloom ) )
}



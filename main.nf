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


// params.each { println "${it}" }
println "${workflow.projectDir}"
println(prettyPrint(toJson( params )))


//////////////////////////////////////////////////////
//  Define the parameters for current testing proces

// include SCENIC from './src/scenic/main.nf' params(params)

include CELLRANGER from './src/cellranger/main.nf' params(params)

include getChannel as getTenXChannel from './src/channels/tenx.nf' params(params)
include BEC_BBKNN from './src/scanpy/bec_bbknn.nf' params(params)



workflow {
    // SCENIC( file( params.sc.scenic.filteredloom ) )
    // CELLRANGER()
    // BEC_BBKNN( CELLRANGER.out )

    data = getTenXChannel( params.global.tenx_folder ).view()
    BEC_BBKNN( data )
}











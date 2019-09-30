//
// Tested on 23/08/2019
// Command: 
// nextflow src/processes/utils/utils_1.test.nf --tenx_folder data/01.count/**/filtered_feature_bc_matrix -resume
//

nextflow.preview.dsl=2

///////////////////////////////////////////
//  Define the parameters for all processes

include 'utils' params(iff: '10x_mtx', off: 'h5ad')

// Make the test workflow 
workflow test {
    get:
        data
    main:
        SC__FILE_CONVERTER( data )
    emit:
        SC__FILE_CONVERTER.out
}

// Uncomment to test
// include getChannel as getTenXChannel from '../channels/tenx.nf'
// test( getTenXChannel( params.tenx_folder ) )
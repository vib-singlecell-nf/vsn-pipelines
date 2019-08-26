//
// Tested on 23/08/2019
// Command: 
// nextflow src/processes/utils/utils_2.test.nf --tenx_folder data/01.count/**/filtered_feature_bc_matrix --project_name tiny -resume
//

nextflow.preview.dsl=2

///////////////////////////////////////////
//  Define the parameters for all processes

include 'utils' params(iff: '10x_mtx', off: 'h5ad') params(project_name: params.project_name)

// Make the test workflow 
workflow test( f ) {
    main:
        SC__FILE_CONCATENATOR( ( f ).collect() )
    emit:
        SC__FILE_CONCATENATOR.out
}

// Uncomment to test
// include getChannel as getTenXChannel from '../../channels/tenx.nf'
// test( getTenXChannel( params.tenx_folder ) )
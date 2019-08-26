//
// Tested on 23/08/2019
// Command: 
// nextflow src/processes/scanpy/transform_2.test.nf --tenx_folder data/01.count/**/filtered_feature_bc_matrix -resume
//
nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Define the parameters for current testing process
include 'transform' params(
    featureScalingMthod: 'zscore_scale',
    featureScalingMaxSD: -1,
    iff: '10x_mtx',
    off: 'h5ad'
)

// Make the test workflow 
workflow test_scripts_scanpy_transform_2( f ) {
    main:
        SC__SCANPY__FEATURE_SCALING( test_scripts_scanpy__adjust( f ) )
    emit:
        SC__SCANPY__FEATURE_SCALING.out
}

// Uncomment to test
include 'adjust.test.nf' params(tenx_folder: params.tenx_folder)
// include getChannel as getTenXChannel from '../../channels/tenx.nf'
// test_scripts_scanpy_transform_2( getTenXChannel( params.tenx_folder ) ).collect().view()
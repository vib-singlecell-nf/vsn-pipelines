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
include QC_FILTER from './src/scanpy/workflows/qc_filter.nf' params(params)
include SC__FILE_CONCATENATOR from './src/utils/processes/utils.nf' params(params.sc.file_concatenator + params.global + params)
include NORMALIZE_TRANSFORM from './src/scanpy/workflows/normalize_transform.nf' params(params)
include HVG_SELECTION from './src/scanpy/workflows/hvg_selection.nf' params(params)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__PCA from './src/scanpy/processes/dim_reduction.nf' params(params.sc.scanpy.dim_reduction.pca + params.global + params)
include DIM_REDUCTION from './src/scanpy/workflows/dim_reduction.nf' params(params)
include CLUSTER_IDENTIFICATION from './src/scanpy/workflows/cluster_identification.nf' params(params)

include BEC_BBKNN from './src/scanpy/workflows/bec_bbknn.nf' params(params)

include SC__H5AD_TO_FILTERED_LOOM from './src/utils/processes/h5ad_to_loom.nf' params(params + params.global)
include SCENIC_append from './src/scenic/main.nf' params(params)

// data channel to start from 10x data:
include getChannel as getTenXChannel from './src/channels/tenx.nf' params(params)


// TODO: Split into workflows to run with -entry

workflow bbknn_scenic {
    // CELLRANGER()

    data = getTenXChannel( params.global.tenx_folder ).view()
    QC_FILTER( data ) // Remove concat
    SC__FILE_CONCATENATOR( QC_FILTER.out.collect() )
    NORMALIZE_TRANSFORM( SC__FILE_CONCATENATOR.out )
    HVG_SELECTION( NORMALIZE_TRANSFORM.out )
    SC__SCANPY__DIM_REDUCTION__PCA( HVG_SELECTION.out )
    CLUSTER_IDENTIFICATION( SC__SCANPY__DIM_REDUCTION__PCA.out )
    scopeloom = BEC_BBKNN( CLUSTER_IDENTIFICATION.out )
    
    filteredloom = SC__H5AD_TO_FILTERED_LOOM( QC_FILTER.out )
    SCENIC_append( filteredloom, scopeloom )

}

workflow single_sample {
    
    data = getTenXChannel( params.global.tenx_folder ).view()
    QC_FILTER( data ) // Remove concat
    NORMALIZE_TRANSFORM( QC_FILTER.out )
    HVG_SELECTION( NORMALIZE_TRANSFORM.out )
    DIM_REDUCTION( HVG_SELECTION.out )
    CLUSTER_IDENTIFICATION( DIM_REDUCTION.out )
    filteredloom = SC__H5AD_TO_FILTERED_LOOM( CLUSTER_IDENTIFICATION.out )
    
}

//
// Version: 0.1.0
// Test: passed
// Command: 
//  nextflow run src/singlecelltxbenchmark/pipelines/bec__bbknn -profile singularity --tenx_folder data/01.count/**/filtered_feature_bc_matrix --project_name tiny
//
/*
 * BEC__BBKNN workflow 
 * Source: https://github.com/Teichlab/bbknn/blob/master/examples/pancreas.ipynb
 * 
 * Steps considered: 
 * - normalize
 * - concatenate the batches
 * - feature selection
 * - log transform
 * - feature scaling
 * - dimensionality reduction (PCA)
 * - batch effect correction using python package bbknn (Park et al. (2018), Fast Batch Alignment of Single Cell Transcriptomes Unifies Multiple Mouse Cell Atlases into an Integrated Landscape)
 */ 

nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:


include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__PCA from '../processes/dim_reduction.nf' params(params.sc.scanpy.dim_reduction.pca + params.global + params)

include '../processes/batch_effect_correct.nf' params(params.sc.scanpy.batch_effect_correct + params.global + params)

include SC__SCANPY__CLUSTERING from '../processes/cluster.nf' params(params.sc.scanpy.clustering + params.global + params)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__UMAP from '../processes/dim_reduction.nf' params(params.sc.scanpy.dim_reduction.umap + params.global + params)
include SC__H5AD_TO_LOOM from '../../utils/processes/h5ad_to_loom.nf' params(params.global + params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow BEC_BBKNN {
    get:
        data
    main:
        SC__SCANPY__DIM_REDUCTION__PCA( data )
        SC__SCANPY__BATCH_EFFECT_CORRECTION( SC__SCANPY__DIM_REDUCTION__PCA.out )
        SC__SCANPY__CLUSTERING( SC__SCANPY__BATCH_EFFECT_CORRECTION.out )
        SC__SCANPY__DIM_REDUCTION__UMAP( SC__SCANPY__CLUSTERING.out )
        SC__H5AD_TO_LOOM(SC__SCANPY__DIM_REDUCTION__UMAP.out )
        // Not using t-SNE as it does not use the neighbour graph (which BBKNN alters) when constructing its dimensionality reduction
    emit:
        SC__H5AD_TO_LOOM.out
}

// Uncomment to test
// workflow {
//     main:
//         bbknn( getTenXChannel( params.tenx_folder ) )
// }

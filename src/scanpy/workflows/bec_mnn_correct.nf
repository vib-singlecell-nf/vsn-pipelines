/*
 * BEC_MNN_CORRECT workflow 
 * - batch effect correction using python package mnnpy (fast and python version of mnnCorrect (Haghverdi et al, 2018)
 */ 

nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include '../processes/batch_effect_correct.nf' params(params.sc.scanpy.batch_effect_correct + params.global + params)
include CLUSTER_IDENTIFICATION from './cluster_identification.nf' params(params + params.global)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__PCA from '../processes/dim_reduction.nf' params(params.sc.scanpy.dim_reduction.pca + params.global + params)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__UMAP from '../processes/dim_reduction.nf' params(params.sc.scanpy.dim_reduction.umap + params.global + params)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__TSNE from '../processes/dim_reduction.nf' params(params.sc.scanpy.dim_reduction.tsne + params.global + params)
include SC__H5AD_TO_LOOM from '../../utils/processes/h5adToLoom.nf' params(params.global + params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow BEC_MNN_CORRECT {
    get:
        data
    main:
        SC__SCANPY__BATCH_EFFECT_CORRECTION( data )
        SC__SCANPY__DIM_REDUCTION__PCA( SC__SCANPY__BATCH_EFFECT_CORRECTION.out )
        SC__SCANPY__DIM_REDUCTION__UMAP( SC__SCANPY__DIM_REDUCTION__PCA.out )
        SC__SCANPY__DIM_REDUCTION__TSNE( SC__SCANPY__DIM_REDUCTION__UMAP.out )
        CLUSTER_IDENTIFICATION( SC__SCANPY__DIM_REDUCTION__TSNE.out )
        SC__H5AD_TO_LOOM( CLUSTER_IDENTIFICATION.out.marker_genes )
    emit:
        SC__H5AD_TO_LOOM.out
}

// Uncomment to test
// workflow {
//     main:
//         bbknn( getTenXChannel( params.tenx_folder ) )
// }

/*
 * BEC_MNN_CORRECT workflow 
 * - batch effect correction using python package mnnpy (fast and python version of mnnCorrect (Haghverdi et al, 2018)
 */ 

nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include '../processes/batch_effect_correct.nf' params(params)
include DIM_REDUCTION_PCA from './dim_reduction_pca' params(params + [method: "pca"])
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__UMAP from '../processes/dim_reduction.nf' params(params + [method: "umap"])
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__TSNE from '../processes/dim_reduction.nf' params(params + [method: "tsne"])
include CLUSTER_IDENTIFICATION from './cluster_identification.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow BEC_MNN_CORRECT {

    get:
        normalizedTransformedData
        data

    main:
        SC__SCANPY__BATCH_EFFECT_CORRECTION( data )
        DIM_REDUCTION_PCA( SC__SCANPY__BATCH_EFFECT_CORRECTION.out )
        SC__SCANPY__DIM_REDUCTION__UMAP( DIM_REDUCTION_PCA.out )
        SC__SCANPY__DIM_REDUCTION__TSNE( SC__SCANPY__DIM_REDUCTION__UMAP.out )
        CLUSTER_IDENTIFICATION(
            normalizedTransformedData,
            SC__SCANPY__DIM_REDUCTION__TSNE.out
        )

    emit:
        data = CLUSTER_IDENTIFICATION.out.marker_genes

}

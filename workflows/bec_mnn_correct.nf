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
include NEIGHBORHOOD_GRAPH from '../src/scanpy/workflows/neighborhood_graph.nf' params(params)
include DIM_REDUCTION_TSNE_UMAP from './dim_reduction' params(params)
include CLUSTER_IDENTIFICATION from './cluster_identification.nf' params(params)


//////////////////////////////////////////////////////
//  Define the workflow 

workflow BEC_MNN_CORRECT {

    take:
        normalizedTransformedData
        data

    main:
        SC__SCANPY__BATCH_EFFECT_CORRECTION( data )
        DIM_REDUCTION_PCA( SC__SCANPY__BATCH_EFFECT_CORRECTION.out )
        NEIGHBORHOOD_GRAPH( DIM_REDUCTION_PCA.out )
        DIM_REDUCTION_TSNE_UMAP( DIM_REDUCTION_PCA.out )
        CLUSTER_IDENTIFICATION(
            normalizedTransformedData,
            DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
            "Post Batch Effect Correction (MNN CORRECT)"
        )

    emit:
        data = CLUSTER_IDENTIFICATION.out.marker_genes

}

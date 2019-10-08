nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__PCA from '../processes/dim_reduction.nf' params(params.sc.scanpy.dim_reduction.pca + params.global + params)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__TSNE from '../processes/dim_reduction.nf' params(params.sc.scanpy.dim_reduction.tsne + params.global + params)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__UMAP from '../processes/dim_reduction.nf' params(params.sc.scanpy.dim_reduction.umap + params.global + params)

//////////////////////////////////////////////////////

workflow DIM_REDUCTION {
    get:
        data
    main:
        SC__SCANPY__DIM_REDUCTION__PCA( data )
        SC__SCANPY__DIM_REDUCTION__TSNE( SC__SCANPY__DIM_REDUCTION__PCA.out )
        SC__SCANPY__DIM_REDUCTION__UMAP( SC__SCANPY__DIM_REDUCTION__TSNE.out )
    emit:
        SC__SCANPY__DIM_REDUCTION__UMAP.out
}

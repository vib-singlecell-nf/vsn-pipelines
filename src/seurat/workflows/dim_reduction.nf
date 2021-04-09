nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SEURAT__DIM_REDUCTION as SC__SEURAT__DIM_REDUCTION__TSNE;
} from '../processes/dim_reduction.nf' params(params + [method: "tsne"])
include {
    SC__SEURAT__DIM_REDUCTION as SC__SEURAT__DIM_REDUCTION__UMAP;
 } from '../processes/dim_reduction.nf' params(params + [method: "umap"])

workflow DIM_REDUCTION_TSNE_UMAP {
    take:
        data

    main:
        SC__SEURAT__DIM_REDUCTION__TSNE( data )
        SC__SEURAT__DIM_REDUCTION__UMAP( SC__SEURAT__DIM_REDUCTION__TSNE.out )

    emit:
        dimred_tsne = SC__SEURAT__DIM_REDUCTION__TSNE.out
        dimred_tsne_umap = SC__SEURAT__DIM_REDUCTION__UMAP.out
}
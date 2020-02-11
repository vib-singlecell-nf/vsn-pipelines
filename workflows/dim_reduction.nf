nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include DIM_REDUCTION_PCA from './dim_reduction_pca' params(params)
// include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__PCA from '../processes/dim_reduction.nf' params(params + [method: "pca"])
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__TSNE from '../processes/dim_reduction.nf' params(params + [method: "tsne"])
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__UMAP from '../processes/dim_reduction.nf' params(params + [method: "umap"])

// reporting:
include GENERATE_REPORT from './create_report.nf' params(params)

//////////////////////////////////////////////////////

workflow DIM_REDUCTION {

    take:
        data

    main:
        dimred_pca = DIM_REDUCTION_PCA( data )
        DIM_REDUCTION_TSNE_UMAP( dimred_pca )

        report = GENERATE_REPORT(
            "DIMENSIONALITY_REDUCTION",
            DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap.map { it -> tuple(it[0], it[1]) },
            file(workflow.projectDir + params.sc.scanpy.dim_reduction.report_ipynb),
            false
        )

    emit:
        dimred_pca_tsne = DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne
        dimred_pca_tsne_umap = DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap
        report

}

workflow DIM_REDUCTION_TSNE_UMAP {

    take:
        data

    main:
        dimred_tsne = SC__SCANPY__DIM_REDUCTION__TSNE( data.map {
                item -> tuple(item[0], item[1], null, null)
        })
        dimred_tsne_umap = SC__SCANPY__DIM_REDUCTION__UMAP( SC__SCANPY__DIM_REDUCTION__TSNE.out.map {
                item -> tuple(item[0], item[1], null, null)
        })

    emit:
        dimred_tsne
        dimred_tsne_umap

}

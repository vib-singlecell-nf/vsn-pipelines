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

    get:
        data

    main:
        DIM_REDUCTION_PCA( data )
        // SC__SCANPY__DIM_REDUCTION__PCA( data )
        // SC__SCANPY__DIM_REDUCTION__TSNE( SC__SCANPY__DIM_REDUCTION__PCA.out )
        SC__SCANPY__DIM_REDUCTION__TSNE( DIM_REDUCTION_PCA.out )
        dimred = SC__SCANPY__DIM_REDUCTION__UMAP( SC__SCANPY__DIM_REDUCTION__TSNE.out )
        report = GENERATE_REPORT(
            SC__SCANPY__DIM_REDUCTION__UMAP.out,
            file(workflow.projectDir + params.sc.scanpy.dim_reduction.report_ipynb),
            "SC_dimensionality_reduction_report"
        )

    emit:
        dimred
        report

}

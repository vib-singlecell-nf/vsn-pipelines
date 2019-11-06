nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// scanpy:
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__PCA from '../processes/dim_reduction.nf' params(params.sc.scanpy.dim_reduction.pca + params.global + params)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__TSNE from '../processes/dim_reduction.nf' params(params.sc.scanpy.dim_reduction.tsne + params.global + params)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__UMAP from '../processes/dim_reduction.nf' params(params.sc.scanpy.dim_reduction.umap + params.global + params)

// reporting:
include GENERATE_REPORT from './create_report.nf' params(params.sc.scanpy.feature_scaling + params)

//////////////////////////////////////////////////////

workflow DIM_REDUCTION {
    get:
        data
    main:
        SC__SCANPY__DIM_REDUCTION__PCA( data )
        SC__SCANPY__DIM_REDUCTION__TSNE( SC__SCANPY__DIM_REDUCTION__PCA.out )
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


nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SEURAT__DIM_REDUCTION as SC__SEURAT__DIM_REDUCTION__TSNE;
} from '../processes/dim_reduction.nf' params(params + [method: "tsne"])
include {
    SC__SEURAT__DIM_REDUCTION as SC__SEURAT__DIM_REDUCTION__UMAP;
} from '../processes/dim_reduction.nf' params(params + [method: "umap"])
include {
    GENERATE_REPORT;
} from './create_report.nf' params(params)

workflow DIM_REDUCTION_TSNE_UMAP {
    take:
        data

    main:
        SC__SEURAT__DIM_REDUCTION__TSNE( data )
        SC__SEURAT__DIM_REDUCTION__UMAP( SC__SEURAT__DIM_REDUCTION__TSNE.out )

        // Drop nPcs from the output, since we don't need it anymore in any of the following steps
        dimred_tsne = SC__SEURAT__DIM_REDUCTION__TSNE.out.map { it -> tuple(it[0], it[1]) }
        dimred_tsne_umap = SC__SEURAT__DIM_REDUCTION__UMAP.out.map { it -> tuple(it[0], it[1]) }

        report = GENERATE_REPORT(
            "DIMENSIONALITY_REDUCTION",
            dimred_tsne_umap,
            file(workflow.projectDir + params.tools.seurat.dim_reduction.report_rmd)
        )

    emit:
        dimred_tsne
        dimred_tsne_umap
        report
}
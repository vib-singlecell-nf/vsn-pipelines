/*
 * BEC_MNNCORRECT workflow 
 * - batch effect correction using python package mnnpy (fast and python version of mnnCorrect (Haghverdi et al, 2018)
 */ 

nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include '../../utils/processes/utils.nf' params(params)
include '../../utils/workflows/utils.nf' params(params)

// scanpy:
include '../processes/batch_effect_correct.nf' params(params)
include DIM_REDUCTION_PCA from './dim_reduction_pca' params(params + [method: "pca"])
include NEIGHBORHOOD_GRAPH from './neighborhood_graph.nf' params(params)
include DIM_REDUCTION_TSNE_UMAP from './dim_reduction' params(params)
include '../processes/dim_reduction.nf' params(params)
include '../processes/cluster.nf' params(params)
include './cluster_identification.nf' params(params) // Don't only import a specific process (the function needs also to be imported)

// reporting:
include GENERATE_DUAL_INPUT_REPORT from './create_report.nf' params(params + params.global)


//////////////////////////////////////////////////////
//  Define the workflow 

workflow BEC_MNNCORRECT {

    take:
        normalizedTransformedData
        data
        // Expects (sampleId, anndata)
        clusterIdentificationPreBatchEffectCorrection

    main:
        SC__SCANPY__BATCH_EFFECT_CORRECTION( data.map { it -> tuple(it[0], it[1], null) } )
        DIM_REDUCTION_PCA( SC__SCANPY__BATCH_EFFECT_CORRECTION.out )
        NEIGHBORHOOD_GRAPH( DIM_REDUCTION_PCA.out )

        // Run dimensionality reduction
        DIM_REDUCTION_TSNE_UMAP( NEIGHBORHOOD_GRAPH.out )

        // Define the parameters for clustering
        def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.sc.scanpy.clustering) )
        CLUSTER_IDENTIFICATION(
            normalizedTransformedData,
            DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
            "Post Batch Effect Correction (MNNCORRECT)"
        )

        marker_genes = CLUSTER_IDENTIFICATION.out.marker_genes.map {
            it -> tuple(
                it[0], // sampleId
                it[1], // data
                !clusteringParams.isParameterExplorationModeOn() ? null : it[2..(it.size()-1)], // Stash params
            )
        }

        SC__PUBLISH_H5AD( 
            marker_genes.map {
                it -> tuple(it[0], it[1], it[2])
            },
            "BEC_MNNCORRECT.output"
        )

        // This will generate a dual report with results from
        // - Pre batch effect correction
        // - Post batch effect correction
        becDualDataPrePost = COMBINE_BY_PARAMS(
            clusterIdentificationPreBatchEffectCorrection,
            SC__PUBLISH_H5AD.out,
            clusteringParams
        )

        mnncorrect_report = GENERATE_DUAL_INPUT_REPORT(
            becDualDataPrePost.map { it -> tuple(it[0], it[1], it[2]) },
            file(workflow.projectDir + params.sc.scanpy.batch_effect_correct.report_ipynb),
            "SC_BEC_MNNCORRECT_report",
            clusteringParams.isParameterExplorationModeOn()
        )

    emit:
        data = DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap
        cluster_report = CLUSTER_IDENTIFICATION.out.report
        mnncorrect_report

}

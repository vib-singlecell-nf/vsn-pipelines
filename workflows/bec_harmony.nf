nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include '../../utils/processes/utils.nf' params(params)
include '../../utils/workflows/utils.nf' params(params)

include SC__HARMONY__HARMONY_MATRIX from './../processes/runHarmony.nf' params(params)
include SC__H5AD_UPDATE_X_PCA from './../../utils/processes/h5adUpdate.nf' params(params)
include DIM_REDUCTION_TSNE_UMAP from './../../scanpy/workflows/dim_reduction.nf' params(params)
include './../../scanpy/processes/cluster.nf' params(params)
include './../../scanpy/workflows/cluster_identification.nf' params(params) // Don't only import a specific process (the function needs also to be imported)

// reporting:
include GENERATE_DUAL_INPUT_REPORT from './../../scanpy/workflows/create_report.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow BEC_HARMONY {

    take:
        normalizedTransformedData
        dimReductionData
        // Expects (sampleId, anndata)
        clusterIdentificationPreBatchEffectCorrection

    main:
        // Run Harmony
        harmony_embeddings = SC__HARMONY__HARMONY_MATRIX( dimReductionData )
        SC__H5AD_UPDATE_X_PCA( 
            dimReductionData.join(harmony_embeddings) 
        )

        // Run dimensionality reduction
        DIM_REDUCTION_TSNE_UMAP( SC__H5AD_UPDATE_X_PCA.out )

        // Run clustering
        // Define the parameters for clustering
        def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.sc.scanpy.clustering) )
        CLUSTER_IDENTIFICATION(
            normalizedTransformedData,
            DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
            "Post Batch Effect Correction (Harmony)"
        )

        SC__PUBLISH_H5AD( 
            CLUSTER_IDENTIFICATION.out.marker_genes.map { it -> tuple(it[0], it[1]) },
            "BEC_HARMONY.output"
        )
        
        // This will generate a dual report with results from
        // - Pre batch effect correction
        // - Post batch effect correction
        becDualDataPrePost = COMBINE_BY_PARAMS(
            clusterIdentificationPreBatchEffectCorrection,
            // Use SC__PUBLISH_H5AD output to avoid "input file name collision"
            SC__PUBLISH_H5AD.out,
            clusteringParams
        )
        harmony_report = GENERATE_DUAL_INPUT_REPORT(
            becDualDataPrePost,
            file(workflow.projectDir + params.sc.harmony.report_ipynb),
            "SC_BEC_HARMONY_report",
            clusteringParams.isBenchmarkMode()
        )

    emit:
        data = CLUSTER_IDENTIFICATION.out.marker_genes
        cluster_report = CLUSTER_IDENTIFICATION.out.report
        harmony_report

}

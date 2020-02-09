/*
 * BEC__BBKNN workflow 
 * Source: https://github.com/Teichlab/bbknn/blob/master/examples/pancreas.ipynb
 * 
 * Steps considered: 
 * - normalize
 * - concatenate the batches
 * - feature selection
 * - log transform
 * - feature scaling
 * - dimensionality reduction (PCA)
 * - batch effect correction using python package bbknn (Park et al. (2018), Fast Batch Alignment of Single Cell Transcriptomes Unifies Multiple Mouse Cell Atlases into an Integrated Landscape)
 */ 

nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include '../../utils/processes/utils.nf' params(params)

// scanpy:
include '../processes/batch_effect_correct.nf' params(params)

include '../processes/cluster.nf' params(params)
include '../processes/dim_reduction.nf' params(params)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__UMAP from '../processes/dim_reduction.nf' params(params + [method: "umap"])
include './cluster_identification.nf' params(params) // Don't only import a specific process (the function needs also to be imported)

// reporting:
include GENERATE_DUAL_INPUT_REPORT from './create_report.nf' params(params + params.global)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow BEC_BBKNN {

    take:
        normalizedTransformedData
        dimReductionData
        // Expects (sampleId, anndata)
        clusterIdentificationPreBatchEffectCorrection

    main:
        SC__SCANPY__BATCH_EFFECT_CORRECTION( 
            dimReductionData.map { 
                it -> tuple(it[0], it[1], it[2]) 
            } 
        )

        // Define the parameters for clustering
        def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.sc.scanpy.clustering) )
        CLUSTER_IDENTIFICATION(
            normalizedTransformedData,
            SC__SCANPY__BATCH_EFFECT_CORRECTION.out,
            "Post Batch Effect Correction (BBKNN)"
        )

        // Define the parameters for dimensionality reduction
        def dimRedParams = SC__SCANPY__DIM_REDUCTION_PARAMS( clean(params.sc.scanpy.dim_reduction.umap) )

        SC__SCANPY__DIM_REDUCTION__UMAP( 
            CLUSTER_IDENTIFICATION.out.marker_genes.map {
                it -> tuple(
                    it[0], // sampleId
                    it[1], // data
                    !clusteringParams.isBenchmarkMode() ? null : it[2..(it.size()-1)], // set runtime process parameters as inert
                )
            }.combine(
                dimRedParams.$()
            )
        )

        SC__PUBLISH_H5AD( 
            SC__SCANPY__DIM_REDUCTION__UMAP.out.map { it -> tuple(it[0], it[1]) },
            "BEC_BBKNN.output"
        )

        // This will generate a dual report with results from
        // - CLUSTER_IDENTIFICATION pre batch effect correction
        // - CLUSTER_IDENTIFICATION post batch effect correction
        if(clusteringParams.isBenchmarkMode()) {
            becDualDataPrePost = clusterIdentificationPreBatchEffectCorrection.concat(
                SC__SCANPY__DIM_REDUCTION__UMAP.out.map { it -> tuple(it[0], it[1], *it[2]) } // get back the parameters stored as inertParameters
            ).map {
                it -> tuple(it[2..(it.size()-1)], it[0], it[1])
            }.groupTuple(
                by: [0, clusteringParams.numParams()-1]
            ).map { 
                it -> tuple(it[1], *it[2]) 
            }
        } else {
            becDualDataPrePost = clusterIdentificationPreBatchEffectCorrection.join(SC__PUBLISH_H5AD.out)
        }
        bbknn_report = GENERATE_DUAL_INPUT_REPORT(
            becDualDataPrePost,
            file(workflow.projectDir + params.sc.scanpy.batch_effect_correct.report_ipynb),
            "SC_BEC_BBKNN_report"
        )

    emit:
        data = SC__SCANPY__DIM_REDUCTION__UMAP.out
        cluster_report = CLUSTER_IDENTIFICATION.out.report
        bbknn_report

}

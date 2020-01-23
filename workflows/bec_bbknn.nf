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

// scanpy:
include '../processes/batch_effect_correct.nf' params(params)

include SC__SCANPY__CLUSTERING from '../processes/cluster.nf' params(params)
include SC__SCANPY__DIM_REDUCTION as SC__SCANPY__DIM_REDUCTION__UMAP from '../processes/dim_reduction.nf' params(params + [method: "umap"])
include CLUSTER_IDENTIFICATION from './cluster_identification.nf' params(params)
include SC__PUBLISH_H5AD from '../../utils/processes/utils.nf' params(params)

// reporting:
include GENERATE_DUAL_INPUT_REPORT from './create_report.nf' params(params + params.global)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow BEC_BBKNN {

    take:
        normalizedTransformedData
        dimReductionData
        clusterIdentifiedWithoutBatchEffectCorrection

    main:
        SC__SCANPY__BATCH_EFFECT_CORRECTION( dimReductionData )
        CLUSTER_IDENTIFICATION(
            normalizedTransformedData,
            SC__SCANPY__BATCH_EFFECT_CORRECTION.out 
        )
        SC__SCANPY__DIM_REDUCTION__UMAP( CLUSTER_IDENTIFICATION.out.marker_genes )
        SC__PUBLISH_H5AD( 
            SC__SCANPY__DIM_REDUCTION__UMAP.out,
            "BEC_BBKNN.output"
        )
        // This will generate a dual report with results from
        // - CLUSTER_IDENTIFICATION without batch effect correction
        // - CLUSTER_IDENTIFICATION with batch effect correction
        bbknn_report = GENERATE_DUAL_INPUT_REPORT(
            clusterIdentifiedWithoutBatchEffectCorrection,
            SC__PUBLISH_H5AD.out,
            file(workflow.projectDir + params.sc.scanpy.batch_effect_correct.report_ipynb),
            "SC_BEC_BBKNN_report"
        )

    emit:
        data = SC__SCANPY__DIM_REDUCTION__UMAP.out
        cluster_report = CLUSTER_IDENTIFICATION.out.report
        bbknn_report

}

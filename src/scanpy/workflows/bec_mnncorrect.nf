/*
 * BEC_MNNCORRECT workflow 
 * - batch effect correction using python package mnnpy (fast and python version of mnnCorrect (Haghverdi et al, 2018)
 */ 

nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    clean;
} from '../../utils/processes/utils.nf' params(params)
include {
    COMBINE_BY_PARAMS;
} from "../../utils/workflows/utils.nf" params(params)
include {
    PUBLISH as PUBLISH_BEC_OUTPUT;
    PUBLISH as PUBLISH_BEC_DIMRED_OUTPUT;
     PUBLISH as PUBLISH_FINAL_HARMONY_OUTPUT;
} from "../../utils/workflows/utils.nf" params(params)

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SCANPY__BATCH_EFFECT_CORRECTION;
} from '../processes/batch_effect_correct.nf' params(params)
include {
    SC__SCANPY__FEATURE_SCALING;
} from '../processes/transform.nf' params(params)
include {
    SC__SCANPY__REGRESS_OUT;
} from '../processes/regress_out.nf' params(params)
include {
    DIM_REDUCTION_PCA;
} from './dim_reduction_pca' params(params + [method: "pca"])
include {
    NEIGHBORHOOD_GRAPH;
} from './neighborhood_graph.nf' params(params)
include {
    DIM_REDUCTION_TSNE_UMAP;
} from './dim_reduction' params(params)
include {
    SC__SCANPY__CLUSTERING_PARAMS;
} from '../processes/cluster.nf' params(params)
include {
    CLUSTER_IDENTIFICATION;
} from './cluster_identification.nf' params(params)
// reporting:
include {
    GENERATE_DUAL_INPUT_REPORT;
} from './create_report.nf' params(params + params.global)


//////////////////////////////////////////////////////
//  Define the workflow 

workflow BEC_MNNCORRECT {

    take:
        normalizedTransformedData
        data
        // Expects (sampleId, anndata)
        clusterIdentificationPreBatchEffectCorrection

    main:
        out = params.getToolParams("scanpy").containsKey("regress_out") 
            ? SC__SCANPY__REGRESS_OUT( hvg ) : data
        SC__SCANPY__BATCH_EFFECT_CORRECTION( 
            out.map { 
                it -> tuple(it[0], it[1], null) 
            }
        )

        PUBLISH_BEC_OUTPUT(
            SC__SCANPY__BATCH_EFFECT_CORRECTION.out,
            "BEC_MNNCORRECT.output",
            "h5ad",
            null,
            false
        )

        SC__SCANPY__FEATURE_SCALING( 
            SC__SCANPY__BATCH_EFFECT_CORRECTION.out.map { 
                it -> tuple(it[0], it[1]) 
            }
        )
        DIM_REDUCTION_PCA( SC__SCANPY__FEATURE_SCALING.out )
        NEIGHBORHOOD_GRAPH( DIM_REDUCTION_PCA.out )

        // Run dimensionality reduction
        DIM_REDUCTION_TSNE_UMAP( NEIGHBORHOOD_GRAPH.out )

        PUBLISH_BEC_DIMRED_OUTPUT(
            DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
            "BEC_HARMONY.dimred_output",
            "h5ad",
            null,
            false
        )

        // Define the parameters for clustering
        def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.getToolParams("scanpy").clustering) )
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

        PUBLISH_FINAL_HARMONY_OUTPUT( 
            marker_genes.map {
                it -> tuple(it[0], it[1], it[2])
            },
            "BEC_MNNCORRECT.final_output",
            "h5ad",
            null,
            clusteringParams.isParameterExplorationModeOn()
        )

        // This will generate a dual report with results from
        // - Pre batch effect correction
        // - Post batch effect correction
        becDualDataPrePost = COMBINE_BY_PARAMS(
            clusterIdentificationPreBatchEffectCorrection,
            // Use PUBLISH output to avoid "input file name collision"
            PUBLISH_FINAL_HARMONY_OUTPUT.out,
            clusteringParams
        )

        mnncorrect_report = GENERATE_DUAL_INPUT_REPORT(
            becDualDataPrePost,
            file(workflow.projectDir + params.getToolParams("scanpy").batch_effect_correct.report_ipynb),
            "SC_BEC_MNNCORRECT_report",
            clusteringParams.isParameterExplorationModeOn()
        )

    emit:
        data = marker_genes
        cluster_report = CLUSTER_IDENTIFICATION.out.report
        mnncorrect_report

}

//
// Version: 
// Test: 
// Command: 
//
/*
 * QC workflow 
 * Source:
 * 
 * Steps considered: 
 * - Aggregating data from multiple motifs (or tracks) SCENIC runs to loom
 */ 

nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include SC__SCENIC__AGGR_MULTI_RUNS_FEATURES from './../processes/aggregateMultiRunsFeatures' params(params)
include SC__SCENIC__AGGR_MULTI_RUNS_REGULONS from './../processes/aggregateMultiRunsRegulons' params(params)
include SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER from './../processes/aucellGeneSigsFromFolder' params(params)
include SC__SCENIC__CONVERT_MULTI_RUNS_FEATURES_TO_REGULONS from './../processes/convertMultiRunsMotifsToRegulons' params(params)
include SC__SCENIC__SAVE_MULTI_RUNS_TO_LOOM from './../processes/saveMultiRunsToLoom' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow AGGREGATE_MULTI_RUNS_TO_LOOM {
    get:
        filteredloom
        ctx
        auc
        type
    main:
        /* Aggregate features (motifs or tracks) from multiple runs */
        ctx_aggr_features = SC__SCENIC__AGGR_MULTI_RUNS_FEATURES(
            ctx.map{ it -> it[1] }.collect(),
            type
        )

        /* Aggregate regulons (motifs or tracks) from multiple runs */
        regulons_folder = SC__SCENIC__AGGR_MULTI_RUNS_REGULONS( 
            auc.collect(),
            type
        )

        /* Run AUCell on aggregated regulons */
        regulons_auc = SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER(
            filteredloom,
            regulons_folder,
            type
        )

        /* Convert aggregated motif enrichment table to regulons */
        aggr_regulons = SC__SCENIC__CONVERT_MULTI_RUNS_FEATURES_TO_REGULONS(
            ctx_aggr_features,
            regulons_folder,
            type
        )

        /* Save multiple motifs (or tracks) SCENIC runs to loom */
        scenic_loom = SC__SCENIC__SAVE_SCENIC_MULTI_RUNS_TO_LOOM( 
            filteredloom,
            aggr_regulons,
            regulons_auc,
            type
        )
    emit:
        scenic_loom
}

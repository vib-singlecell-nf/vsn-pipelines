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

include AGGR_MULTI_RUNS_FEATURES as AGGR_FEATURES from './../processes/multiruns/aggregateFeatures' params(params)
include AGGR_MULTI_RUNS_REGULONS as AGGR_REGULONS from './../processes/multiruns/aggregateRegulons' params(params)
include AUCELL_FROM_FOLDER as AUCELL from './../processes/multiruns/aucellFromFolder' params(params)
include CONVERT_MULTI_RUNS_FEATURES_TO_REGULONS as FEATURES_TO_REGULONS from './../processes/multiruns/convertMotifsToRegulons' params(params)
include SAVE_MULTI_RUNS_TO_LOOM as SAVE_TO_LOOM from './../processes/multiruns/saveToLoom' params(params)

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
        ctx_aggr_features = AGGR_FEATURES(
            ctx.map{ it -> it[1] }.collect(),
            type
        )

        /* Aggregate regulons (motifs or tracks) from multiple runs */
        regulons_folder = AGGR_REGULONS( 
            auc.collect(),
            type
        )

        /* Run AUCell on aggregated regulons */
        regulons_auc = AUCELL(
            filteredloom,
            regulons_folder,
            type
        )

        /* Convert aggregated motif enrichment table to regulons */
        aggr_regulons = FEATURES_TO_REGULONS(
            ctx_aggr_features,
            regulons_folder,
            type
        )

        /* Save multiple motifs (or tracks) SCENIC runs to loom */
        scenic_loom = SAVE_TO_LOOM( 
            filteredloom,
            aggr_regulons,
            regulons_auc,
            type
        )
    emit:
        scenic_loom
}

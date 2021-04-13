nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include {
    AGGR_MULTI_RUNS_FEATURES as AGGR_FEATURES;
 } from './../processes/multiruns/aggregateFeatures' params(params)
include {
    AGGR_MULTI_RUNS_REGULONS as AGGR_REGULONS;
} from './../processes/multiruns/aggregateRegulons' params(params)
include {
    AUCELL_FROM_FOLDER as AUCELL;
} from './../processes/multiruns/aucellFromFolder' params(params)
include {
    CONVERT_MULTI_RUNS_FEATURES_TO_REGULONS as FEATURES_TO_REGULONS;
} from './../processes/multiruns/convertMotifsToRegulons' params(params)
include {
    SAVE_MULTI_RUNS_TO_LOOM as SAVE_TO_LOOM;
} from './../processes/multiruns/saveToLoom' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow AGGREGATE_MULTI_RUNS_TO_LOOM {

    take:
        filteredLoom
        ctx
        auc
        type

    main:
        /* Aggregate features (motifs or tracks) from multiple runs */
        ctxAggrFeatures = AGGR_FEATURES(
            ctx
                // (sampleId, loom, resultFromStep, runId)
                .map { it -> tuple(it[0], it[2]) }
                .groupTuple(),
            type
        )

        /* Aggregate regulons (motifs or tracks) from multiple runs */
        regulonsFolder = AGGR_REGULONS( 
            auc
                // (sampleId, loom, resultFromStep, runId)
                .map { it -> tuple(it[0], it[2]) }
                .groupTuple(),
            type
        )

        /* Run AUCell on aggregated regulons */
        regulonAuc = AUCELL(
            filteredLoom.join(regulonsFolder),
            type
        )

        /* Convert aggregated motif enrichment table to regulons */
        aggrRegulons = FEATURES_TO_REGULONS(
            ctxAggrFeatures.join(regulonsFolder),
            type
        )

        /* Save multiple motifs (or tracks) SCENIC runs to loom */
        scenic_loom = SAVE_TO_LOOM( 
            filteredLoom.join(aggrRegulons).join(regulonAuc),
            type
        )

    emit:
        scenic_loom

}

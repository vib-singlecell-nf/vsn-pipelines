//
// Version:
// Test:
// Command: 
//  nextflow run src/singlecelltxbenchmark/pipelines/bec__bbknn -profile singularity --tenx_folder data/01.count/**/filtered_feature_bc_matrix --project_name tiny
//
/*
 * SCENIC workflow 
 * Source:
 * 
 * Steps considered: 

 */ 
import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2


//////////////////////////////////////////////////////
//  Define the parameters for current testing proces
include SC__SCENIC__GRNBOOST2WITHOUTDASK                        from './processes/grnboost2withoutDask'  params(params)
include SC__SCENIC__CISTARGET as SC__SCENIC__CISTARGET__MOTIF   from './processes/cistarget'             params(params)
include SC__SCENIC__CISTARGET as SC__SCENIC__CISTARGET__TRACK   from './processes/cistarget'             params(params)
include SC__SCENIC__AUCELL as SC__SCENIC__AUCELL__MOTIF         from './processes/aucell'                params(params)
include SC__SCENIC__AUCELL as SC__SCENIC__AUCELL__TRACK         from './processes/aucell'                params(params)
include SC__SCENIC__AGGR_MULTI_RUNS_FEATURES as SC__SCENIC__AGGR_MULTI_RUNS_FEATURES__MOTIF from './processes/aggregateMultiRunsFeatures' params(params)
include SC__SCENIC__AGGR_MULTI_RUNS_FEATURES as SC__SCENIC__AGGR_MULTI_RUNS_FEATURES__TRACK from './processes/aggregateMultiRunsFeatures' params(params)
include SC__SCENIC__AGGR_MULTI_RUNS_REGULONS as SC__SCENIC__AGGR_MULTI_RUNS_REGULONS__MOTIF from './processes/aggregateMultiRunsRegulons' params(params)
include SC__SCENIC__AGGR_MULTI_RUNS_REGULONS as SC__SCENIC__AGGR_MULTI_RUNS_REGULONS__TRACK from './processes/aggregateMultiRunsRegulons' params(params)
include SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER as SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER__MOTIF from './processes/aucellGeneSigsFromFolder' params(params)
include SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER as SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER__TRACK from './processes/aucellGeneSigsFromFolder' params(params)
include SC__SCENIC__SAVE_SCENIC_MULTI_RUNS_TO_LOOM as SC__SCENIC__SAVE_SCENIC_MULTI_RUNS_TO_LOOM_MOTIF from './processes/saveScenicMultiRunsToLoom' params(params)
include SC__SCENIC__SAVE_SCENIC_MULTI_RUNS_TO_LOOM as SC__SCENIC__SAVE_SCENIC_MULTI_RUNS_TO_LOOM_TRACK from './processes/saveScenicMultiRunsToLoom' params(params)
include SC__SCENIC__MERGESCENICLOOMS                            from './processes/mergeScenicLooms'      params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

/*
 * SCENIC workflow
 */ 

// Create channel for the different runs
if(params.sc.scenic.containsKey("numRuns")) {
    runs = Channel.from( 1..params.sc.scenic.numRuns )
} else {
    runs = Channel.from( 1..1 )
}

workflow SCENIC {
    get:
        filteredloom
    main:
        // filteredloom = file(params.sc.scenic.filteredloom)
        /* GRN */
        tfs = file(params.sc.scenic.grn.TFs)
        grn = SC__SCENIC__GRNBOOST2WITHOUTDASK( runs, filteredloom, tfs )

        /* cisTarget 
            motif analysis
        */
        // channel for SCENIC databases resources:
        motifDB = Channel
            .fromPath( params.sc.scenic.cistarget.mtfDB )
            .collect() // use all files together in the ctx command
        motifANN = file(params.sc.scenic.cistarget.mtfANN)
        ctx_mtf = SC__SCENIC__CISTARGET__MOTIF( runs, filteredloom, grn, motifDB, motifANN, 'mtf' )

        /* cisTarget 
            track analysis
        */
        trackDB = Channel
            .fromPath( params.sc.scenic.cistarget.trkDB )
            .collect() // use all files together in the ctx command
        trackANN = file(params.sc.scenic.cistarget.trkANN)
        ctx_trk = SC__SCENIC__CISTARGET__TRACK( runs, filteredloom, grn, trackDB, trackANN, 'trk' )

        /* AUCell, motif regulons */
        auc_mtf = SC__SCENIC__AUCELL__MOTIF( runs, filteredloom, ctx_mtf, 'mtf' )

        /* AUCell, track regulons */
        auc_trk = SC__SCENIC__AUCELL__TRACK( runs, filteredloom, ctx_trk, 'trk' )

        // visualize and merge
        if(params.sc.scenic.containsKey("numRuns") && params.sc.scenic.numRuns > 1) {
            if(params.sc.scenic.numRuns > 2 && params.global.qsubaccount.length() == 0)
                throw new Exception("Consider to run SCENIC in multi-runs mode as jobs. Specify the qsubaccount parameter accordingly.")
            // Aggregate features (motifs and tracks)
            /* Aggregate motifs from multiple runs */
            aggr_features_mtf = SC__SCENIC__AGGR_MULTI_RUNS_FEATURES__MOTIF(
                ctx_mtf.collect(),
                'mtf'
            )
            /* Aggregate tracks from multiple runs */
            aggr_features_trk = SC__SCENIC__AGGR_MULTI_RUNS_FEATURES__TRACK(
                ctx_trk.collect(),
                'trk'
            )

            // Aggregate regulons (motifs and tracks)
            /* Aggregate motif regulons from multiple runs */
            regulons_folder_mtf = SC__SCENIC__AGGR_MULTI_RUNS_REGULONS__MOTIF( 
                auc_mtf.collect(),
                'mtf'
            )
            /* Aggregate track regulons from multiple runs */
            regulons_folder_trk = SC__SCENIC__AGGR_MULTI_RUNS_REGULONS__TRACK(
                auc_trk.collect(),
                'trk'
            )

            // Run AUCell on aggregated regulons
            /* Aggregate motif regulons from multiple runs */
            regulons_auc_mtf = SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER__MOTIF(
                filteredloom,
                regulons_folder_mtf,
                'mtf'
            )
            /* Aggregate track regulons from multiple runs */
            regulons_auc_trk = SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER__TRACK(
                filteredloom,
                regulons_folder_trk,
                'trk'
            )

            // Save to loom
            /* Save multiple motif SCENIC runs to loom*/
            scenic_loom_mtf = SC__SCENIC__SAVE_SCENIC_MULTI_RUNS_TO_LOOM_MOTIF( 
                filteredloom,
                aggr_features_mtf,
                regulons_folder_mtf,
                regulons_auc_mtf,
                'mtf'
            )
            /* Save multiple track SCENIC runs to loom*/
            scenic_loom_trk = SC__SCENIC__SAVE_SCENIC_MULTI_RUNS_TO_LOOM_TRACK( 
                filteredloom,
                aggr_features_trk,
                regulons_folder_trk,
                regulons_auc_trk,
                'trk'
            )
            SC__SCENIC__MERGESCENICLOOMS(
                scenic_loom_mtf,
                scenic_loom_trk
            )
        } else {
            SC__SCENIC__MERGESCENICLOOMS( auc_mtf, auc_trk )
        }

}

// Uncomment to test
workflow {
    main:
        SCENIC( file( params.sc.scenic.filteredloom ) )
}


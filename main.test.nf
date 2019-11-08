//
// Tests
//  

// Test 1: SC__SCENIC__GRNBOOST2WITHOUTDASK (from processes/)
// Time: ~2min
// Command: 
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf -profile singularity --test SC__SCENIC__GRNBOOST2WITHOUTDASK

// Test 2: SC__SCENIC__CISTARGET (from processes/)
// Time: ~10min
// Command:
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf -profile singularity --test SC__SCENIC__CISTARGET

// Test 3: SC__SCENIC__AUCELL (from processes/)
// Time: ~1min
// Command:
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf -profile singularity --test SC__SCENIC__AUCELL

// Test 4: SC__SCENIC__AGGR_MULTI_RUNS_FEATURES (from processes/)
// Time: ~1min
// Command:
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf -profile singularity --test SC__SCENIC__AGGR_MULTI_RUNS_FEATURES

// Test 5: SC__SCENIC__AGGR_MULTI_RUNS_REGULONS (from processes/)
// Time: <1min
// Command:
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf -profile singularity --test SC__SCENIC__AGGR_MULTI_RUNS_REGULONS

// Test 6: SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER (from processes/)
// Time: ~?min
// Command:
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf -profile singularity --test SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER

// Test 7: SC__SCENIC__SAVE_SCENIC_MULTI_RUNS_TO_LOOM (from processes/)
// Time: ~?min
// Command:
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf --test SC__SCENIC__SAVE_SCENIC_MULTI_RUNS_TO_LOOM


nextflow.preview.dsl=2

///////////////////////////////////////////
//  Define the parameters for all processes

include SC__SCENIC__GRNBOOST2WITHOUTDASK from './processes/grnboost2withoutDask' params(params)
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
include SC__SCENIC__PUBLISH_LOOM            from './processes/scenicLoomHandler'     params(params)
include SC__SCENIC__VISUALIZE               from './processes/scenicLoomHandler'     params(params)

// Create channel for the different runs
runs = Channel.from( 1..params.sc.scenic.numRuns )

// Make the test workflow 
workflow test_SC__SCENIC__GRNBOOST2WITHOUTDASK {
    get:
        loom
    main:
        tfs = file(params.sc.scenic.grn.TFs)
        SC__SCENIC__GRNBOOST2WITHOUTDASK( runs, loom, tfs )
    emit:
        SC__SCENIC__GRNBOOST2WITHOUTDASK.out
}

// Make the test workflow 
workflow test_SC__SCENIC__CISTARGET {
    get:
        filteredloom
        grn
    main:
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
    emit:
        ctx_mtf
        ctx_trk
}

// Make the test workflow 
workflow test_SC__SCENIC__AUCELL {
    get:
        filteredloom
        ctx_mtf
        ctx_trk
    main:
        /* AUCell, motif regulons */
        auc_mtf = SC__SCENIC__AUCELL__MOTIF( runs, filteredloom, ctx_mtf, 'mtf' )

        /* AUCell, track regulons */
        auc_trk = SC__SCENIC__AUCELL__TRACK( runs, filteredloom, ctx_trk, 'trk' )
    emit:
        auc_mtf
        auc_trk
}

// Make the test workflow 
workflow test_SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER {
    get:
        filteredloom
        regulons_folder_mtf
        regulons_folder_trk
    main:
        /* Aggregate motif regulons from multiple runs */
        regulons_auc_mtf = SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER__MOTIF( filteredloom, regulons_folder_mtf, 'mtf' )

        /* Aggregate track regulons from multiple runs */
        regulons_auc_trk = SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER__TRACK( filteredloom, regulons_folder_trk, 'trk' )
    emit:
        regulons_auc_mtf
        regulons_auc_trk
}

workflow {
    main:
        switch(params.test) {
            case "SC__SCENIC__GRNBOOST2WITHOUTDASK":
                test_SC__SCENIC__GRNBOOST2WITHOUTDASK( file( params.sc.scenic.filteredloom ) )
            break;
            case "SC__SCENIC__CISTARGET":
                grn = Channel.fromPath(params.sc.scenic.scenicoutdir + "/grnboost2withoutDask/run_*/run_*__adj.tsv")
                test_SC__SCENIC__CISTARGET( file( params.sc.scenic.filteredloom ), grn )
            break;
            case "SC__SCENIC__AUCELL":
                ctx_mtf = Channel.fromPath(params.sc.scenic.scenicoutdir + "/cistarget/run_*/run_*__reg_mtf.csv")
                ctx_trk = Channel.fromPath(params.sc.scenic.scenicoutdir + "/cistarget/run_*/run_*__reg_trk.csv")
                test_SC__SCENIC__AUCELL( file( params.sc.scenic.filteredloom ), ctx_mtf, ctx_trk )
            break;
            case "SC__SCENIC__AGGR_MULTI_RUNS_FEATURES":
                /* Aggregate motifs from multiple runs */
                reg_mtf = Channel.fromPath(params.sc.scenic.scenicoutdir + "/cistarget/run_*/run_*__reg_mtf.csv")
                SC__SCENIC__AGGR_MULTI_RUNS_FEATURES__MOTIF( reg_mtf.collect(), 'mtf' )
                if(params.sc.scenic.cistarget.trkDB) {
                    /* Aggregate tracks from multiple runs */
                    reg_trk = Channel.fromPath(params.sc.scenic.scenicoutdir + "/cistarget/run_*/run_*__reg_trk.csv")
                    SC__SCENIC__AGGR_MULTI_RUNS_FEATURES__TRACK( reg_trk.collect(), 'trk' )
                }
            break;
            case "SC__SCENIC__AGGR_MULTI_RUNS_REGULONS":
                /* Aggregate motif regulons from multiple runs */
                auc_mtf_looms = Channel.fromPath(params.sc.scenic.scenicoutdir + "/aucell/run_*/run_*__auc_mtf.loom")
                SC__SCENIC__AGGR_MULTI_RUNS_REGULONS__MOTIF( auc_mtf_looms.collect(), 'mtf' )
                if(params.sc.scenic.cistarget.trkDB) {
                    /* Aggregate track regulons from multiple runs */
                    auc_trk_looms = Channel.fromPath(params.sc.scenic.scenicoutdir + "/aucell/run_*/run_*__auc_trk.loom")
                    SC__SCENIC__AGGR_MULTI_RUNS_REGULONS__TRACK( auc_trk_looms.collect(), 'trk' )
                }
            break;
            case "SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER":
                /* Aggregate motif regulons from multiple runs */
                regulons_folder_mtf = file(params.sc.scenic.scenicoutdir + "/multi_runs_regulons_mtf")
                SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER__MOTIF( file(params.sc.scenic.filteredloom), regulons_folder_mtf, 'mtf' )
                if(params.sc.scenic.cistarget.trkDB) {
                    /* Aggregate track regulons from multiple runs */
                    regulons_folder_trk = file(params.sc.scenic.scenicoutdir + "/multi_runs_regulons_trk")
                    SC__SCENIC__AUCELL_GENESIGS_FROM_FOLDER__TRACK( file(params.sc.scenic.filteredloom), regulons_folder_trk, 'trk' )
                }
            break;
            case "SC__SCENIC__SAVE_SCENIC_MULTI_RUNS_TO_LOOM":
                filteredloom = file(params.sc.scenic.filteredloom)
                aggr_features_mtf = file(params.sc.scenic.scenicoutdir + "/multi_runs_cistarget/multi_runs_features_mtf.csv.gz")
                regulons_folder_mtf = file(params.sc.scenic.scenicoutdir + "/multi_runs_regulons_mtf")
                regulons_auc_mtf = file(params.sc.scenic.scenicoutdir + "/multi_runs_aucell/multi_runs_regulons_auc_mtf.tsv")
                
                /* Save multiple motif SCENIC runs to loom*/
                SC__SCENIC__SAVE_SCENIC_MULTI_RUNS_TO_LOOM_MOTIF( 
                    filteredloom,
                    aggr_features_mtf,
                    regulons_folder_mtf,
                    regulons_auc_mtf,
                    'mtf' 
                )
                if(params.sc.scenic.cistarget.trkDB) {
                    regulons_folder_trk = file(params.sc.scenic.scenicoutdir + "/multi_runs_regulons_trk")
                    aggr_features_trk = file(params.sc.scenic.scenicoutdir + "/multi_runs_cistarget/multi_runs_features_trk.csv.gz")
                    regulons_auc_trk = file(params.sc.scenic.scenicoutdir + "/multi_runs_aucell/multi_runs_regulons_auc_trk.tsv")
                    /* Save multiple track SCENIC runs to loom*/
                    SC__SCENIC__SAVE_SCENIC_MULTI_RUNS_TO_LOOM_TRACK( 
                        filteredloom,
                        aggr_features_trk,
                        regulons_folder_trk,
                        regulons_auc_trk,
                        'trk' 
                    )
                }
            break;
            case "SC__SCENIC__VISUALIZE_PUBLISH":
                /* Aggregate motif regulons from multiple runs */
                multi_runs_aucell_mtf_loom = file(params.sc.scenic.scenicoutdir + "/multi_runs_looms/multi_runs_regulons_auc_mtf.loom")
                SC__SCENIC__PUBLISH_LOOM( SC__SCENIC__VISUALIZE( multi_runs_aucell_mtf_loom ) )
            break;
            default:
                throw new Exception("The test parameters should be specified.")
            break;
        }
}
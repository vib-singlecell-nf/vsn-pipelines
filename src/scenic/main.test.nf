//
// Tests
// CURRENTLY DEPRECATED 
// TODO: NEEDS TO BE UPDATED USING NEW API
//  

// Test 1: GRNBOOST2WITHOUTDASK (from processes/)
// Time: ~2min
// Command: 
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf -profile singularity --test GRNBOOST2WITHOUTDASK

// Test 2: CISTARGET (from processes/)
// Time: ~10min
// Command:
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf -profile singularity --test CISTARGET

// Test 3: AUCELL (from processes/)
// Time: ~1min
// Command:
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf -profile singularity --test AUCELL

// Test 4: AGGR_MULTI_RUNS_FEATURES (from processes/)
// Time: ~1min
// Command:
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf -profile singularity --test AGGR_MULTI_RUNS_FEATURES

// Test 5: AGGR_MULTI_RUNS_REGULONS (from processes/)
// Time: <1min
// Command:
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf -profile singularity --test AGGR_MULTI_RUNS_REGULONS

// Test 6: AUCELL_FROM_FOLDER (from processes/)
// Time: ~?min
// Command:
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf -profile singularity --test AUCELL_FROM_FOLDER

// Test 7: SAVE_SCENIC_MULTI_RUNS_TO_LOOM (from processes/)
// Time: ~?min
// Command:
//  nextflow -C conf/multi_runs.config,conf/test.config run main.test.nf --test SAVE_SCENIC_MULTI_RUNS_TO_LOOM


nextflow.enable.dsl=2

///////////////////////////////////////////
//  Define the parameters for all processes

include {
    GRNBOOST2_WITHOUT_DASK;
} from './processes/grnboost2withoutDask' params(params)
include {
    CISTARGET as CISTARGET__MOTIF;
    CISTARGET as CISTARGET__TRACK;
}   from './processes/cistarget'             params(params)
include {
    AUCELL as AUCELL__MOTIF;
    AUCELL as AUCELL__TRACK;
} from './processes/aucell'                params(params)
include {
    AGGR_MULTI_RUNS_FEATURES as AGGR_MULTI_RUNS_FEATURES__MOTIF;
    AGGR_MULTI_RUNS_FEATURES as AGGR_MULTI_RUNS_FEATURES__TRACK;
} from './processes/multiruns/aggregateFeatures' params(params)
include {
    AGGR_MULTI_RUNS_REGULONS as AGGR_MULTI_RUNS_REGULONS__MOTIF;
    AGGR_MULTI_RUNS_REGULONS as AGGR_MULTI_RUNS_REGULONS__TRACK;
} from './processes/multiruns/aggregateMultiRunsRegulons' params(params)
include {
    AUCELL_FROM_FOLDER as AUCELL_FROM_FOLDER__MOTIF;
    AUCELL_FROM_FOLDER as AUCELL_FROM_FOLDER__TRACK;
} from './processes/aucellFromFolder' params(params)
include {
    SAVE_SCENIC_MULTI_RUNS_TO_LOOM as SAVE_SCENIC_MULTI_RUNS_TO_LOOM_MOTIF;
    SAVE_SCENIC_MULTI_RUNS_TO_LOOM as SAVE_SCENIC_MULTI_RUNS_TO_LOOM_TRACK;
} from './processes/saveScenicMultiRunsToLoom' params(params)
include {
    MERGE_MOTIF_TRACK_LOOMS;
    PUBLISH_LOOM;
    VISUALIZE;
} from './processes/loomHandler'     params(params)

// Create channel for the different runs
runs = Channel.from( 1..params.tools.scenic.numRuns )

// Make the test workflow 
workflow test_GRNBOOST2WITHOUTDASK {

    take:
        loom

    main:
        tfs = file(params.tools.scenic.grn.TFs)
        GRNBOOST2WITHOUTDASK( runs, loom, tfs )

    emit:
        GRNBOOST2WITHOUTDASK.out

}

// Make the test workflow 
workflow test_CISTARGET {

    take:
        filteredloom
        grn

    main:
        // channel for SCENIC databases resources:
        motifDB = Channel
            .fromPath( params.tools.scenic.cistarget.mtfDB )
            .collect() // use all files together in the ctx command
        motifANN = file(params.tools.scenic.cistarget.mtfANN)
        ctx_mtf = CISTARGET__MOTIF( runs, filteredloom, grn, motifDB, motifANN, 'mtf' )

        /* cisTarget 
            track analysis
        */
        trackDB = Channel
            .fromPath( params.tools.scenic.cistarget.trkDB )
            .collect() // use all files together in the ctx command
        trackANN = file(params.tools.scenic.cistarget.trkANN)
        ctx_trk = CISTARGET__TRACK( runs, filteredloom, grn, trackDB, trackANN, 'trk' )

    emit:
        ctx_mtf
        ctx_trk

}

// Make the test workflow 
workflow test_AUCELL {

    take:
        filteredloom
        ctx_mtf
        ctx_trk

    main:
        /* AUCell, motif regulons */
        auc_mtf = AUCELL__MOTIF( runs, filteredloom, ctx_mtf, 'mtf' )

        /* AUCell, track regulons */
        auc_trk = AUCELL__TRACK( runs, filteredloom, ctx_trk, 'trk' )

    emit:
        auc_mtf
        auc_trk

}

// Make the test workflow 
workflow test_SINGLE_RUN_BY_ID {

    take:
        runId

    main:
        filteredloom = file( params.tools.scenic.filteredloom )
        tfs = file(params.tools.scenic.grn.TFs)
        run = Channel.from( runId..runId )
        grn = GRNBOOST2WITHOUTDASK( run, filteredloom, tfs )
        // channel for SCENIC databases resources:
        motifDB = Channel
            .fromPath( params.tools.scenic.cistarget.mtfDB )
            .collect() // use all files together in the ctx command
        motifANN = file(params.tools.scenic.cistarget.mtfANN)
        ctx_mtf = CISTARGET__MOTIF( run, filteredloom, grn, motifDB, motifANN, 'mtf' )
        /* AUCell, motif regulons */
        auc_mtf = AUCELL__MOTIF( run, filteredloom, ctx_mtf, 'mtf' )

    emit:
        auc_mtf

}

// Make the test workflow 
workflow test_AUCELL_FROM_FOLDER {

    take:
        filteredloom
        regulons_folder_mtf
        regulons_folder_trk

    main:
        /* Aggregate motif regulons from multiple runs */
        regulons_auc_mtf = AUCELL_FROM_FOLDER__MOTIF( filteredloom, regulons_folder_mtf, 'mtf' )

        /* Aggregate track regulons from multiple runs */
        regulons_auc_trk = AUCELL_FROM_FOLDER__TRACK( filteredloom, regulons_folder_trk, 'trk' )

    emit:
        regulons_auc_mtf
        regulons_auc_trk

}

workflow {

    main:
        switch(params.test) {
            case "SC__SCENIC_SINGLE_RUN_BY_ID":
                test_SINGLE_RUN_BY_ID( params.runId )
            break;
            case "GRNBOOST2WITHOUTDASK":
                test_GRNBOOST2WITHOUTDASK( file( params.tools.scenic.filteredloom ) )
            break;
            case "CISTARGET":
                grn = Channel.fromPath(params.tools.scenic.scenicoutdir + "/grnboost2withoutDask/run_*/run_*__adj.tsv")
                test_CISTARGET( file( params.tools.scenic.filteredloom ), grn )
            break;
            case "AUCELL":
                ctx_mtf = Channel.fromPath(params.tools.scenic.scenicoutdir + "/cistarget/run_*/run_*__reg_mtf.csv")
                ctx_trk = Channel.fromPath(params.tools.scenic.scenicoutdir + "/cistarget/run_*/run_*__reg_trk.csv")
                test_AUCELL( file( params.tools.scenic.filteredloom ), ctx_mtf, ctx_trk )
            break;
            case "AGGR_MULTI_RUNS_FEATURES":
                /* Aggregate motifs from multiple runs */
                reg_mtf = Channel.fromPath(params.tools.scenic.scenicoutdir + "/cistarget/run_*/run_*__reg_mtf.csv")
                AGGR_MULTI_RUNS_FEATURES__MOTIF( reg_mtf.collect(), 'mtf' )
                if(params.tools.scenic.cistarget.trkDB) {
                    /* Aggregate tracks from multiple runs */
                    reg_trk = Channel.fromPath(params.tools.scenic.scenicoutdir + "/cistarget/run_*/run_*__reg_trk.csv")
                    AGGR_MULTI_RUNS_FEATURES__TRACK( reg_trk.collect(), 'trk' )
                }
            break;
            case "AGGR_MULTI_RUNS_REGULONS":
                /* Aggregate motif regulons from multiple runs */
                auc_mtf_looms = Channel.fromPath(params.tools.scenic.scenicoutdir + "/aucell/run_*/run_*__auc_mtf.loom")
                AGGR_MULTI_RUNS_REGULONS__MOTIF( auc_mtf_looms.collect(), 'mtf' )
                if(params.tools.scenic.cistarget.trkDB) {
                    /* Aggregate track regulons from multiple runs */
                    auc_trk_looms = Channel.fromPath(params.tools.scenic.scenicoutdir + "/aucell/run_*/run_*__auc_trk.loom")
                    AGGR_MULTI_RUNS_REGULONS__TRACK( auc_trk_looms.collect(), 'trk' )
                }
            break;
            case "AUCELL_FROM_FOLDER":
                /* Aggregate motif regulons from multiple runs */
                regulons_folder_mtf = file(params.tools.scenic.scenicoutdir + "/multi_runs_regulons_mtf")
                AUCELL_FROM_FOLDER__MOTIF( file(params.tools.scenic.filteredloom), regulons_folder_mtf, 'mtf' )
                if(params.tools.scenic.cistarget.trkDB) {
                    /* Aggregate track regulons from multiple runs */
                    regulons_folder_trk = file(params.tools.scenic.scenicoutdir + "/multi_runs_regulons_trk")
                    AUCELL_FROM_FOLDER__TRACK( file(params.tools.scenic.filteredloom), regulons_folder_trk, 'trk' )
                }
            break;
            case "SAVE_SCENIC_MULTI_RUNS_TO_LOOM_MOTIF":
                filteredloom = file(params.tools.scenic.filteredloom)
                aggr_features_mtf = file(params.tools.scenic.scenicoutdir + "/multi_runs_cistarget/multi_runs_features_mtf.csv.gz")
                regulons_folder_mtf = file(params.tools.scenic.scenicoutdir + "/multi_runs_regulons_mtf")
                regulons_auc_mtf = file(params.tools.scenic.scenicoutdir + "/multi_runs_aucell/multi_runs_regulons_auc_mtf.tsv")
                
                /* Save multiple motif SCENIC runs to loom*/
                SAVE_SCENIC_MULTI_RUNS_TO_LOOM_MOTIF( 
                    filteredloom,
                    aggr_features_mtf,
                    regulons_folder_mtf,
                    regulons_auc_mtf,
                    'mtf' 
                )
            break;
            case "SAVE_SCENIC_MULTI_RUNS_TO_LOOM_TRACK":
                filteredloom = file(params.tools.scenic.filteredloom)
                regulons_folder_trk = file(params.tools.scenic.scenicoutdir + "/multi_runs_regulons_trk")
                aggr_features_trk = file(params.tools.scenic.scenicoutdir + "/multi_runs_cistarget/multi_runs_features_trk.csv.gz")
                regulons_auc_trk = file(params.tools.scenic.scenicoutdir + "/multi_runs_aucell/multi_runs_regulons_auc_trk.tsv")
                /* Save multiple track SCENIC runs to loom*/
                SAVE_SCENIC_MULTI_RUNS_TO_LOOM_TRACK( 
                    filteredloom,
                    aggr_features_trk,
                    regulons_folder_trk,
                    regulons_auc_trk,
                    'trk' 
                )
            break;
            case "MERGE_MOTIF_TRACK_LOOMS":
                scenic_loom_mtf = file( params.tools.scenic.scenicoutdir + "/multi_runs_looms/multi_runs_regulons_auc_mtf.loom" )
                scenic_loom_trk = file( params.tools.scenic.scenicoutdir + "/multi_runs_looms/multi_runs_regulons_auc_trk.loom" )
                MERGE_MOTIF_TRACK_LOOMS(
                    scenic_loom_mtf,
                    scenic_loom_trk
                )
            break;
            case "VISUALIZE_PUBLISH":
                /* Aggregate motif regulons from multiple runs */
                scenic_loom = file( params.tools.scenic.scenicoutdir + "/" + params.tools.scenic.scenicOutputLoom )
                PUBLISH_LOOM( VISUALIZE( scenic_loom ) )
            break;
            default:
                throw new Exception("The test parameters should be specified.")
            break;
        }

}

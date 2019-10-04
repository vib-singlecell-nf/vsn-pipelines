//
// Tests
//  

// Test 1: SC__SCENIC__GRNBOOST2WITHOUTDASK (from processes/)
// Time: ~3min
// Command: 
//  nextflow -C conf/test.config,scenic.config run main.test.nf -profile singularity --test SC__SCENIC__GRNBOOST2WITHOUTDASK

// Test 2: SC__SCENIC__CISTARGET (from processes/)
// Time: ~7min
// Command:
//  nextflow -C conf/test.config,scenic.config run main.test.nf -profile singularity --test SC__SCENIC__CISTARGET

// Test 3: SC__SCENIC__AUCELL (from processes/)
// Time: ~1min
// Command:
//  nextflow -C conf/test.config,scenic.config run main.test.nf -profile singularity --test SC__SCENIC__AUCELL

// Test 4: SC__SCENIC__MULTI_RUNS_AGGR_REGULONS (from processes/)
// Time: <1min
// Command:
//  nextflow -C conf/test.config,scenic.config run main.test.nf -profile singularity --test SC__SCENIC__MULTI_RUNS_AGGR_REGULONS


nextflow.preview.dsl=2

///////////////////////////////////////////
//  Define the parameters for all processes

include SC__SCENIC__GRNBOOST2WITHOUTDASK from './processes/grnboost2withoutDask' params(params)
include SC__SCENIC__CISTARGET as SC__SCENIC__CISTARGET__MOTIF   from './processes/cistarget'             params(params)
include SC__SCENIC__CISTARGET as SC__SCENIC__CISTARGET__TRACK   from './processes/cistarget'             params(params)
include SC__SCENIC__AUCELL as SC__SCENIC__AUCELL__MOTIF         from './processes/aucell'                params(params)
include SC__SCENIC__AUCELL as SC__SCENIC__AUCELL__TRACK         from './processes/aucell'                params(params)
include SC__SCENIC__MULTI_RUNS_AGGR_REGULONS as SC__SCENIC__MULTI_RUNS_AGGR_REGULONS__MOTIF from './processes/multiRunsAggrRegulons' params(params)
include SC__SCENIC__MULTI_RUNS_AGGR_REGULONS as SC__SCENIC__MULTI_RUNS_AGGR_REGULONS__TRACK from './processes/multiRunsAggrRegulons' params(params)

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
workflow test_SC__SCENIC__MULTI_RUNS_AGGR_REGULONS {
    get:
        mtf_auc_looms
        trk_auc_looms
    main:
        /* Aggregate motif regulons from multiple runs */
        aggr_regulons_mtf = SC__SCENIC__MULTI_RUNS_AGGR_REGULONS__MOTIF( mtf_auc_looms, 'mtf' )

        /* Aggregate track regulons from multiple runs */
        aggr_regulons_trk = SC__SCENIC__MULTI_RUNS_AGGR_REGULONS__TRACK( trk_auc_looms, 'trk' )
    emit:
        aggr_regulons_mtf
        aggr_regulons_trk
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
            case "SC__SCENIC__MULTI_RUNS_AGGR_REGULONS":
                aucell_mtf_looms = Channel.fromPath(params.sc.scenic.scenicoutdir + "/aucell/run_*/run_*__auc_mtf.loom")
                aucell_trk_looms = Channel.fromPath(params.sc.scenic.scenicoutdir + "/aucell/run_*/run_*__auc_trk.loom")
                test_SC__SCENIC__MULTI_RUNS_AGGR_REGULONS(aucell_mtf_looms.collect(), aucell_trk_looms.collect())
            break;
            default:
                throw new Exception("The test parameters should be specified.")
            break;
        }
}
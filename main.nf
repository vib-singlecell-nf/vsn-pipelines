//
// Version:
// Test:
// Command: 
//
/*
 * SCENIC workflow 
 * Source:
 * 
 * Steps considered: 

 */ 
nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces
include SC__SCENIC__GRNBOOST2_WITHOUT_DASK                                                  from './processes/grnboost2withoutDask'  params(params)
include SC__SCENIC__CISTARGET as SC__SCENIC__CISTARGET__MOTIF                               from './processes/cistarget'             params(params)
include SC__SCENIC__CISTARGET as SC__SCENIC__CISTARGET__TRACK                               from './processes/cistarget'             params(params)
include SC__SCENIC__AUCELL as SC__SCENIC__AUCELL__MOTIF                                     from './processes/aucell'                params(params)
include SC__SCENIC__AUCELL as SC__SCENIC__AUCELL__TRACK                                     from './processes/aucell'                params(params)
include AGGREGATE_MULTI_RUNS_TO_LOOM as AGGREGATE_MULTI_RUNS_TO_LOOM__MOTIF   from './workflows/aggregateMultiRuns'    params(params)
include AGGREGATE_MULTI_RUNS_TO_LOOM as AGGREGATE_MULTI_RUNS_TO_LOOM__TRACK   from './workflows/aggregateMultiRuns'    params(params)
include SC__SCENIC__PUBLISH_LOOM                                                            from './processes/scenicLoomHandler'     params(params)
include SC__SCENIC__MERGE_MOTIF_TRACK_LOOMS                                                 from './processes/scenicLoomHandler'     params(params)
include SC__SCENIC__APPEND_SCENIC_LOOM                                                      from './processes/scenicLoomHandler'     params(params)
include SC__SCENIC__VISUALIZE                                                               from './processes/scenicLoomHandler'     params(params)

// reporting:
include './processes/reports.nf' params(params + params.global)

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
        /* GRN */
        tfs = file(params.sc.scenic.grn.TFs)
        grn = SC__SCENIC__GRNBOOST2_WITHOUT_DASK( runs, filteredloom, tfs )

        /* cisTarget motif analysis */
        // channel for SCENIC databases resources:
        motifDB = Channel
            .fromPath( params.sc.scenic.cistarget.mtfDB )
            .collect() // use all files together in the ctx command
        motifANN = file(params.sc.scenic.cistarget.mtfANN)
        ctx_mtf = SC__SCENIC__CISTARGET__MOTIF( grn, filteredloom, motifDB, motifANN, 'mtf' )

        /* cisTarget track analysis */
        if(params.sc.scenic.cistarget.trkDB) {
            trackDB = Channel
                .fromPath( params.sc.scenic.cistarget.trkDB )
                .collect() // use all files together in the ctx command
            trackANN = file(params.sc.scenic.cistarget.trkANN)
            ctx_trk = SC__SCENIC__CISTARGET__TRACK( grn, filteredloom, trackDB, trackANN, 'trk' )
        }

        /* AUCell, motif regulons */
        auc_mtf = SC__SCENIC__AUCELL__MOTIF( ctx_mtf, filteredloom, 'mtf' )

        if(params.sc.scenic.cistarget.trkDB) {
            /* AUCell, track regulons */
            auc_trk = SC__SCENIC__AUCELL__TRACK( ctx_trk, filteredloom, 'trk' )
        }

        // multi-runs aggregation:
        if(params.sc.scenic.containsKey("numRuns") && params.sc.scenic.numRuns > 1) {
            if(params.sc.scenic.numRuns > 2 && params.global.qsubaccount.length() == 0)
                throw new Exception("Consider to run SCENIC in multi-runs mode as jobs. Specify the qsubaccount parameter accordingly.")
            
            scenic_loom_mtf = AGGREGATE_MULTI_RUNS_TO_LOOM__MOTIF(
                filteredloom,
                ctx_mtf,
                auc_mtf,
                'mtf'
            )
            if(params.sc.scenic.cistarget.trkDB) {
                scenic_loom_trk = AGGREGATE_MULTI_RUNS_TO_LOOM__TRACK(
                    filteredloom,
                    ctx_trk,
                    auc_trk,
                    'trk'
                )
                SC__SCENIC__MERGE_MOTIF_TRACK_LOOMS(
                    scenic_loom_mtf,
                    scenic_loom_trk
                )
                out = SC__SCENIC__VISUALIZE(SC__SCENIC__MERGE_MOTIF_TRACK_LOOMS.out)
            } else {
                out = SC__SCENIC__VISUALIZE(scenic_loom_mtf)
            }
        } else {
            if(params.sc.scenic.cistarget.trkDB) {
                out = SC__SCENIC__VISUALIZE(
                    SC__SCENIC__MERGE_MOTIF_TRACK_LOOMS(
                        auc_mtf,
                        auc_trk
                    ))
            } else {
                out = SC__SCENIC__VISUALIZE(auc_mtf)
            }
        }
        SC__SCENIC__PUBLISH_LOOM(out)
    emit:
        out
}


workflow SCENIC_append {
    get:
        filteredloom
        scopeloom
    main:
        scenicloom = SCENIC( filteredloom )
        SC__SCENIC__APPEND_SCENIC_LOOM( scopeloom, scenicloom )
        report_notebook = SC__SCENIC__GENERATE_REPORT(
            file(workflow.projectDir + params.sc.scenic.report_ipynb),
            SC__SCENIC__APPEND_SCENIC_LOOM.out,
            "SCENIC_report"
        )
        SC__SCENIC__REPORT_TO_HTML(report_notebook)
    emit:
        SC__SCENIC__APPEND_SCENIC_LOOM.out
}


// Uncomment to test
workflow {
    main:
        SCENIC( file( params.sc.scenic.filteredloom ) )
}


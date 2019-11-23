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
include GRNBOOST2_WITHOUT_DASK                                    from './processes/grnboost2withoutDask'  params(params)
include CISTARGET as CISTARGET__MOTIF                             from './processes/cistarget'             params(params)
include CISTARGET as CISTARGET__TRACK                             from './processes/cistarget'             params(params)
include AUCELL as AUCELL__MOTIF                                   from './processes/aucell'                params(params)
include AUCELL as AUCELL__TRACK                                   from './processes/aucell'                params(params)
include AGGREGATE_MULTI_RUNS_TO_LOOM as MULTI_RUNS_TO_LOOM__MOTIF from './workflows/aggregateMultiRuns'    params(params)
include AGGREGATE_MULTI_RUNS_TO_LOOM as MULTI_RUNS_TO_LOOM__TRACK from './workflows/aggregateMultiRuns'    params(params)
include PUBLISH_LOOM                                              from './processes/loomHandler'     params(params)
include MERGE_MOTIF_TRACK_LOOMS                                   from './processes/loomHandler'     params(params)
include APPEND_SCENIC_LOOM                                        from './processes/loomHandler'     params(params)
include VISUALIZE                                                 from './processes/loomHandler'     params(params)

// reporting:
include './processes/reports.nf' params(params)

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
        // Expects (sampleId, loom)
        filteredLoom
    main:
        /* GRN */
        filteredLoom.view()
        tfs = file(params.sc.scenic.grn.TFs)
        grn = GRNBOOST2_WITHOUT_DASK( filteredLoom.combine(runs), tfs )

        /* cisTarget motif analysis */
        // channel for SCENIC databases resources:
        motifDB = Channel
            .fromPath( params.sc.scenic.cistarget.mtfDB )
            .collect() // use all files together in the ctx command
        motifANN = file(params.sc.scenic.cistarget.mtfANN)
        ctx_mtf = CISTARGET__MOTIF( grn, motifDB, motifANN, 'mtf' )

        /* cisTarget track analysis */
        if(params.sc.scenic.cistarget.trkDB) {
            trackDB = Channel
                .fromPath( params.sc.scenic.cistarget.trkDB )
                .collect() // use all files together in the ctx command
            trackANN = file(params.sc.scenic.cistarget.trkANN)
            ctx_trk = CISTARGET__TRACK( grn, trackDB, trackANN, 'trk' )
        }

        /* AUCell, motif regulons */
        auc_mtf = AUCELL__MOTIF( ctx_mtf, 'mtf' )

        if(params.sc.scenic.cistarget.trkDB) {
            /* AUCell, track regulons */
            auc_trk = AUCELL__TRACK( ctx_trk, 'trk' )
        }

        // multi-runs aggregation:
        if(params.sc.scenic.containsKey("numRuns") && params.sc.scenic.numRuns > 1) {
            if(params.sc.scenic.numRuns > 2 && params.global.qsubaccount.length() == 0)
                throw new Exception("Consider to run SCENIC in multi-runs mode as jobs. Specify the qsubaccount parameter accordingly.")
            
            scenic_loom_mtf = MULTI_RUNS_TO_LOOM__MOTIF(
                filteredLoom,
                ctx_mtf,
                auc_mtf,
                'mtf'
            )
            if(params.sc.scenic.cistarget.trkDB) {
                scenic_loom_trk = MULTI_RUNS_TO_LOOM__TRACK(
                    filteredLoom,
                    ctx_trk,
                    auc_trk,
                    'trk'
                )
                MERGE_MOTIF_TRACK_LOOMS(
                    scenic_loom_mtf.join(scenic_loom_trk)
                )
                out = VISUALIZE(MERGE_MOTIF_TRACK_LOOMS.out)
            } else {
                out = VISUALIZE(scenic_loom_mtf)
            }
        } else {
            if(params.sc.scenic.cistarget.trkDB) {
                out = VISUALIZE(
                    MERGE_MOTIF_TRACK_LOOMS(
                        auc_mtf
                            .map { it -> tuple(it[0], it[2]) }
                            .join(auc_trk.map { it -> tuple(it[0], it[2]) })
                    ))
            } else {
                out = VISUALIZE(
                    auc_mtf.map { it -> tuple(it[0], it[2]) }
                )
            }
        }
        PUBLISH_LOOM(out)
    emit:
        out
}


workflow SCENIC_append {
    get:
        filteredLoom
        scopeLoom
    main:
        scenicLoom = SCENIC( filteredLoom )
        APPEND_SCENIC_LOOM( scopeLoom.join(scenicLoom) )
        report_notebook = GENERATE_REPORT(
            file(workflow.projectDir + params.sc.scenic.report_ipynb),
            APPEND_SCENIC_LOOM.out,
            "SCENIC_report"
        )
        REPORT_TO_HTML(report_notebook)
    emit:
        APPEND_SCENIC_LOOM.out
}


// Uncomment to test
workflow {
    main:
        SCENIC( Channel.of( tuple("foobar", file(params.sc.scenic.filteredLoom)) ) )
}

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
include './../utils/processes/config'
resolveParams(params, true)

isAppendOnlyMode = params.sc.scenic.containsKey("existingScenicLoom")

def ALLOWED_GENOME_ASSEMBLIES = ['dm6','hg19','hg38', 'mm10']

//////////////////////////////////////////////////////
//  Sanity checks
if(!isAppendOnlyMode && !params.global.containsKey("genome"))
    throw new Exception("params.global.genome is required.")

if(!isAppendOnlyMode && !params.global.genome.containsKey("assembly"))
    throw new Exception("params.global.genome.assembly is required. Choose of the profiles: " + ALLOWED_GENOME_ASSEMBLIES.join(', '))

if(!isAppendOnlyMode && params.global.genome.assembly == '')
    throw new Exception("params.global.genome.assembly cannot be empty. Choose of the profiles: " + ALLOWED_GENOME_ASSEMBLIES.join(', '))

if(!isAppendOnlyMode && !(params.global.genome.assembly in ALLOWED_GENOME_ASSEMBLIES))
    throw new Exception("The given genome assembly "+ params.global.genome.assembly + " is not implemented. Choose of the profiles: " + ALLOWED_GENOME_ASSEMBLIES.join(', '))

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces

include './../channels/file.nf' params(params)

include ARBORETO_WITH_MULTIPROCESSING                             from './processes/arboreto_with_multiprocessing'  params(params)
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

workflow scenic {

    take:
        // Expects (sampleId, loom)
        filteredLoom

    main:
        /* GRN */
        tfs = file(params.sc.scenic.grn.tfs)
        grn = ARBORETO_WITH_MULTIPROCESSING( filteredLoom.combine(runs), tfs )

        /* cisTarget motif analysis */
        // channel for SCENIC databases resources:
        motifsDb = Channel
            .fromPath( params.sc.scenic.cistarget.motifsDb )
            .collect() // use all files together in the ctx command
        motifsAnnotation = file(params.sc.scenic.cistarget.motifsAnnotation)
        ctx_mtf = CISTARGET__MOTIF( grn, motifsDb, motifsAnnotation, 'mtf' )

        /* cisTarget track analysis */
        if(params.sc.scenic.cistarget.tracksDb) {
            tracksDb = Channel
                .fromPath( params.sc.scenic.cistarget.tracksDb )
                .collect() // use all files together in the ctx command
            tracksAnnotation = file(params.sc.scenic.cistarget.tracksAnnotation)
            ctx_trk = CISTARGET__TRACK( grn, tracksDb, tracksAnnotation, 'trk' )
        }

        /* AUCell, motif regulons */
        auc_mtf = AUCELL__MOTIF( ctx_mtf, 'mtf' )

        if(params.sc.scenic.cistarget.tracksDb) {
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
            if(params.sc.scenic.cistarget.tracksDb) {
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
            if(params.sc.scenic.cistarget.tracksDb) {
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


workflow scenic_append {

    take:
        filteredLoom
        scopeLoom

    main:
        if(params.sc.scenic.containsKey("existingScenicLoom")) {
            scenicLoom = getChannelFromFilePath(
                params.sc.scenic.existingScenicLoom,
                params.sc.scenic.sampleSuffixWithExtension
            ).view {
            """
---------------------------------------------------------------------------
\u001B[32m Existing SCENIC loom detected \u001B[0m
\u001B[32m SCENIC won't run and this loom will be used as input to APPEND_SCENIC_LOOM \u001B[0m
---------------------------------------------------------------------------
            """
            }
        } else {
            scenicLoom = scenic( filteredLoom ).out
        }
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
        if(!("filteredLoom" in params.sc.scenic))
            throw new Exception("The given filteredLoom required parameter does not exist in the params.sc.scenic scope.")
        scenic( Channel.of( tuple("foobar", file(params.sc.scenic.filteredLoom)) ) )

}

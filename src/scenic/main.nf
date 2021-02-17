nextflow.enable.dsl=2

include {
    resolveParams;
} from './../utils/processes/config'

resolveParams(params, true)

def isAppendOnlyMode = params.tools.scenic.containsKey("existingScenicLoom")
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

include {
    getChannelFromFilePath;
} from './../channels/file.nf' params(params)

include {
    ARBORETO_WITH_MULTIPROCESSING;
} from './processes/arboreto_with_multiprocessing' params(params)
include {
    ADD_PEARSON_CORRELATION;
} from './processes/add_correlation' params(params)
include {
    CISTARGET as CISTARGET__MOTIF;
    CISTARGET as CISTARGET__TRACK;
} from './processes/cistarget' params(params)
include {
    AUCELL as AUCELL__MOTIF;
    AUCELL as AUCELL__TRACK;
} from './processes/aucell' params(params)
include {
    AGGREGATE_MULTI_RUNS_TO_LOOM as MULTI_RUNS_TO_LOOM__MOTIF;
    AGGREGATE_MULTI_RUNS_TO_LOOM as MULTI_RUNS_TO_LOOM__TRACK;
} from './workflows/aggregateMultiRuns' params(params)
include {
    PUBLISH_LOOM;
    MERGE_MOTIF_TRACK_LOOMS;
    APPEND_SCENIC_LOOM;
    VISUALIZE;
} from './processes/loomHandler' params(params)

// reporting:
include {
    GENERATE_REPORT;
    REPORT_TO_HTML;
} from './processes/reports.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

/*
 * SCENIC workflow
 */ 

// Create channel for the different runs
if(params.tools.scenic.containsKey("numRuns")) {
    runs = Channel.from( 1..params.tools.scenic.numRuns )
} else {
    runs = Channel.from( 1..1 )
}

workflow scenic {

    take:
        // Expects (sampleId, loom)
        filteredLoom

    main:
        /* GRN */
        tfs = file(params.tools.scenic.grn.tfs)
        grn = ARBORETO_WITH_MULTIPROCESSING( filteredLoom.combine(runs), tfs )
        grn_with_correlation = ADD_PEARSON_CORRELATION(grn)

        /* cisTarget motif analysis */
        // channel for SCENIC databases resources:
        motifsDb = Channel
            .fromPath( params.tools.scenic.cistarget.motifsDb )
            .collect() // use all files together in the ctx command
        motifsAnnotation = file(params.tools.scenic.cistarget.motifsAnnotation)
        ctx_mtf = CISTARGET__MOTIF( grn_with_correlation, motifsDb, motifsAnnotation, 'mtf' )

        /* cisTarget track analysis */
        if(params.tools.scenic.cistarget.tracksDb) {
            tracksDb = Channel
                .fromPath( params.tools.scenic.cistarget.tracksDb )
                .collect() // use all files together in the ctx command
            tracksAnnotation = file(params.tools.scenic.cistarget.tracksAnnotation)
            ctx_trk = CISTARGET__TRACK( grn_with_correlation, tracksDb, tracksAnnotation, 'trk' )
        }

        /* AUCell, motif regulons */
        auc_mtf = AUCELL__MOTIF( ctx_mtf, 'mtf' )

        if(params.tools.scenic.cistarget.tracksDb) {
            /* AUCell, track regulons */
            auc_trk = AUCELL__TRACK( ctx_trk, 'trk' )
        }

        // multi-runs aggregation:
        if(params.tools.scenic.containsKey("numRuns") && params.tools.scenic.numRuns > 1) {
            scenic_loom_mtf = MULTI_RUNS_TO_LOOM__MOTIF(
                filteredLoom,
                ctx_mtf,
                auc_mtf,
                'mtf'
            )
            if(params.tools.scenic.cistarget.tracksDb) {
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
            if(params.tools.scenic.cistarget.tracksDb) {
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
        if(params.tools.scenic.containsKey("existingScenicLoom")) {
            scenicLoom = getChannelFromFilePath(
                params.tools.scenic.existingScenicLoom,
                params.tools.scenic.sampleSuffixWithExtension
            )
            if(!params.containsKey('quiet')) {
                Channel.from('').view {
            """
---------------------------------------------------------------------------
\u001B[32m Existing SCENIC loom detected \u001B[0m
\u001B[32m SCENIC won't run and this loom will be used as input to APPEND_SCENIC_LOOM \u001B[0m
---------------------------------------------------------------------------
            """
                }
            }
        } else {
            scenicLoom = scenic( filteredLoom ).out
        }
        APPEND_SCENIC_LOOM(
            scopeLoom.map {
                // Extract only sampleId, path
                it -> tuple(it[0], it[1])
            }.join(
                scenicLoom.map {
                    // Extract only sampleId, path
                    it -> tuple(it[0], it[1])
                }
            ).ifEmpty{
                throw new Exception("Cannot append SCENIC loom to SCope loom because the IDs do not match.")
            }
        )
        if(!params.tools.scenic.skipReports) {
            report_notebook = GENERATE_REPORT(
                file(workflow.projectDir + params.tools.scenic.report_ipynb),
                APPEND_SCENIC_LOOM.out,
                "SCENIC_report"
            )
            REPORT_TO_HTML(report_notebook)
        }

    emit:
        APPEND_SCENIC_LOOM.out

}


// Uncomment to test
workflow {

    main:
        if(!("filteredLoom" in params.tools.scenic))
            throw new Exception("The given filteredLoom required parameter does not exist in the params.tools.scenic scope.")
        scenic( Channel.of( tuple(params.global.project_name, file(params.tools.scenic.filteredLoom)) ) )

}

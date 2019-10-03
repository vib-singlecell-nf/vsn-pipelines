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
include SC__SCENIC__MERGESCENICLOOMS                            from './processes/mergeScenicLooms'      params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

/*
 * SCENIC workflow
 */ 

// Create channel for the different runs
runs = Channel.from( 1..params.sc.scenic.numRuns )

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

        //visualize and merge
        SC__SCENIC__MERGESCENICLOOMS( auc_mtf, auc_trk )

}

// // Uncomment to test
// workflow {
//     main:
//         SCENIC( file( params.sc.scenic.filteredloom ) )
// }


//
// Tests
//  

// Test 1: SC__SCENIC__GRNBOOST2WITHOUTDASK (from processes/)
// Command: 
//  nextflow -C conf/test.config,scenic.config run main.test.nf -profile singularity --test SC__SCENIC__GRNBOOST2WITHOUTDASK
//

nextflow.preview.dsl=2

///////////////////////////////////////////
//  Define the parameters for all processes

include SC__SCENIC__GRNBOOST2WITHOUTDASK from './processes/grnboost2withoutDask' params(params)

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


workflow {
    main:
        switch(params.test) {
            case "SC__SCENIC__GRNBOOST2WITHOUTDASK":
                test_SC__SCENIC__GRNBOOST2WITHOUTDASK( file( params.sc.scenic.filteredloom ) )
            break;
            default:
                throw new Exception("The test parameters should be specified.")
            break;
        }
}
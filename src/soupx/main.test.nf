nextflow.enable.dsl=2

include { 
    INIT;
} from './../utils/workflows/utils' params(params)

INIT(params)

///////////////////////////////////////////
// How to run ?
// nextflow -C nextflow.config run [path-to-root-vsn]/src/soupx/main.test.nf --test <TEST>
// <TEST>:
// - SOUPX_CORRECT

///////////////////////////////////////////
// Define the parameters for all processes

include {
    SC__FILE_CONVERTER;
} from './../utils/processes/utils' params(params)

// Uncomment to test
include {
    getDataChannel;
} from './../channels/channels.nf' params(params)

workflow {

    main:
        switch(params.test) {
            case "SOUPX_CORRECT":
                include {
                    SOUPX_CORRECT;
                } from "./workflows/soupX"
                getDataChannel \
                    | SOUPX_CORRECT
            break;
            default:
                throw new Exception("The test parameters should be specified.")
            break;
        }

}

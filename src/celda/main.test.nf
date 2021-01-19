nextflow.enable.dsl=2

include { 
    INIT;
} from './../utils/workflows/utils' params(params)

INIT(params)

///////////////////////////////////////////
// How to run ?
// nextflow -C nextflow.config run [path-to-root-vsn]/src/celda/main.test.nf --test <TEST>
// <TEST>:
// - DECONTX_FILTER

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
            case "DECONTX_FILTER":
                include {
                    DECONTX_FILTER;
                } from "./workflows/decontX"
                getDataChannel \
                    | SC__FILE_CONVERTER \
                    | DECONTX_FILTER
            break;
            case "DECONTX_CORRECT":
                include {
                    DECONTX_CORRECT;
                } from "./workflows/decontX"
                getDataChannel \
                    | SC__FILE_CONVERTER \
                    | DECONTX_CORRECT
            break;
            default:
                throw new Exception("The test parameters should be specified.")
            break;
        }

}

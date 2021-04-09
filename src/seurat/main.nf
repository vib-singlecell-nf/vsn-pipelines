nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include {
    getDataChannel
} from '../channels/channels' params(params)
include {
    SC__FILE_CONVERTER;
} from '../utils/processes/utils.nf' params(params)
include {
    single_sample as SINGLE_SAMPLE
} from './workflows/single_sample.nf' params(params)

//////////////////////////////////////////////////////
// Define the workflow

workflow single_sample {
    
    main:
        getDataChannel | \
            SC__FILE_CONVERTER | \
            SINGLE_SAMPLE
}
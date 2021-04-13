nextflow.enable.dsl=2

import static groovy.json.JsonOutput.*



//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include { 
    INIT;
    getDataChannel;
} from '../utils/workflows/utils' params(params)
INIT(params)
include {
    SC__FILE_CONVERTER
} from '../utils/processes/utils' params(params)
include {
    getDataChannel
} from '../channels/channels' params(params)
include {
    QC_FILTER;
} from './workflows/qc_filter.nf' params(params)
include {
    SINGLE_SAMPLE
} from './workflows/single_sample.nf' params(params)


workflow single_sample {

    main:
        // run the pipeline
        getDataChannel | \
            SC__FILE_CONVERTER | \
            SINGLE_SAMPLE

}

workflow single_sample_qc {

    main:
        getDataChannel | \
            SC__FILE_CONVERTER | \
            QC_FILTER

}

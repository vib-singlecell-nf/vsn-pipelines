nextflow.preview.dsl=2

import static groovy.json.JsonOutput.*

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include '../utils/workflows/utils.nf' params(params)
INIT()
include '../src/utils/processes/utils.nf' params(params)
include SINGLE_SAMPLE from './workflows/single_sample.nf' params(params)


workflow single_sample {

    main:
        // run the pipeline
        getDataChannel | \
            SC__FILE_CONVERTER | \
            SINGLE_SAMPLE

}

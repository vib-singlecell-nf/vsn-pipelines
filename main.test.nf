nextflow.preview.dsl=2

///////////////////////////////////////////
//  Define the parameters for all processes

include './main.nf' params(params)
include SC__FILE_CONVERTER from '../utils/processes/utils.nf' params(params)
include getDataChannel from '../channels/channels.nf' params(params)

// Make the test workflow 
workflow test_single_sample {

    take:
        data

    main:
        single_sample( data )

}

workflow test_single_sample_scrublet {

    take:
        data

    main:
        single_sample_scrublet( data )

}

workflow {

    main:
        switch(params.test) {
            case "single_sample":
                getDataChannel | SC__FILE_CONVERTER | test_single_sample
            break;
            case "single_sample_scrublet":
                getDataChannel | SC__FILE_CONVERTER | test_single_sample_scrublet
            break;
            default:
                throw new Exception("The test parameters should be specified.")
            break;
        }

}

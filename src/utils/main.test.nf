//
// Tests
//  All the tests have to be run from the root of the project

// Test 1: SC__FILE_CONVERTER (from processes/)
// Command: 
//  nextflow -C nextflow_tiny.test.config,src/utils/conf/test.config run src/utils/main.test.nf --test SC__FILE_CONVERTER
//

nextflow.preview.dsl=2

///////////////////////////////////////////
//  Define the parameters for all processes

include SC__FILE_CONVERTER from './processes/utils' params(params.sc.file_converter + params.global + params)

// Uncomment to test
include getChannel as getTenXChannel from '../channels/tenx.nf'

// Make the test workflow 
workflow test_SC__FILE_CONVERTER {
    get:
        data
    main:
        SC__FILE_CONVERTER( data )
    emit:
        SC__FILE_CONVERTER.out
}

workflow test_SC__FILE_CONCATENATOR {
    get:
        data
    main:
        SC__FILE_CONCATENATOR( ( data ).collect() )
    emit:
        SC__FILE_CONCATENATOR.out
}


workflow {
    main:
        switch(params.test) {
            case "SC__FILE_CONVERTER":
                test_SC__FILE_CONVERTER( getTenXChannel( params.global.tenx_folder ) )
            break;
            case "SC__FILE_CONCATENATOR":
                test_SC__FILE_CONCATENATOR( getTenXChannel( params.global.tenx_folder ) )
            break;
            default:
                throw new Exception("The test parameters should be specified.")
            break;
        }
}
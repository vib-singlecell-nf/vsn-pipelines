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

include SC__FILE_CONVERTER from './processes/utils' params(params)

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

include FILTER_BY_CELL_METADATA from './workflows/filterByCellMetadata' params(params)
include SC__ANNOTATE_BY_CELL_METADATA from './processes/h5adAnnotate' params(params)

workflow {
    main:
        switch(params.test) {
            case "SC__FILE_CONVERTER":
                test_SC__FILE_CONVERTER( getTenXChannel( params.global.tenx_folder ) )
            break;
            case "SC__FILE_CONCATENATOR":
                test_SC__FILE_CONCATENATOR( getTenXChannel( params.global.tenx_folder ) )
            break;
            case "FILTER_BY_CELL_METADATA":
                if(params.sc.cell_filter) {
                    data = getTenXChannel( params.global.tenx_folder )
                    SC__FILE_CONVERTER( data )    
                    FILTER_BY_CELL_METADATA( SC__FILE_CONVERTER.out )
                }
            break;
            case "FILTER_AND_ANNOTATE_BY_CELL_METADATA":
                if(params.sc.cell_filter && params.sc.cell_annotate) {
                    data = getTenXChannel( params.global.tenx_folder )
                    SC__FILE_CONVERTER( data )    
                    FILTER_BY_CELL_METADATA( SC__FILE_CONVERTER.out )
                    SC__ANNOTATE_BY_CELL_METADATA( FILTER_BY_CELL_METADATA.out )
                }
            break;
            default:
                throw new Exception("The test parameters should be specified.")
            break;
        }
}
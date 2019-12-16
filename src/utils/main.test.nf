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

workflow {

    main:
        switch(params.test) {
            case "SC__FILE_CONVERTER":
                include SC__FILE_CONVERTER from './processes/utils' params(params)
                test_SC__FILE_CONVERTER( getTenXChannel( params.global.tenx_folder ) )
            break;
            case "SC__FILE_CONCATENATOR":
                test_SC__FILE_CONCATENATOR( getTenXChannel( params.global.tenx_folder ) )
            break;
            case "FILTER_BY_CELL_METADATA":
                // Imports
                include FILTER_BY_CELL_METADATA from './workflows/filterByCellMetadata' params(params)
                // Run 
                if(params.sc.cell_filter) {
                    data = getTenXChannel( params.global.tenx_folder )
                    SC__FILE_CONVERTER( data )    
                    FILTER_BY_CELL_METADATA( SC__FILE_CONVERTER.out )
                }
            break;
            case "FILTER_AND_ANNOTATE_BY_CELL_METADATA":
                // Imports
                include FILTER_BY_CELL_METADATA from './workflows/filterByCellMetadata' params(params)
                include SC__ANNOTATE_BY_CELL_METADATA from './processes/h5adAnnotate' params(params)
                // Run 
                if(params.sc.cell_filter && params.sc.cell_annotate) {
                    data = getTenXChannel( params.global.tenx_folder )
                    SC__FILE_CONVERTER( data )
                    FILTER_BY_CELL_METADATA( SC__FILE_CONVERTER.out )
                    SC__ANNOTATE_BY_CELL_METADATA( FILTER_BY_CELL_METADATA.out )
                }
            break;
            case "GET_METADATA_FROM_SRA":
                // Imports
                include getChannel as getSRAChannel from './../channels/sra' params(params)
                include SRA_TO_METADATA from './processes/sra' params(params)
                // Run
                sra = getSRAChannel( params.data.sra )
                db = file(params.utils.sra_metadata.sraDbOutDir + "/SRAmetadb.sqlite")
                SRA_TO_METADATA( sra, db )
            break;
            case "DOWNLOAD_FROM_SRA":
                // Imports
                include DOWNLOAD_FROM_SRA from './workflows/downloadFromSRA' params(params)
                include SC__CELLRANGER__PREPARE_FOLDER from './../cellranger/processes/utils.nf'
                include SC__CELLRANGER__COUNT   from './../cellranger/processes/count'    params(params)
                // Run 
                DOWNLOAD_FROM_SRA(
                    tuple('SRP125768', ["w1118_15d_*"]) //["DGRP-551_*d_r*","w1118_*d_r*"]
                )
            break;
            case "DOWNLOAD_FROM_SRA_AND_RUN_CELLRANGER":
                // Imports
                include DOWNLOAD_FROM_SRA from './workflows/downloadFromSRA' params(params)
                include SC__CELLRANGER__PREPARE_FOLDER from './../cellranger/processes/utils.nf'
                include SC__CELLRANGER__COUNT   from './../cellranger/processes/count'    params(params)
                // Run 
                DOWNLOAD_FROM_SRA(
                    tuple('SRP125768', ["w1118_15d_r1"]) //["DGRP-551_*d_r*","w1118_*d_r*"] //,"w1118_30d_*" //"w1118_15d_*" 
                )
                SC__CELLRANGER__PREPARE_FOLDER( DOWNLOAD_FROM_SRA.out.groupTuple() ).view()
                SC__CELLRANGER__COUNT( file("/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/resources/refs/flybase/r6.16/cellranger/2.0.0/flybase_r6.16"), SC__CELLRANGER__PREPARE_FOLDER.out )
            break;
            default:
                throw new Exception("The test parameters should be specified.")
            break;
        }

}

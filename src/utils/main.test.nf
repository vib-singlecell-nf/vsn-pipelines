nextflow.enable.dsl=2

include { 
    INIT;
} from './workflows/utils' params(params)

INIT(params)

///////////////////////////////////////////
// How to run ?
// nextflow -C nextflow.config run [path-to-root-vsn]/src/utils/main.test.nf --test TEST
// TEST:
// - SC__FILE_CONVERTER
// - SC__FILE_CONCATENATOR
// - SC__ANNOTATE_BY_SAMPLE_METADATA
// - ANNOTATE_BY_CELL_METADATA
// - FILTER_BY_CELL_METADATA
// - FILTER_AND_ANNOTATE_BY_CELL_METADATA
// - GET_METADATA_FROM_SRA
// - GET_METADATA_FROM_SRA_WEB
// - DOWNLOAD_FROM_SRA_AND_RUN_CELLRANGER
// - UPDATE_FEATURE_NOMENCLATURE

///////////////////////////////////////////
// Define the parameters for all processes

include {
    SC__FILE_CONVERTER;
} from './processes/utils' params(params)

// Uncomment to test
include {
    getDataChannel;
} from '../channels/channels.nf' params(params)

// Make the test workflow 
workflow test_SC__FILE_CONVERTER {

    take:
        data

    main:
        SC__FILE_CONVERTER( data )

    emit:
        SC__FILE_CONVERTER.out

}

workflow test_SC__FILE_CONCATENATOR {

    take:
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
                getDataChannel | test_SC__FILE_CONVERTER
            break;
            case "SC__FILE_CONCATENATOR":
                getDataChannel | test_SC__FILE_CONCATENATOR
            break;
            case "SC__ANNOTATE_BY_SAMPLE_METADATA":
                // Imports
                include {
                    SC__ANNOTATE_BY_SAMPLE_METADATA;
                } from './processes/h5adAnnotate' params(params)
                // Run
                if(params.hasUtilsParams("sample_annotate")) {
                    getDataChannel | \
                        SC__FILE_CONVERTER | \
                        SC__ANNOTATE_BY_SAMPLE_METADATA
                } else {
                    throw new Exception("Missing sample_annotate params in config.")
                }
                break;
            case "ANNOTATE_BY_CELL_METADATA":
                // Imports
                include {
                    STATIC__ANNOTATE_BY_CELL_METADATA as ANNOTATE_BY_CELL_METADATA;
                } from './workflows/annotateByCellMetadata' params(params)
                // Run
                if(params.hasUtilsParams("cell_annotate")) {
                    getDataChannel | \
                        SC__FILE_CONVERTER
                    ANNOTATE_BY_CELL_METADATA( 
                        SC__FILE_CONVERTER.out, 
                        null
                    )
                } else {
                    throw new Exception("Missing cell_annotate params in config.")
                }
            break;
            case "FILTER_BY_CELL_METADATA":
                // Imports
                include {
                    FILTER_BY_CELL_METADATA;
                } from './workflows/filterByCellMetadata' params(params)
                // Run 
                if(params.hasUtilsParams("cell_filter")) {
                    getDataChannel | \
                        SC__FILE_CONVERTER
                    FILTER_BY_CELL_METADATA(
                        SC__FILE_CONVERTER.out,
                        null
                    )
                }
            break;
            case "FILTER_AND_ANNOTATE_BY_CELL_METADATA":
                // Imports
                include {
                    STATIC__ANNOTATE_BY_CELL_METADATA as ANNOTATE_BY_CELL_METADATA;
                } from './workflows/annotateByCellMetadata' params(params)
                include {
                    FILTER_BY_CELL_METADATA;
                } from './workflows/filterByCellMetadata' params(params)
                // Run
                if(params.hasUtilsParams("cell_annotate")) {
                    getDataChannel | \
                        SC__FILE_CONVERTER

                    ANNOTATE_BY_CELL_METADATA(
                        SC__FILE_CONVERTER.out,
                        null
                    )
                    FILTER_BY_CELL_METADATA(
                        ANNOTATE_BY_CELL_METADATA.out,
                        null
                    )
                }
            break;
            case "GET_METADATA_FROM_SRA":
                // Imports
                include {
                    getChannel as getSRAChannel;
                } from './../channels/sra' params(params)
                include {
                    SRA_TO_METADATA;
                } from './processes/sra' params(params)
                // Run
                sra = getSRAChannel( params.data.sra )
                db = file(params.getUtilsParams("sra_metadata").sraDbOutDir + "/SRAmetadb.sqlite")
                SRA_TO_METADATA( sra, db )
            break;
            case "GET_METADATA_FROM_SRA_WEB":
                // Imports
                include {
                    getChannel as getSRAChannel;
                } from './../channels/sra' params(params)
                include {
                    SRA_TO_METADATA;
                } from './processes/sra' params(params)
                // Run
                sra = getSRAChannel( params.data.sra )
                SRA_TO_METADATA( sra, file('NO_FILE') )
            break;
            case "DOWNLOAD_FROM_SRA":
                // Imports
                include {
                    DOWNLOAD_FROM_SRA;
                } from './workflows/downloadFromSRA' params(params)
                include {
                    SC__CELLRANGER__PREPARE_FOLDER;
                } from './../cellranger/processes/utils.nf'
                include {
                    SC__CELLRANGER__COUNT;
                } from './../cellranger/processes/count'    params(params)
                // Run 
                DOWNLOAD_FROM_SRA(
                    tuple('SRP162698', ["10x, sample 1", "10x, sample 2"])
                )
            break;
            case "DOWNLOAD_FROM_SRA_AND_RUN_CELLRANGER":
                // Imports
                include {
                    DOWNLOAD_FROM_SRA;
                } from './workflows/downloadFromSRA' params(params)
                include {
                    SC__CELLRANGER__PREPARE_FOLDER;
                } from './../cellranger/processes/utils.nf'
                include {
                    SC__CELLRANGER__COUNT;
                } from './../cellranger/processes/count'    params(params)
                // Run 
                DOWNLOAD_FROM_SRA(
                    tuple('SRP125768', ["w1118_15d_r1"]) //["DGRP-551_*d_r*","w1118_*d_r*"] //,"w1118_30d_*" //"w1118_15d_*" 
                )
                SC__CELLRANGER__PREPARE_FOLDER( DOWNLOAD_FROM_SRA.out.groupTuple() ).view()
                SC__CELLRANGER__COUNT( file("/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/resources/refs/flybase/r6.16/cellranger/2.0.0/flybase_r6.16"), SC__CELLRANGER__PREPARE_FOLDER.out )
            break;
            case "UPDATE_FEATURE_NOMENCLATURE":
                // Imports
                include {
                    SC__UTILS__EXTRACT_FEATURE_METADATA;
                } from './processes/h5adExtractMetadata' params(params)
                include {
                    FLYBASER__CONVERT_FBGN_TO_GENE_SYMBOL;
                } from './../flybaser/processes/convertNomenclature' params(params)
                include {
                    SC__UTILS__UPDATE_FEATURE_METADATA_INDEX;
                } from './processes/h5adUpdateMetadata' params(params)
                // Run
                test = Channel.of(tuple('sample-id', 'path-to-h5ad'))
                SC__UTILS__EXTRACT_FEATURE_METADATA( test )
                FLYBASER__CONVERT_FBGN_TO_GENE_SYMBOL( SC__UTILS__EXTRACT_FEATURE_METADATA.out )
                SC__UTILS__UPDATE_FEATURE_METADATA_INDEX( test.join(FLYBASER__CONVERT_FBGN_TO_GENE_SYMBOL.out) )
            break;
            case "SC__H5AD_BEAUTIFY":
                // Imports
                include {
                    SC__H5AD_BEAUTIFY;
                } from './processes/h5adUpdate' params(params)
                // Run
                if(params.hasUtilsParams("file_cleaner")) {
                    getDataChannel | \
                    map {
                        it -> tuple(it[0], it[1], null)
                    } | \
                        SC__H5AD_BEAUTIFY 
                } else {
                    throw new Exception("Missing file_cleaner params in config.")
                }
                break;
            default:
                throw new Exception("The test parameters should be specified.")
            break;
        }

}

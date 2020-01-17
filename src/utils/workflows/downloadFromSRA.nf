//
// Version: 
// Test: 
// Command: 
//
/*
 * QC workflow 
 * Source:
 * 
 */ 

nextflow.preview.dsl=2

import java.nio.file.Files
import java.nio.file.Paths

//////////////////////////////////////////////////////
//  process imports:

include GET_SRA_DB from './../processes/sra' params(params)
include SRA_TO_METADATA from './../processes/sra' params(params)
include DOWNLOAD_FASTQS_FROM_SRA_ACC_ID from './../../sratoolkit/processes/downloadFastQ' params(params)
include NORMALIZE_SRA_FASTQS from './../processes/sra' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

// dataParams = params.data.sra
utilsParams = params.utils

if(!utilsParams.containsKey("sra_metadata"))
    throw new Exception("DOWNLOAD_FROM_SRA workflow requires sra_metadata.config")

workflowParams = params.utils.sra_metadata

workflow DOWNLOAD_FROM_SRA {

    take:
        // Expects (sraProjectId, sampleFilters)
        sra

    main:
        db = file('NO_FILE')
        sraDbFile = workflowParams.sraDb != '' ? file(workflowParams.sraDb): file(workflowParams.sraDbOutDir + "/SRAmetadb.sqlite")
        if(sraDbFile.exists() 
            && sraDbFile.canRead()
            && !workflowParams.sraDbForceDownload) {
            println("Local SRA database detected ${sraDbFile}!")
            db = sraDbFile
        } else {
            if(workflowParams.sraDbForceDownload
                || workflowParams.sraDb == '') {
                println("Downloading SRA database to ${sraDbFile}...")
                db = GET_SRA_DB()
                println("Done!")
            }
        }
        // Get metadata for the given SRA Project ID and keep only the samples that passes the given sampleFilters
        metadata = SRA_TO_METADATA( 
            sra,
            db
        ).splitCsv(
            header:true,
            sep: '\t'
        ).map {
            row -> tuple(row.run_accession, row.sample_name)
        }
        // Download and compress all the SRA runs defined in the metadata
        data = DOWNLOAD_FASTQS_FROM_SRA_ACC_ID( 
            metadata 
        ).join(
            metadata
        ).map {
            // Put sample as primary key
            run -> tuple(run[2], run[1])
        }
        out = NORMALIZE_SRA_FASTQS( data )

    emit:
        out

}

// workflow test {
//     Channel
//         .fromFilePairs('work/**/SRR*_{1,2}.fastq.gz')
// }

nextflow.enable.dsl=2

//////////////////////////////////////////////////////
// process imports:
include { SCTK__BARCODE_CORRECTION; } from './processes/barcode_correction.nf'
include { SCTK__BARCODE_10X_SCATAC_FASTQ; } from './processes/barcode_10x_scatac_fastqs.nf'
include { SCTK__EXTRACT_AND_CORRECT_BIORAD_BARCODE; } from './processes/extract_and_correct_biorad_barcode.nf'
include { BAP__BIORAD_DEBARCODE; } from './../bap/workflows/bap_debarcode.nf'

include {
    SIMPLE_PUBLISH as PUBLISH_BC_STATS;
    SIMPLE_PUBLISH as PUBLISH_BR_BC_STATS;
} from '../../src/utils/processes/utils.nf'

//////////////////////////////////////////////////////
//  Define the workflow


/* Barcode correction */
workflow barcode_correction {
    take:
        data

    main:

        // gather barcode whitelists from params into a channel:
        wl = Channel.empty()
        wl_cnt = 0
        params.tools.singlecelltoolkit.barcode_correction.whitelist.each { k, v ->
            if(v != '') {
                wl = wl.mix( Channel.of(tuple(k, file(v)) ))
                wl_cnt = wl_cnt + 1
            }
        }

        /* TO DO: fix ability to skip barcode correction */
        if(wl_cnt == 0) {
            if(!params.containsKey('quiet')) {
                println("No whitelist files were found in 'params.tools.singlecelltoolkit.barcode_correction.whitelist'. Skipping barcode correction for standard-type samples.")
            }
            // run barcode demultiplexing on each read+barcode:
            fastq_dex = SCTK__BARCODE_10X_SCATAC_FASTQ(data)
        } else {
            // join wl to the data channel:
            data_wl = wl.cross( data.map { it -> tuple(it[1], it[0], it[2], it[3], it[4]) } ) // technology, sampleId, R1, R2, R3
                        .map { it -> tuple(it[1][1], it[1][0],           // sampleId, technology
                                           it[1][2], it[1][3], it[1][4], // R1, R2, R3
                                           it[0][1]                      // whitelist
                                           ) }

            // run barcode correction against a whitelist:
            fastq_bc_corrected = SCTK__BARCODE_CORRECTION(data_wl)
            PUBLISH_BC_STATS(fastq_bc_corrected.map { it -> tuple(it[0], it[5]) }, '.corrected.bc_stats.log', 'reports/barcode')

            // run barcode demultiplexing on each read+barcode:
            fastq_dex = SCTK__BARCODE_10X_SCATAC_FASTQ(
                fastq_bc_corrected.map { it -> tuple(*it[0..4]) }
            )
        }

    emit:
        fastq_dex
}


workflow biorad_bc {

    take:
        data_biorad

    main:

        /* run BioRad barcode correction and debarcoding separately: */
        // using BAP:
        //fastq_dex_br = BAP__BIORAD_DEBARCODE(data.biorad.map{ it -> tuple(it[0], it[2], it[4]) })

        // using singlecelltoolkit:
        fastq_dex_br = SCTK__EXTRACT_AND_CORRECT_BIORAD_BARCODE(data_biorad.map{ it -> tuple(it[0], it[1], it[2], it[4]) })
        PUBLISH_BR_BC_STATS(fastq_dex_br.map { it -> tuple(it[0], it[3]) }, '.corrected.bc_stats.log', 'reports/barcode')

    emit:
        fastq_dex_br.map { it -> tuple(*it[0..2]) }

}


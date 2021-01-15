nextflow.enable.dsl=2

//////////////////////////////////////////////////////
// process imports:
include { SC__SINGLECELLTOOLKIT__BARCODE_CORRECTION; } from './../../src/singlecelltoolkit/processes/barcode_correction.nf' params(params)
include { SC__SINGLECELLTOOLKIT__DEBARCODE_10X_FASTQ; } from './../../src/singlecelltoolkit/processes/debarcode_10x_scatac_fastqs.nf' params(params)
include { SC__TRIMGALORE__TRIM; } from './../../src/trimgalore/processes/trim.nf' params(params)

// workflow imports:
include { BWA_MAPPING_PE; } from './../../src/bwamaptools/main.nf' params(params)
include { BAM_TO_FRAGMENTS; } from './../../src/sinto/main.nf' params(params)
include { BAP__BIORAD_DEBARCODE; } from './../../src/bap/workflows/bap_debarcode.nf' params(params)

include {
    PUBLISH as PUBLISH_BC_STATS;
    PUBLISH as PUBLISH_FASTQS_PE1;
    PUBLISH as PUBLISH_FASTQS_PE2;
    PUBLISH as PUBLISH_FASTQS_TRIMLOG_PE1;
    PUBLISH as PUBLISH_FASTQS_TRIMLOG_PE2;
} from "../../src/utils/workflows/utils.nf" params(params)


//////////////////////////////////////////////////////
//  Define the workflow 

workflow ATAC_PREPROCESS_WITH_METADATA {

    take:
        metadata

    main:
        // import metadata
        data = Channel.from(metadata)
                      .splitCsv(
                          header:true,
                          sep: '\t'
                          )
                      .map {
                          row -> tuple(
                              row.sample_name,
                              row.technology,
                              row.fastq_PE1_path,
                              row.fastq_barcode_path,
                              row.fastq_PE2_path
                              )
                      }
                      .branch {
                        biorad:   it[1] == 'biorad'
                        standard: true // capture all other technology types here
                      }

        // run biorad barcode correction and debarcoding separately:
        fastq_dex_br = BAP__BIORAD_DEBARCODE(data.biorad.map{ it -> tuple(it[0], it[2], it[4]) })

        /* Barcode correction */
        // gather barcode whitelists from params into a channel:
        wl = Channel.empty()
        params.tools.singlecelltoolkit.barcode_correction.whitelist.each { k, v ->
            if(v != '') {
                wl = wl.mix( Channel.of(tuple(k, file(v)) ))
            }
        }

        // join wl to the data channel:
        data_wl = wl.cross( data.standard.map { it -> tuple(it[1], it[0], it[2], it[3], it[4]) } ) // technology, sampleId, R1, R2, R3
                .map { it -> tuple(it[1][1], it[1][0],           // sampleId, technology
                                   it[1][2], it[1][3], it[1][4], // R1, R2, R3
                                   it[0][1]                      // whitelist
                                   ) }

        // run barcode correction against a whitelist:
        fastq_bc_corrected = SC__SINGLECELLTOOLKIT__BARCODE_CORRECTION(data_wl.map{ it -> tuple(it[0], it[3], it[5]) } )
        PUBLISH_BC_STATS(fastq_bc_corrected.map { it -> tuple(it[0], it[2]) }, 'corrected.bc_stats', 'log', 'fastq', false)


        // run barcode demultiplexing on each read+barcode:
        fastq_dex = SC__SINGLECELLTOOLKIT__DEBARCODE_10X_FASTQ(
            data.standard.join(fastq_bc_corrected).map { it -> tuple(it[0], it[2], it[5], it[4]) }
        )

        // concatenate the read channels:
        fastq_dex = fastq_dex.concat(fastq_dex_br)

        // run adapter trimming:
        fastq_dex_trim = SC__TRIMGALORE__TRIM(fastq_dex)
        // publish fastq output:
        PUBLISH_FASTQS_PE1(fastq_dex_trim, 'R1.fastq', 'gz', 'fastq', false)
        PUBLISH_FASTQS_PE2(fastq_dex_trim.map{ it -> tuple(it[0], it[2]) }, 'R2.fastq', 'gz', 'fastq', false)
        PUBLISH_FASTQS_TRIMLOG_PE1(fastq_dex_trim.map{ it -> tuple(it[0], it[3]) }, 'R1.trimming_report', 'txt', 'fastq', false)
        PUBLISH_FASTQS_TRIMLOG_PE2(fastq_dex_trim.map{ it -> tuple(it[0], it[4]) }, 'R2.trimming_report', 'txt', 'fastq', false)

        // map with bwa mem:
        bam = BWA_MAPPING_PE(fastq_dex_trim.map { it -> tuple(it[0..2]) })

        // generate a fragments file:
        fragments = BAM_TO_FRAGMENTS(bam)

    emit:
        bam
        fragments

}


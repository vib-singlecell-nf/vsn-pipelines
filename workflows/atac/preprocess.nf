nextflow.preview.dsl=2

//////////////////////////////////////////////////////
// process imports:
//include { SC__SNAPTOOLS__DEX_FASTQ; } from './../../src/snaptools/processes/dex_fastq.nf' params(params)
include { SC__SINGLECELLTOOLKIT__DEBARCODE_10X_FASTQ; } from './../../src/singlecelltoolkit/processes/debarcode_10x_scatac_fastqs.nf' params(params)
include { SC__TRIMGALORE__TRIM; } from './../../src/trimgalore/processes/trim.nf' params(params)

// workflow imports:
include { BWA_MAPPING_PE; } from './../../src/bwamaptools/main.nf' params(params)
include { BAM_TO_FRAGMENTS; } from './../../src/sinto/main.nf' params(params)
include { BAP__BIORAD_DEBARCODE; } from './../../src/bap/workflows/bap_debarcode.nf' params(params)


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
                        standard: it[1] == 'standard'
                        biorad:   it[1] == 'biorad'
                      }

        // run biorad debarcoding
        fastq_dex_br = BAP__BIORAD_DEBARCODE(data.biorad.map{ it -> tuple(it[0], it[2], it[4]) })

        // run barcode demultiplexing on each read+barcode:
        fastq_dex = SC__SINGLECELLTOOLKIT__DEBARCODE_10X_FASTQ(data.standard.map{ it -> tuple(it[0], it[2], it[3], it[4]) })

        // concatenate the read channels:
        fastq_dex = fastq_dex.concat(fastq_dex_br)

        // run adapter trimming:
        fastq_dex_trim = SC__TRIMGALORE__TRIM(fastq_dex)

        // map with bwa mem:
        bam = BWA_MAPPING_PE(fastq_dex_trim.map { it -> tuple(it[0..2]) })

        // generate a fragments file:
        fragments = BAM_TO_FRAGMENTS(bam)

        //out = bam.map { it -> tuple(it[0], [it[1], it[2]]) }
        //            .join(fragments.map { it -> tuple(it[0], [it[1], it[2]]) })

    emit:
        bam
        fragments

}


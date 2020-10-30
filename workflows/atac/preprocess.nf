nextflow.preview.dsl=2

//////////////////////////////////////////////////////
// process imports:
include { SC__SNAPTOOLS__DEX_FASTQ; } from './../../src/snaptools/processes/dex_fastq.nf' params(params)
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
                              row.fastq_type,
                              row.fastq_path,
                              row.fastq_path_barcode
                              )
                      }
                      .branch {
                        standard: it[1] == 'standard'
                        biorad:   it[1] == 'biorad'
                      }

        // run biorad debarcoding
        fastq_dex_br = BAP__BIORAD_DEBARCODE(
            data.biorad.map{ it -> tuple(it[0], it[2], it[3]) }
                       .groupTuple(sort: { it[0] } )
                       .map { it -> tuple(it[0], it[2][0], it[2][1]) }
        )

        // run barcode demultiplexing on each read+barcode:
        fastq_dex = SC__SNAPTOOLS__DEX_FASTQ(data.standard.map{ it -> tuple(it[0], it[2], it[3], it[4]) })

        // group parired fastqs by sampleID, keeping the R1, R2 order:
        fastq_dex_paired = fastq_dex.map { it -> tuple(it[0], [it[1], it[2]]) }
                     .groupTuple(sort: { it[0] } )
                     .map { it -> tuple(it[0], it[1][0][1], it[1][1][1]) }

        // concatenate the read channels:
        fastq_dex_paired = fastq_dex_paired.concat(fastq_dex_br)

        // run adapter trimming:
        fastq_dex_trim = SC__TRIMGALORE__TRIM(fastq_dex_paired)

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


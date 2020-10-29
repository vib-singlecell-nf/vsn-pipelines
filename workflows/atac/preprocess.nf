nextflow.preview.dsl=2

//////////////////////////////////////////////////////
// process imports:
include { SC__SNAPTOOLS__DEX_FASTQ; } from './../../src/snaptools/processes/dex_fastq.nf' params(params)
include { SC__TRIMGALORE__TRIM; } from './../../src/trimgalore/processes/trim.nf' params(params)

// workflow imports:
include { BWA_MAPPING_PE; } from './../../src/bwamaptools/main.nf' params(params)
include { BAM_TO_FRAGMENTS; } from './../../src/sinto/main.nf' params(params)


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
                              row.fastq_type,
                              row.fastq_path,
                              row.fastq_path_barcode
                              )
                      }

        // run barcode demultiplexing on each read+barcode:
        fastq_dex = SC__SNAPTOOLS__DEX_FASTQ(data)

        // group parired fastqs by sampleID, keeping the R1, R2 order:
        fastq_dex_paired = fastq_dex.map { it -> tuple(it[0], [it[1], it[2]]) }
                     .groupTuple(sort: { it[0] } )
                     .map { it -> tuple(it[0], it[1][0][1], it[1][1][1]) }

        // run adapter trimming:
        fastq_dex_trim = SC__TRIMGALORE__TRIM(fastq_dex_paired)

        // map with bwa mem:
        mapped = BWA_MAPPING_PE(fastq_dex_trim.map { it -> tuple(it[0..2]) })

        // generate a fragments file:
        fragments = BAM_TO_FRAGMENTS(mapped)

    emit:
        fragments

}


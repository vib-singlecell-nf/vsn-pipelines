nextflow.enable.dsl=2

//////////////////////////////////////////////////////
// process imports:
include { SCTK__BARCODE_CORRECTION; } from './../../src/singlecelltoolkit/processes/barcode_correction.nf' params(params)
include { SCTK__BARCODE_10X_SCATAC_FASTQ; } from './../../src/singlecelltoolkit/processes/barcode_10x_scatac_fastqs.nf' params(params)
include { SCTK__EXTRACT_AND_CORRECT_BIORAD_BARCODE; } from './../../src/singlecelltoolkit/processes/extract_and_correct_biorad_barcode.nf' params(params)
include { TRIMGALORE__TRIM; } from './../../src/trimgalore/processes/trim.nf' params(params)
include {
    FASTP__ADAPTER_TRIMMING as FASTP__TRIM;
} from './../../src/fastp/processes/adapter_trimming.nf' params(params)

// workflow imports:
include { BWA_MAPPING_PE; } from './../../src/bwamaptools/main.nf' params(params)
include { BAM_TO_FRAGMENTS; } from './../../src/sinto/main.nf' params(params)
include { BAP__BIORAD_DEBARCODE; } from './../../src/bap/workflows/bap_debarcode.nf' params(params)

include {
    SIMPLE_PUBLISH as PUBLISH_BC_STATS;
    SIMPLE_PUBLISH as PUBLISH_BR_BC_STATS;
    SIMPLE_PUBLISH as PUBLISH_FASTQS_TRIMLOG_PE1;
    SIMPLE_PUBLISH as PUBLISH_FASTQS_TRIMLOG_PE2;
    SIMPLE_PUBLISH as PUBLISH_FASTQS_TRIMLOG_FASTP;
    SIMPLE_PUBLISH as PUBLISH_FRAGMENTS;
    SIMPLE_PUBLISH as PUBLISH_FRAGMENTS_INDEX;
} from "../../src/utils/processes/utils.nf" params(params)


//////////////////////////////////////////////////////
//  Define the workflow 

workflow ATAC_PREPROCESS {

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


        /* Barcode correction */
        // gather barcode whitelists from params into a channel:
        wl = Channel.empty()
        wl_cnt = 0
        params.tools.singlecelltoolkit.barcode_correction.whitelist.each { k, v ->
            if(v != '') {
                wl = wl.mix( Channel.of(tuple(k, file(v)) ))
                wl_cnt = wl_cnt + 1
            }
        }

        if(wl_cnt == 0) {
            if(!params.containsKey('quiet')) {
                println("No whitelist files were found in 'params.tools.singlecelltoolkit.barcode_correction.whitelist'. Skipping barcode correction for standard-type samples.")
            }
            // run barcode demultiplexing on each read+barcode:
            fastq_dex = SCTK__BARCODE_10X_SCATAC_FASTQ(
                data.standard.map { it -> tuple(it[0], it[2], it[3], it[4]) }
            )
        } else {
            // join wl to the data channel:
            data_wl = wl.cross( data.standard.map { it -> tuple(it[1], it[0], it[2], it[3], it[4]) } ) // technology, sampleId, R1, R2, R3
                        .map { it -> tuple(it[1][1], it[1][0],           // sampleId, technology
                                           it[1][2], it[1][3], it[1][4], // R1, R2, R3
                                           it[0][1]                      // whitelist
                                           ) }

            // run barcode correction against a whitelist:
            fastq_bc_corrected = SCTK__BARCODE_CORRECTION(data_wl.map{ it -> tuple(it[0], it[3], it[5]) } )
            PUBLISH_BC_STATS(fastq_bc_corrected.map { it -> tuple(it[0], it[2]) }, '.corrected.bc_stats.log', 'fastq')


            // run barcode demultiplexing on each read+barcode:
            fastq_dex = SCTK__BARCODE_10X_SCATAC_FASTQ(
                data.standard.join(fastq_bc_corrected).map { it -> tuple(it[0], it[2], it[5], it[4]) }
            )

        }

        /* run BioRad barcode correction and debarcoding separately: */
        // using BAP:
        //fastq_dex_br = BAP__BIORAD_DEBARCODE(data.biorad.map{ it -> tuple(it[0], it[2], it[4]) })
        // using singlecelltoolkit:
        fastq_dex_br = SCTK__EXTRACT_AND_CORRECT_BIORAD_BARCODE(data.biorad.map{ it -> tuple(it[0], it[2], it[4]) })
        PUBLISH_BR_BC_STATS(fastq_dex_br.map { it -> tuple(it[0], it[3]) }, '.corrected.bc_stats.log', 'reports')


        // concatenate the read channels:
        fastq_dex = fastq_dex.concat(fastq_dex_br.map{ it -> tuple(it[0], it[1],it[2])})

        // run adapter trimming:
        switch(params.atac_preprocess_tools.adapter_trimming_method) {
            case 'Trim_Galore':
                fastq_dex_trim = TRIMGALORE__TRIM(fastq_dex);
                PUBLISH_FASTQS_TRIMLOG_PE1(fastq_dex_trim.map{ it -> tuple(it[0], it[3]) }, '.R1.trimming_report.txt', 'reports');
                PUBLISH_FASTQS_TRIMLOG_PE2(fastq_dex_trim.map{ it -> tuple(it[0], it[4]) }, '.R2.trimming_report.txt', 'reports');
                break;
            case 'fastp':
                fastq_dex_trim = FASTP__TRIM(fastq_dex);
                PUBLISH_FASTQS_TRIMLOG_FASTP(fastq_dex_trim.map{ it -> tuple(it[0], it[3]) }, '.fastp.trimming_report.html', 'reports');
                break;
        }

        // map with bwa mem:
        bam = BWA_MAPPING_PE(fastq_dex_trim.map { it -> tuple(it[0..2]) },
                             params.atac_preprocess_tools.mark_duplicates_method)

        // generate a fragments file:
        fragments = BAM_TO_FRAGMENTS(bam)
        // publish fragments output:
        PUBLISH_FRAGMENTS(fragments.map{ it -> tuple(it[0..1]) }, '.sinto.fragments.tsv.gz', 'fragments')
        PUBLISH_FRAGMENTS_INDEX(fragments.map{ it -> tuple(it[0],it[2]) }, '.sinto.fragments.tsv.gz.tbi', 'fragments')

    emit:
        bam
        fragments

}


nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:
include { 
    BAP__BIORAD_DEBARCODE as BIORAD_DEBARCODE;
    BAP__MERGE_FASTQS as MERGE_FASTQS;
} from './../processes/biorad_debarcode.nf' params(params)

include {
    PUBLISH as PUBLISH_BC_STATS_BR;
} from "../../utils/workflows/utils.nf" params(params)


//////////////////////////////////////////////////////
// Define the workflow

workflow BAP__BIORAD_DEBARCODE {

    take:
        data // a channel of [val(sampleId), path(fastq_PE1), path(fastq_PE2)]

    main:

        BIORAD_DEBARCODE(data) |
            MERGE_FASTQS

        PUBLISH_BC_STATS_BR(BIORAD_DEBARCODE.out.map { it -> tuple(it[0], it[3]) }, 'corrected.bc_stats', 'log', 'fastq', false)

    emit:
        MERGE_FASTQS.out

}


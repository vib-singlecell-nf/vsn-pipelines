nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:
include { 
    SC__BAP__BIORAD_DEBARCODE;
    SC__BAP__MERGE_FASTQS;
} from './../processes/biorad_debarcode.nf' params(params)



//////////////////////////////////////////////////////
// Define the workflow

workflow BAP__BIORAD_DEBARCODE {

    take:
        data // a channel of [val(sampleId), path(fastq_PE1), path(fastq_PE2)]

    main:

        SC__BAP__BIORAD_DEBARCODE(data) |
            SC__BAP__MERGE_FASTQS

    emit:
        SC__BAP__MERGE_FASTQS.out

}


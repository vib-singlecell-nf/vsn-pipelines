nextflow.enable.dsl=2

include { 
    BCL2FASTQ__DEMULTIPLEX 
} from "./processes/demultiplex" params(params)

workflow demultiplex {
    // include { 
    //     BCL2FASTQ__DEMULTIPLEX 
    // } from "./processes/demultiplex" params(params)

    take:
        data
        // data format: tuple(runPath, sampleSheet)

    main:
        BCL2FASTQ__DEMULTIPLEX(data)

}
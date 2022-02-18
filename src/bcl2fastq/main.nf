nextflow.enable.dsl=2

include { 
    BCL2FASTQ__DEMULTIPLEX 
} from "./processes/demultiplex" params(params)

workflow demultiplex {

    take:
        data

    main:
        BCL2FASTQ__DEMULTIPLEX(data)

}
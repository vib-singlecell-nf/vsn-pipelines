nextflow.enable.dsl=2

include { 
    BCL2FASTQ__DEMULTIPLEX 
} from "./processes/demultiplex" params(params)

include {
    PROCESS_SAMPLESHEET
} from "./processes/process_samplesheet" params(params)

workflow demultiplex {

    take:
        data

    main:
        BCL2FASTQ__DEMULTIPLEX(data)

    emit:
        fastqs = BCL2FASTQ__DEMULTIPLEX.out.fastqs
        stats = BCL2FASTQ__DEMULTIPLEX.out.stats
        sampleSheet = BCL2FASTQ__DEMULTIPLEX.out.sampleSheet
}

workflow process_samplesheet {

    take:
        data

    main:
        PROCESS_SAMPLESHEET(data)

    emit:
        processedSampleSheet = PROCESS_SAMPLESHEET.out.processedSampleSheet
}
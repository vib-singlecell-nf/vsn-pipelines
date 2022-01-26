nextflow.enable.dsl=2

workflow demultiplex {
    include { BCL2FASTQ__DEMULTIPLEX } from "../process/demultiplex.nf" params(params)

    BCL2FASTQ__DEMULTIPLEX

}
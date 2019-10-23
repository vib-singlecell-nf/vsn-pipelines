nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces

include FASTP__CLEAN_AND_FASTQC from '../src/fastp/processes/clean_and_fastqc.nf' params(params.fastp + params.fastp.clean_and_fastqc + params)

//////////////////////////////////////////////////////
// Define the input data

params.genome = '/ddn1/vol1/staging/leuven/res_00001/genomes/homo_sapiens/hg38_iGenomes/iGenomes_Raw/Sequence/WholeGenomeFasta/genome.fa'
params.annotation = '/ddn1/vol1/staging/leuven/res_00001/genomes/homo_sapiens/hg38_iGenomes/iGenomes_Raw/Annotation/Archives/archive-2015-08-14-08-18-15/Genes/genes.gtf'


//////////////////////////////////////////////////////
//  Define the workflow 

workflow nemesh {

    /*
    * Create a channel for input read files
    */
    Channel
        .fromFilePairs( params.reads, size: 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .set { data }

    data.subscribe { println it }
    // selectedBarcodesByCustom.subscribe { println it }

    FASTP__CLEAN_AND_FASTQC( data )
}
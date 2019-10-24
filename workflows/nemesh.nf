nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces

include FASTP__CLEAN_AND_FASTQC from '../src/fastp/processes/clean_and_fastqc.nf' params(params.fastp + params.fastp.clean_and_fastqc + params)
include PICARD__FASTQ_TO_BAM from '../src/picard/processes/fastq_to_bam.nf' params(params)
include DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLBARCODE from '../src/dropseqtools/processes/tag_bam_with_read_sequence_extended.nf' params(params.dropseqtools.tag_unaligned_bam_with_cellbarcode + params)
include DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLMOLECULAR from '../src/dropseqtools/processes/tag_bam_with_read_sequence_extended.nf' params(params.dropseqtools.tag_unaligned_bam_with_cellmolecular + params)
include DROP_SEQ_TOOLS__FILTER_UNALIGNED_TAGGED_BAM from '../src/dropseqtools/processes/filter_bam.nf' params(params.dropseqtools.filter_unaligned_tagged_bam + params)
include DROP_SEQ_TOOLS__TRIM_SMART_UNALIGNED_TAGGED_FILTERED_BAM from '../src/dropseqtools/processes/trim_starting_sequence.nf' params(params.dropseqtools.trim_smart_unaligned_tagged_filtered_bam + params)
include DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART from '../src/dropseqtools/processes/polya_trimmer.nf' params(params.dropseqtools.trim_polya_unaligned_tagged_trimmed_smart + params)
include PICARD__BAM_TO_FASTQ from '../src/picard/processes/sam_to_fastq.nf' from params(params)

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
    PICARD__FASTQ_TO_BAM( FASTP__CLEAN_AND_FASTQC.out.fastq )
    DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLBARCODE( PICARD__FASTQ_TO_BAM.out.bam )
    DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLMOLECULAR( DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLBARCODE.out.bam )
    DROP_SEQ_TOOLS__FILTER_UNALIGNED_TAGGED_BAM( DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLMOLECULAR.out.bam )
    DROP_SEQ_TOOLS__TRIM_SMART_UNALIGNED_TAGGED_FILTERED_BAM( DROP_SEQ_TOOLS__FILTER_UNALIGNED_TAGGED_BAM.out.bam )
    DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART( DROP_SEQ_TOOLS__TRIM_SMART_UNALIGNED_TAGGED_FILTERED_BAM.out.bam )
    PICARD__BAM_TO_FASTQ( DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART.out.bam )
}
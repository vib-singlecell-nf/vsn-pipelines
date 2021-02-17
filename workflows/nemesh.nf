nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces

include {
    FASTP__CLEAN_AND_FASTQC;
} from '../src/fastp/processes/clean_and_fastqc.nf' params(params)
include PICARD__FASTQ_TO_BAM from '../src/picard/processes/fastq_to_bam.nf' params(params)
include {
    SC__DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLBARCODE;
    SC__DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLMOLECULAR;
} from '../src/dropseqtools/processes/tag_bam_with_read_sequence_extended.nf' params(params)
include {
    SC__DROP_SEQ_TOOLS__FILTER_UNALIGNED_TAGGED_BAM;
} from '../src/dropseqtools/processes/filter_bam.nf' params(params)
include {
    SC__DROP_SEQ_TOOLS__TRIM_SMART_UNALIGNED_TAGGED_FILTERED_BAM;
} from '../src/dropseqtools/processes/trim_starting_sequence.nf' params(params)
include {
    SC__DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART;
} from '../src/dropseqtools/processes/polya_trimmer.nf' params(params)
include {
    PICARD__BAM_TO_FASTQ;
} from '../src/picard/processes/sam_to_fastq.nf' params(params)
include {
    GZIP;
} from '../src/dropseqtools/processes/gzip.nf' params(params)
include {
    SC__STAR__BUILD_INDEX;
} from '../src/star/processes/build_genome.nf' params(params)
include {
    SC__STAR__LOAD_GENOME;
} from '../src/star/processes/load_genome.nf' params(params)
include {
    SC__STAR__MAP_COUNT;
} from '../src/star/processes/map_count.nf' params(params)
include {
    PICARD__SORT_SAM;
} from '../src/picard/processes/sort_sam.nf' params(params)
include {
    PICARD__CREATE_SEQUENCE_DICTIONARY;
} from '../src/picard/processes/create_sequence_dictionary.nf' params(params)
include {
    PICARD__MERGE_BAM_ALIGNMENT;
} from '../src/picard/processes/merge_bam_alignment.nf' params(params)
include {
    FORMAT_GTF_IGENOMES;
} from '../src/utils/processes/gtf.nf' params(params)
include {
    SC__DROP_SEQ_TOOLS__CONVERT_TO_REFFLAT;
} from '../src/dropseqtools/processes/convert_to_ref_flat.nf' params(params)
include {
    SC__DROP_SEQ_TOOLS__TAG_READ_WITH_GENE_EXON;
} from '../src/dropseqtools/processes/tag_read_with_gene_exon.nf' params(params)
include {
    SC__DROP_SEQ_TOOLS__DETECT_REPAIR_BARCODE_SYNTHESIS_ERRORS;
} from '../src/dropseqtools/processes/detect_bead_synthesis_errors.nf' params(params)
include {
    SC__DROP_SEQ_TOOLS__BAM_TAG_HISTOGRAM;
} from '../src/dropseqtools/processes/bam_tag_histogram.nf' params(params)
include {
    SC__DROPLET_UTILS__BARCODE_SELECTION;
} from '../src/dropletutils/processes/barcode_selection.nf' params(params)
include {
    SC__DROP_SEQ_TOOLS__DIGITAL_EXPRESSION;
} from '../src/dropseqtools/processes/digital_expression.nf' params(params)

//////////////////////////////////////////////////////
// Define the input data

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

    // Check if custom selected barcodes file has been specified
    if (params.getToolParams("nemesh").custom_selected_barcodes) {
        Channel
            .fromPath(params.getToolParams("nemesh").custom_selected_barcodes)
            .map {
                path -> tuple(path.baseName.split('\\.')[0], params.getToolParams("nemesh").custom_selected_barcodes, path)
            }
            .set { selectedBarcodesByCustom }
        selectedBarcodesByCustom.subscribe { println it }
    }

    FASTP__CLEAN_AND_FASTQC( data )
    PICARD__FASTQ_TO_BAM( FASTP__CLEAN_AND_FASTQC.out.fastq, path( params.global.tmpDir ) )
    SC__DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLBARCODE( PICARD__FASTQ_TO_BAM.out.bam )
    SC__DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLMOLECULAR( SC__DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLBARCODE.out.bam )
    SC__DROP_SEQ_TOOLS__FILTER_UNALIGNED_TAGGED_BAM( SC__DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLMOLECULAR.out.bam )
    SC__DROP_SEQ_TOOLS__TRIM_SMART_UNALIGNED_TAGGED_FILTERED_BAM( SC__DROP_SEQ_TOOLS__FILTER_UNALIGNED_TAGGED_BAM.out.bam )
    SC__DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART( SC__DROP_SEQ_TOOLS__TRIM_SMART_UNALIGNED_TAGGED_FILTERED_BAM.out.bam )
    PICARD__BAM_TO_FASTQ( SC__DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART.out.bam, path( params.global.tmpDir ) )
    GZIP( PICARD__BAM_TO_FASTQ.out.fastq )
    SC__STAR__BUILD_INDEX( file(params.global.genome_annotation), path(params.global.genome) )
    // STAR_index = file("")
    // STAR__LOAD( STAR_index )
    SC__STAR__LOAD_GENOME( SC__STAR__BUILD_INDEX.out )
    SC__STAR__MAP_COUNT(
        SC__STAR__BUILD_INDEX.out,
        SC__STAR__LOAD_GENOME.out,
        GZIP.out.fastq_gz
    )
    PICARD__SORT_SAM( SC__STAR__MAP_COUNT.out.bam, path( params.global.tmpDir ) )
    PICARD__CREATE_SEQUENCE_DICTIONARY( file(params.global.genome), path( params.global.tmpDir ) )
    PICARD__MERGE_BAM_ALIGNMENT( 
        SC__DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART.out.bam,
        PICARD__SORT_SAM.out,
        file( params.global.genome ),
        PICARD__CREATE_SEQUENCE_DICTIONARY.out,
        file( params.global.tmpDir )
    )
    // FORMAT_GTF( file(params.global.genome_annotation) )
    FORMAT_GTF_IGENOMES( file(params.global.genome_annotation) )
    SC__DROP_SEQ_TOOLS__CONVERT_TO_REFFLAT(
        FORMAT_GTF_IGENOMES.out,
        PICARD__CREATE_SEQUENCE_DICTIONARY.out
    )
    SC__DROP_SEQ_TOOLS__TAG_READ_WITH_GENE_EXON(
        PICARD__MERGE_BAM_ALIGNMENT.out,
        SC__DROP_SEQ_TOOLS__CONVERT_TO_REFFLAT.out
    )
    SC__DROP_SEQ_TOOLS__DETECT_REPAIR_BARCODE_SYNTHESIS_ERRORS( SC__DROP_SEQ_TOOLS__TAG_READ_WITH_GENE_EXON.out )
    FINAL_BAM = SC__DROP_SEQ_TOOLS__DETECT_REPAIR_BARCODE_SYNTHESIS_ERRORS.out.bam
    SC__DROP_SEQ_TOOLS__BAM_TAG_HISTOGRAM( FINAL_BAM )
    SC__DROPLET_UTILS__BARCODE_SELECTION( SC__DROP_SEQ_TOOLS__BAM_TAG_HISTOGRAM.out )
    a = FINAL_BAM.combine(SC__DROPLET_UTILS__BARCODE_SELECTION.out.selectedCellBarcodesByKnee, by: 0)
    b = FINAL_BAM.combine(SC__DROPLET_UTILS__BARCODE_SELECTION.out.selectedCellBarcodesByInflection, by: 0)

    if (params.getToolParams("nemesh").custom_selected_barcodes) {
        c = FINAL_BAM.combine(selectedBarcodesByCustom, by: 0)
        SC__DROP_SEQ_TOOLS__DIGITAL_EXPRESSION(
            a.mix(b,c)
        )
    } else {
        SC__DROP_SEQ_TOOLS__DIGITAL_EXPRESSION(
            a.mix(b)
        )
    }

}
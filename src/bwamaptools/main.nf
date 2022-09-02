nextflow.enable.dsl=2

import java.nio.file.Paths

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include {
    BWAMAPTOOLS__BWA_MEM_PE as BWA_MEM_PE;
} from './processes/mapping.nf' params(params)
include {
    BWAMAPTOOLS__MAPPING_SUMMARY as MAPPING_SUMMARY;
} from './processes/mapping_summary.nf' params(params)
include {
    PICARD__MARK_DUPLICATES_AND_SORT;
} from './../../src/gatk/processes/mark_duplicates.nf' params(params)
include {
    PICARD__ESTIMATE_LIBRARY_COMPLEXITY;
} from './../../src/gatk/processes/estimate_library_complexity.nf' params(params)
include {
    GATK__MARK_DUPLICATES_SPARK;
} from './../../src/gatk/processes/mark_duplicates_spark.nf' params(params)
include {
    SIMPLE_PUBLISH as PUBLISH_BAM;
    SIMPLE_PUBLISH as PUBLISH_BAM_INDEX;
    SIMPLE_PUBLISH as PUBLISH_MAPPING_SUMMARY;
    SIMPLE_PUBLISH as PUBLISH_MARKDUPS_METRICS;
    SIMPLE_PUBLISH as PUBLISH_LIBRARY_METRICS;
} from "../utils/processes/utils.nf" params(params)

//////////////////////////////////////////////////////
// Define the workflow

workflow get_bwa_index {

    take:
        fasta_path

    main:

        bwa_fasta = Channel.fromPath(fasta_path)

        bwa_index_path = Paths.get(
                                   Paths.get(fasta_path).getParent().toString(),
                                   "*.{amb,ann,bwt,fai,flat,gdx,pac,sa,0123,bwt.2bit.64}"
                                   )
        bwa_index = Channel.fromPath(bwa_index_path,
                                     glob: true,
                                     type: 'file',
                                     )
                           .ifEmpty { exit 1, "ERROR: Could not find bwa indices from: ${bwa_index_path}." }
                           .collect()
                           .toList()

        data_channel = bwa_fasta.combine(bwa_index)

    emit:
        data_channel

}


workflow BWA_MAPPING_PE {

    take:
        data // a channel of [val(unique_sampleId), val(sampleId), path(fastq_PE1), path(fastq_PE2)]
        // unique_sampleId is used to label the read group field "SM" and (part of) "LB",
        // while sampleId represents each split fastq file for a unique sample.

    main:
        /*
           1) create a channel linking bwa index files from genome.fa in params, and
           2) combine this channel with the items in the data channel
        */
        bwa_inputs = get_bwa_index(params.tools.bwamaptools.bwa_fasta).combine(data)

        aligned_bam = BWA_MEM_PE(bwa_inputs)


        // publish output:
        PUBLISH_BAM(aligned_bam.map{it -> tuple(it[0], it[1])}, '.bwa.out.possorted.bam', 'bam')
        PUBLISH_BAM_INDEX(aligned_bam.map{it -> tuple(it[0], it[2])}, '.bwa.out.possorted.bam.bai', 'bam')
        MAPPING_SUMMARY(aligned_bam)
        PUBLISH_MAPPING_SUMMARY(MAPPING_SUMMARY.out, '.mapping_stats.tsv', 'reports/mapping_stats')

    emit:
        aligned_bam
}


workflow MARK_DUPLICATES {

    take:
        data // a channel of [val(sampleId), path(bam) ]
        mark_duplicates_method

    main:

        switch(mark_duplicates_method) {
            case 'MarkDuplicates':
                dup_marked_bam = PICARD__MARK_DUPLICATES_AND_SORT(data)
                PUBLISH_MARKDUPS_METRICS(dup_marked_bam.map{it -> tuple(it[0], it[3])}, '.mark_duplicates_metrics.txt', 'reports/mark_duplicates')
                break
            case 'MarkDuplicatesSpark':
                dup_marked_bam = GATK__MARK_DUPLICATES_SPARK(data)
                break
        }

        MAPPING_SUMMARY(dup_marked_bam.map { it -> tuple(it[0..2]) })
        PICARD__ESTIMATE_LIBRARY_COMPLEXITY(data)

        // publish output:
        PUBLISH_BAM(dup_marked_bam.map{it -> tuple(it[0], it[1])}, '.bwa.out.possorted.bam', 'bam')
        PUBLISH_BAM_INDEX(dup_marked_bam.map{it -> tuple(it[0], it[2])}, '.bwa.out.possorted.bam.bai', 'bam')
        PUBLISH_LIBRARY_METRICS(PICARD__ESTIMATE_LIBRARY_COMPLEXITY.out, '.library_complexity_metrics.txt', 'reports/mark_duplicates')
        PUBLISH_MAPPING_SUMMARY(MAPPING_SUMMARY.out, '.mapping_stats.tsv', 'reports/mapping_stats')

    emit:
        dup_marked_bam.map { it -> tuple(it[0..2]) }

}


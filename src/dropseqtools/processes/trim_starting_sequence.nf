nextflow.preview.dsl=2

process SC__DROP_SEQ_TOOLS__TRIM_SMART_UNALIGNED_TAGGED_FILTERED_BAM {

    container params.sc.dropseqtools.container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path('*.unaligned_tagged_trimmed_smart.bam'), emit: bam
        tuple file('*.adapter_trimming_report.txt'), emit: report

    script:
        processParams = params.sc.dropseqtools.trim_smart_unaligned_tagged_filtered_bam
        """
        TrimStartingSequence \
            INPUT=${bam} \
            OUTPUT=${sample}.unaligned_tagged_trimmed_smart.bam \
            OUTPUT_SUMMARY=${sample}.adapter_trimming_report.txt \
            SEQUENCE=${processParams.adapterSequence} \
            MISMATCHES=${processParams.mismatches} \
            NUM_BASES=${processParams.numBases}
        """

}
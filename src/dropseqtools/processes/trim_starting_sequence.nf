nextflow.enable.dsl=2

process SC__DROP_SEQ_TOOLS__TRIM_SMART_UNALIGNED_TAGGED_FILTERED_BAM {

    container params.getToolParams("dropseqtools").container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path('*.unaligned_tagged_trimmed_smart.bam'), emit: bam
        tuple file('*.adapter_trimming_report.txt'), emit: report

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.getToolParams("dropseqtools").trim_smart_unaligned_tagged_filtered_bam)
		processParams = sampleParams.local
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

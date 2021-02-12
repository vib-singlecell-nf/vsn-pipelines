nextflow.enable.dsl=2

process SC__DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART {

    container params.getToolParams("dropseqtools").container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path('*.unaligned_tagged_polyA_filtered.bam'), emit: bam
    tuple file('*.polyA_trimming_report.txt'), emit: report

    script:
    def sampleParams = params.parseConfig(sampleId, params.global, params.getToolParams("dropseqtools").trim_polya_unaligned_tagged_trimmed_smart)
		processParams = sampleParams.local
    """
    PolyATrimmer \
        INPUT=${bam} \
        OUTPUT=${sample}.unaligned_tagged_polyA_filtered.bam \
        OUTPUT_SUMMARY=${sample}.polyA_trimming_report.txt \
        MISMATCHES=${processParams.mismatches} \
        NUM_BASES=${processParams.numBases}
    """

}

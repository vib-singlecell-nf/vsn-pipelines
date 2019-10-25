nextflow.preview.dsl=2

process SC__DROP_SEQ_TOOLS__TRIM_SMART_UNALIGNED_TAGGED_FILTERED_BAM {
    
    container params.sc.dropseqtools.container
    publishDir "${params.outdir}/01.clean", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_trimmed_smart.bam'), emit: bam
        tuple file('*.adapter_trimming_report.txt'), emit: report
    script:
        """
		TrimStartingSequence \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_trimmed_smart.bam \
			OUTPUT_SUMMARY=${sample}.adapter_trimming_report.txt \
			SEQUENCE=${params.sc.dropseqtools.trim_smart_unaligned_tagged_filtered_bam.adapterSequence} \
			MISMATCHES=${params.sc.dropseqtools.trim_smart_unaligned_tagged_filtered_bam.mismatches} \
			NUM_BASES=${params.sc.dropseqtools.trim_smart_unaligned_tagged_filtered_bam.numBases}
        """
}
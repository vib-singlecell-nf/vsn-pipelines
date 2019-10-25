nextflow.preview.dsl=2

process SC__DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART {
    
    container params.sc.dropseqtools.container
    publishDir "${params.outdir}/01.clean", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_polyA_filtered.bam'), emit: bam
        tuple file('*.polyA_trimming_report.txt'), emit: report
    script:
        """
		PolyATrimmer \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_polyA_filtered.bam \
			OUTPUT_SUMMARY=${sample}.polyA_trimming_report.txt \
			MISMATCHES=${params.sc.dropseqtools.trim_polya_unaligned_tagged_trimmed_smart.mismatches} \
			NUM_BASES=${params.sc.dropseqtools.trim_polya_unaligned_tagged_trimmed_smart.numBases}
        """
}
nextflow.preview.dsl=2

process DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_polyA_filtered.bam'), emit: bam
        tuple file('*.polyA_trimming_report.txt'), emit: report
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		PolyATrimmer \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_polyA_filtered.bam \
			OUTPUT_SUMMARY=${sample}.polyA_trimming_report.txt \
			MISMATCHES=${params.mismatches} \
			NUM_BASES=${params.numBases}
        """
}
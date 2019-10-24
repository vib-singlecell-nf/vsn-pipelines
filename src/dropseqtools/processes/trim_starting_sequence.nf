nextflow.preview.dsl=2

process DROP_SEQ_TOOLS__TRIM_SMART_UNALIGNED_TAGGED_FILTERED_BAM {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_trimmed_smart.bam'), emit: bam
        tuple file('*.adapter_trimming_report.txt'), emit: report
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		TrimStartingSequence \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_trimmed_smart.bam \
			OUTPUT_SUMMARY=${sample}.adapter_trimming_report.txt \
			SEQUENCE=${params.adapterSequence} \
			MISMATCHES=${params.mismatches} \
			NUM_BASES=${params.numBases}
        """
}
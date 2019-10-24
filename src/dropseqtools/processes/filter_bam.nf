nextflow.preview.dsl=2

process DROP_SEQ_TOOLS__FILTER_UNALIGNED_TAGGED_BAM {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_filtered.bam'), emit: bam
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		FilterBAM \
			TAG_REJECT=${params.tagReject} \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_filtered.bam
        """
}
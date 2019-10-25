nextflow.preview.dsl=2

process SC__DROP_SEQ_TOOLS__FILTER_UNALIGNED_TAGGED_BAM {

    container params.sc.dropseqtools.container
    publishDir "${params.outdir}/01.clean", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_filtered.bam'), emit: bam
    script:
        """
		FilterBAM \
			TAG_REJECT=${params.sc.dropseqtools.filter_unaligned_tagged_bam.tagReject} \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_filtered.bam
        """
}
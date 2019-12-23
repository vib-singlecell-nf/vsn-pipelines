nextflow.preview.dsl=2

process SC__DROP_SEQ_TOOLS__FILTER_UNALIGNED_TAGGED_BAM {

    container params.sc.dropseqtools.container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path('*.unaligned_tagged_filtered.bam'), emit: bam

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sc.dropseqtools.filter_unaligned_tagged_bam)
		processParams = sampleParams.local
        """
        FilterBAM \
            TAG_REJECT=${processParams.tagReject} \
            INPUT=${bam} \
            OUTPUT=${sample}.unaligned_tagged_filtered.bam
        """

}
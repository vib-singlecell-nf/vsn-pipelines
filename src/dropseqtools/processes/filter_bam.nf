nextflow.enable.dsl=2

process SC__DROP_SEQ_TOOLS__FILTER_UNALIGNED_TAGGED_BAM {

    container params.tools.dropseqtools.container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path('*.unaligned_tagged_filtered.bam'), emit: bam

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.dropseqtools.filter_unaligned_tagged_bam)
		processParams = sampleParams.local
        """
        FilterBAM \
            TAG_REJECT=${processParams.tagReject} \
            INPUT=${bam} \
            OUTPUT=${sample}.unaligned_tagged_filtered.bam
        """

}

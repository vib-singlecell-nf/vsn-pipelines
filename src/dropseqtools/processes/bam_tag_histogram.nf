nextflow.enable.dsl=2

process SC__DROP_SEQ_TOOLS__BAM_TAG_HISTOGRAM {
    
    container params.getToolParams("dropseqtools").container
    publishDir "${params.global.outdir}/03.count", mode: 'symlink'
    label 'compute_resources__default'

    input:
        tuple val(sample), path(bam)
    
    output:
	    tuple val(sample), path("*.cell_readcounts.txt.gz")
    
    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.getToolParams("dropseqtools").bam_tag_histogram)
		processParams = sampleParams.local
        """
        BAMTagHistogram \
            I=${bam} \
            O=${sample}.cell_readcounts.txt.gz \
            TAG=${processParams.tag}
        """

}

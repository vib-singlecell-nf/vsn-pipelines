nextflow.preview.dsl=2

process SC__DROP_SEQ_TOOLS__BAM_TAG_HISTOGRAM {
    
    container params.sc.dropseqtools.container
    publishDir "${params.global.outdir}/03.count", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple val(sample), path(bam)
    
    output:
	    tuple val(sample), path("*.cell_readcounts.txt.gz")
    
    script:
        processParams = params.sc.dropseqtools.bam_tag_histogram
        """
        BAMTagHistogram \
            I=${bam} \
            O=${sample}.cell_readcounts.txt.gz \
            TAG=${processParams.tag}
        """

}
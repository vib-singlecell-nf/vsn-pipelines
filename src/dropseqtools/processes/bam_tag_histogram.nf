nextflow.preview.dsl=2

process SC__DROP_SEQ_TOOLS__BAM_TAG_HISTOGRAM {
    
    container params.sc.dropseqtools.container
    publishDir "${params.outdir}/03.count", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
		tuple val(sample), file("*.cell_readcounts.txt.gz")
    script:
        """
		BAMTagHistogram \
			I=${bam} \
			O=${sample}.cell_readcounts.txt.gz \
			TAG=${params.sc.dropseqtools.bam_tag_histogram.tag}
        """
}
nextflow.preview.dsl=2

process DROP_SEQ_TOOLS__BAM_TAG_HISTOGRAM {
    publishDir "${params.outdir}/03.count", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
		tuple val(sample), file("*.cell_readcounts.txt.gz")
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		BAMTagHistogram \
			I=${bam} \
			O=${sample}.cell_readcounts.txt.gz \
			TAG=${params.dropseqtools.bam_tag_histogram.tag}
        """
}
nextflow.enable.dsl=2

process PICARD__FASTQ_TO_BAM {

    container params.getToolParams("picard").container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        file(reads)
        file(tmpDir)
    
    output:
        tuple val(sample), path('*.unaligned.bam'), emit: bam
    
    script:
        sample = reads[0].toString() - ~/(_R1)?(\.clean)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        java -Djava.io.tmpdir=$tmpDir -jar \
            /picard.jar \
                FastqToSam \
                    FASTQ=${reads[0]} \
                    FASTQ2=${reads[1]} \
                    O=${sample}.unaligned.bam \
                    SAMPLE_NAME=${sample}
        """

}

nextflow.preview.dsl=2

process PICARD__FASTQ_TO_BAM {

    container params.picard.container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    file(reads)
    file(tmpDir)
    
    output:
    tuple val(sample), file('*.unaligned.bam'), emit: bam
    
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
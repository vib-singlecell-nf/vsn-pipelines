nextflow.preview.dsl=2

process PICARD__BAM_TO_FASTQ {

    container params.picard.container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    tuple val(sample), file(bam)
    file(tmpDir)

    output:
    tuple val(sample), file('*.unaligned_tagged_polyA_filtered.fastq'), emit: fastq

    script:
    """
    java -Djava.io.tmpdir=$tmpDir -jar \
        /picard.jar \
            SamToFastq \
                INPUT=${bam} \
                FASTQ=${sample}.unaligned_tagged_polyA_filtered.fastq
    """

}
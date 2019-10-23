nextflow.preview.dsl=2

process PICARD__FASTQ_TO_BAM {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        file(reads)
    output:
        tuple val(sample), file('*.unaligned.bam'), emit: bam
    script:
        sample = reads[0].toString() - ~/(_R1)?(\.clean)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        module load Java/1.8.0_192
        java -Djava.io.tmpdir=$DWMAX/tmp -jar \
            /ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/software/genius/picard/2.20.6/bin/picard.jar \
            FastqToSam \
                    FASTQ=${reads[0]} \
                    FASTQ2=${reads[1]} \
                    O=${sample}.unaligned.bam \
                    SAMPLE_NAME=${sample}
        """
}
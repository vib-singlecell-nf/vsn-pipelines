nextflow.enable.dsl=2

process PICARD__BAM_TO_FASTQ {

    container params.tools.picard.container
    publishDir "${params.global.outdir}/01.clean", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        tuple val(sample), path(bam)
        file(tmpDir)

    output:
        tuple val(sample), path('*.unaligned_tagged_polyA_filtered.fastq'), emit: fastq

    script:
        """
        java -Djava.io.tmpdir=$tmpDir -jar \
            /picard.jar \
                SamToFastq \
                    INPUT=${bam} \
                    FASTQ=${sample}.unaligned_tagged_polyA_filtered.fastq
        """

}

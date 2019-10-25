nextflow.preview.dsl=2

process PICARD__SORT_SAM {

    container params.picard.container
    publishDir "${params.outdir}/02.map", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
        file(tmpDir)
    output:
        tuple val(sample), file("*.STAR_aligned_sorted.bam")
    script:
        """
        java -Djava.io.tmpdir=$tmpDir -jar \
            /picard.jar \
                SortSam \
                    I=${bam} \
                    O=${sample}.STAR_aligned_sorted.bam \
                    SO=${params.picard.sort_sam.so}
        """
}
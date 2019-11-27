nextflow.preview.dsl=2

process PICARD__SORT_SAM {

    container params.picard.container
    publishDir "${params.global.outdir}/02.map", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    tuple val(sample), path(bam)
    file(tmpDir)

    output:
    tuple val(sample), path("*.STAR_aligned_sorted.bam")
    script:
    processParams = params.picard.sort_sam
    """
    java -Djava.io.tmpdir=$tmpDir -jar \
        /picard.jar \
            SortSam \
                I=${bam} \
                O=${sample}.STAR_aligned_sorted.bam \
                SO=${processParams.so}
    """

}
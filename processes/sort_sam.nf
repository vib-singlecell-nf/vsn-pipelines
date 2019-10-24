nextflow.preview.dsl=2

process PICARD__SORT_SAM {
    publishDir "${params.outdir}/02.map", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file("*.STAR_aligned_sorted.bam")
    script:
        """
        module load Java/1.8.0_192
        java -Djava.io.tmpdir=$DWMAX/tmp -jar \
            /ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/software/genius/picard/2.20.6/bin/picard.jar \
                SortSam \
                    I=${bam} \
                    O=${sample}.STAR_aligned_sorted.bam \
                    SO=${params.picard.sort_sam.so}
        """
}
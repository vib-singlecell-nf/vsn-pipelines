nextflow.preview.dsl=2

process PICARD__MERGE_BAM_ALIGNMENT {

    container params.picard.container
    publishDir "${params.global.outdir}/02.map", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=24:00:00 -A ${params.global.qsubaccount}"

    input:
    tuple val(sample), path(unmappedBam)
    tuple val(sample), path(mappedBam)
    file(genome)
    file(dict)
    file(tmpDir)

    output:
    tuple val(sample), path("*.merged.bam")

    script:
    processParams = params.picard.merge_bam_alignment
    """
    java -Djava.io.tmpdir=$tmpDir -jar \
        /picard.jar \
            MergeBamAlignment \
                REFERENCE_SEQUENCE=${genome} \
                UNMAPPED_BAM=${unmappedBam} \
                ALIGNED_BAM=${mappedBam} \
                OUTPUT=${sample}.merged.bam \
                INCLUDE_SECONDARY_ALIGNMENTS=${processParams.includeSecondaryAlignments} \
                PAIRED_RUN=${processParams.pairedRun}
    """

}
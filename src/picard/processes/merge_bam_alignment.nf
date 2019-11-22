nextflow.preview.dsl=2

process PICARD__MERGE_BAM_ALIGNMENT {

    container params.picard.container
    publishDir "${params.outdir}/02.map", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(unmappedBam)
        tuple val(sample), file(mappedBam)
        file(genome)
        file(dict)
        file(tmpDir)
    output:
        tuple val(sample), file("*.merged.bam")
    script:
        """
        java -Djava.io.tmpdir=$tmpDir -jar \
            /picard.jar \
                MergeBamAlignment \
                    REFERENCE_SEQUENCE=${genome} \
                    UNMAPPED_BAM=${unmappedBam} \
                    ALIGNED_BAM=${mappedBam} \
                    OUTPUT=${sample}.merged.bam \
                    INCLUDE_SECONDARY_ALIGNMENTS=${params.picard.merge_bam_alignment.includeSecondaryAlignments} \
                    PAIRED_RUN=${params.picard.merge_bam_alignment.pairedRun}
        """    
}
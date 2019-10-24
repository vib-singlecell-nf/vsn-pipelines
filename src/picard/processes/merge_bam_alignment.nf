nextflow.preview.dsl=2

process PICARD__MERGE_BAM_ALIGNMENT {
    publishDir "${params.outdir}/02.map", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(unmappedBam)
        tuple val(sample), file(mappedBam)
        file(genome)
        file(dict)
    output:
        tuple val(sample), file("*.merged.bam")
    script:
        """
        module load Java/1.8.0_192
        java -Djava.io.tmpdir=$DWMAX/tmp -jar \
            /ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/software/genius/picard/2.20.6/bin/picard.jar \
                MergeBamAlignment \
                    REFERENCE_SEQUENCE=${genome} \
                    UNMAPPED_BAM=${unmappedBam} \
                    ALIGNED_BAM=${mappedBam} \
                    OUTPUT=${sample}.merged.bam \
                    INCLUDE_SECONDARY_ALIGNMENTS=${params.picard.merge_bam_alignment.includeSecondaryAlignments} \
                    PAIRED_RUN=${params.picard.merge_bam_alignment.pairedRun}
        """    
}
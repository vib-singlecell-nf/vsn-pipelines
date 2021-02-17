nextflow.enable.dsl=2

process PICARD__MERGE_BAM_ALIGNMENT {

    container params.tools.picard.container
    publishDir "${params.global.outdir}/02.map", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        tuple val(sample), path(unmappedBam)
        tuple val(sample), path(mappedBam)
        file(genome)
        file(dict)
        file(tmpDir)

    output:
        tuple val(sample), path("*.merged.bam")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.picard.merge_bam_alignment)
		processParams = sampleParams.local
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

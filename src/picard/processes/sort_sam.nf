nextflow.enable.dsl=2

process PICARD__SORT_SAM {

    container params.tools.picard.container
    publishDir "${params.global.outdir}/02.map", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        tuple val(sample), path(bam)
        file(tmpDir)

    output:
        tuple val(sample), path("*.STAR_aligned_sorted.bam")
    
    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.picard.sort_sam)
		processParams = sampleParams.local
        """
        java -Djava.io.tmpdir=$tmpDir -jar \
            /picard.jar \
                SortSam \
                    I=${bam} \
                    O=${sample}.STAR_aligned_sorted.bam \
                    SO=${processParams.so}
        """

}

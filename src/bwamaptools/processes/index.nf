nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.bwamaptools

process SC__BWAMAPTOOLS__INDEX_BAM {

    container toolParams.container
    label 'compute_resources__default'

    input:
        tuple val(sampleId),
              path(bam)

    output:
        tuple val(sampleId),
              path("*.bai")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        processParams = sampleParams.local
        """
        samtools index ${bam}
        """
}

process SC__BWAMAPTOOLS__INDEX_BED {

    container toolParams.container
    label 'compute_resources__default'

    input:
        tuple val(sampleId),
              path(bed)

    output:
        tuple val(sampleId),
              path("*.tbi")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        processParams = sampleParams.local
        """
        tabix -p bed ${bed}
        """
}


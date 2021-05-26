nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.gatk

process PICARD__ESTIMATE_LIBRARY_COMPLEXITY {

    container toolParams.container
    label 'compute_resources__default','compute_resources__24hqueue'

    input:
        tuple val(sampleId),
              path(bam)

    output:
        tuple val(sampleId),
              path("${sampleId}.picard_library_complexity_metrics.txt")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.estimate_library_complexity)
        processParams = sampleParams.local
        """
        gatk EstimateLibraryComplexity \
            -I ${bam} \
            -O ${sampleId}.picard_library_complexity_metrics.txt \
            --BARCODE_TAG ${processParams.barcode_tag} \
        """
}


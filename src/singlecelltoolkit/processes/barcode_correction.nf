nextflow.enable.dsl=2

//binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/singlecelltoolkit/bin/" : ""

toolParams = params.tools.singlecelltoolkit

process SC__SINGLECELLTOOLKIT__BARCODE_CORRECTION {

    container toolParams.container
    label 'compute_resources__default'

    input:
        tuple val(sampleId),
              path(fastq_bc),
              path(bc_whitelist)

    output:
        tuple val(sampleId),
              path("${sampleId}_bc_corrected.fastq.gz"),
              path("${sampleId}_bc_corrected.fastq.gz.corrected.bc_stats.tsv")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        processParams = sampleParams.local
        """
        correct_barcode_in_fastq.sh \
            ${bc_whitelist} \
            ${fastq_bc} \
            ${sampleId}_bc_corrected.fastq.gz
        """
}


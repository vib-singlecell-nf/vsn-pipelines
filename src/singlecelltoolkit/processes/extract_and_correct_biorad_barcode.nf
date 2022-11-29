nextflow.enable.dsl=2

//binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/singlecelltoolkit/bin/" : ""

toolParams = params.tools.singlecelltoolkit

process SCTK__EXTRACT_AND_CORRECT_BIORAD_BARCODE {

    container toolParams.container
    label 'compute_resources__sctk_barcode'

    input:
        tuple val(sampleId),
              val(technology),
              path(fastq_PE1),
              path(fastq_PE2)

    output:
        tuple val(sampleId),
              path("${sampleId}_dex_R1.fastq.gz"),
              path("${sampleId}_dex_R2.fastq.gz"),
              path("${sampleId}_dex.corrected_bc_stats.tsv")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        //processParams = sampleParams.local
        """
        extract_and_correct_biorad_barcode_in_fastq.sh \
            ${fastq_PE1} \
            ${fastq_PE2} \
            ${sampleId}_dex
        """
}


nextflow.enable.dsl=2

//binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/singlecelltoolkit/bin/" : ""

toolParams = params.tools.singlecelltoolkit

process SCTK__BARCODE_CORRECTION {

    container toolParams.container
    label 'compute_resources__sctk_barcode'

    input:
        tuple val(sampleId),
              val(technology),
              path(fastq_PE1),
              path(fastq_bc),
              path(fastq_PE2),
              path(bc_whitelist)

    output:
        tuple val(sampleId),
              val(technology),
              path(fastq_PE1),
              path("${sampleId}_bc_corrected.fastq.gz"),
              path(fastq_PE2),
              path("${sampleId}_bc_corrected.fastq.gz.corrected.bc_stats.tsv")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.barcode_correction)
        processParams = sampleParams.local
        """
        correct_barcode_in_fastq.sh \
            ${bc_whitelist} \
            ${fastq_bc} \
            ${sampleId}_bc_corrected.fastq.gz \
            ${processParams.max_mismatches} \
            ${processParams.min_frac_bcs_to_find}
        """
}


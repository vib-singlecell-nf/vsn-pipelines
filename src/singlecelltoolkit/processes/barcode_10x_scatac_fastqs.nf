nextflow.enable.dsl=2

//binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/singlecelltoolkit/bin/" : ""

toolParams = params.tools.singlecelltoolkit

process SCTK__BARCODE_10X_SCATAC_FASTQ {

    container toolParams.container
    label 'compute_resources__sctk__barcode_10x_scatac_fastq_5cpus'

    input:
        tuple val(sampleId),
              val(technology),
              path(fastq_PE1),
              path(fastq_bc),
              path(fastq_PE2)

    output:
        tuple val(sampleId),
              path("${sampleId}_dex_R1.fastq.gz"),
              path("${sampleId}_dex_R2.fastq.gz")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.barcode_10x_scatac_fastqs)
        processParams = sampleParams.local
        def max_threads = (task.cpus > 5) ? 5 : task.cpus
        """
        export compress_fastq_threads="${max_threads}"
        barcode_10x_scatac_fastqs.sh \
            ${fastq_PE1} \
            ${fastq_bc} \
            ${fastq_PE2} \
            ${sampleId}_dex \
            false \
            true \
            ${processParams.uncorrected_bc_tag}_${processParams.barcode_quality_tag}
        """
}


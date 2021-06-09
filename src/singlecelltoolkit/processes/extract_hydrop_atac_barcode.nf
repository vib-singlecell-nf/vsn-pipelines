nextflow.enable.dsl=2

//binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/singlecelltoolkit/bin/" : ""

toolParams = params.tools.singlecelltoolkit

process SCTK__EXTRACT_HYDROP_ATAC_BARCODE {

    container toolParams.container
    label 'compute_resources__default'

    input:
        tuple val(sampleId),
              val(technology),
              path(fastq_R1),
              path(fastq_R2),
              path(fastq_R3)

    output:
        tuple val(sampleId),
              val(technology),
              path(fastq_R1),
              path("${sampleId}_hydrop_barcode_R2.fastq.gz"),
              path(fastq_R3)

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        //processParams = sampleParams.local
        """
        extract_hydrop_atac_barcode_from_R2_fastq.sh \
            ${fastq_R2} \
            ${sampleId}_hydrop_barcode_R2.fastq.gz \
            igzip
        """
}


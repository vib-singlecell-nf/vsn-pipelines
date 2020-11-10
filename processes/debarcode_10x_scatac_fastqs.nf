nextflow.preview.dsl=2

//binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/singlecelltoolkit/bin/" : ""

toolParams = params.sc.atac.singlecelltoolkit

process SC__SINGLECELLTOOLKIT__DEBARCODE_10X_FASTQ {

    container toolParams.container
    publishDir "${params.global.outdir}/fastq/barcode_demultiplexed", mode: 'symlink'
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId),
              path(fastq_PE1),
              path(fastq_bc),
              path(fastq_PE2)

    output:
        tuple val(sampleId),
              path("${sampleId}_dex_R1.fastq.gz"),
              path("${sampleId}_dex_R2.fastq.gz")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        processParams = sampleParams.local
        """
        debarcode_10x_scatac_fastqs.sh \
            ${fastq_PE1} \
            ${fastq_bc} \
            ${fastq_PE2} \
            ${sampleId}_dex
        """
}


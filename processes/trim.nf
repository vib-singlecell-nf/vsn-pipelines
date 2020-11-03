nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

process SC__TRIMGALORE__TRIM {

    container params.sc.atac.trimgalore.container
    publishDir "${params.global.outdir}/fastq/trimgalore", mode: 'symlink'
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        tuple val(sampleId),
              path(fastq_PE1),
              path(fastq_PE2)

    output:
        tuple val(sampleId),
              path("${sampleId}.R1.dex_val_1.fq.gz"),
              path("${sampleId}.R2.dex_val_2.fq.gz"),
              path("*trimming_report.txt")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sc.atac.trimgalore.trim)
        processParams = sampleParams.local
        """
        trim_galore \
            -j ${task.cpus} \
            -o . \
            ${fastq_PE1} \
            ${fastq_PE2} \
            --paired \
            --gzip
        """
}


nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.fastp

process FASTP__ADAPTER_TRIMMING {

    container toolParams.container
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        tuple val(sampleId),
              path(fastq_PE1),
              path(fastq_PE2)

    output:
        tuple val(sampleId),
              path("${sampleId}_dex_R1_val_1.fq.gz"),
              path("${sampleId}_dex_R2_val_2.fq.gz"),
              path("${sampleId}_fastp.html")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        processParams = sampleParams.local
        def max_threads = (task.cpus > 6) ? 6 : task.cpus
        """
        fastp \
            --in1 ${fastq_PE1} \
            --in2 ${fastq_PE2} \
            --out1 ${sampleId}_dex_R1_val_1.fq.gz \
            --out2 ${sampleId}_dex_R2_val_2.fq.gz \
            --detect_adapter_for_pe \
            --html ${sampleId}_fastp.html \
            --thread ${max_threads}
        """
}


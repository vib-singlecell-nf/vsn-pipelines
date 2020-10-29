nextflow.preview.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.sc.atac.bwamaptools

process SC__BWAMAPTOOLS__BWA_MEM_PE {

    container toolParams.container
    label 'compute_resources__bwa_mem'

    input:
        tuple path(bwa_index),
              val(sampleId),
              path(fastq_PE1),
              path(fastq_PE2)

    output:
        tuple val(sampleId),
              path("${sampleId}.bwa.out.possorted.bam")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        processParams = sampleParams.local
        """
        bwa mem \
            -t ${task.cpus} \
            ${toolParams.index} \
            ${fastq_PE1} \
            ${fastq_PE2} \
            | samtools view -bS - \
            | samtools sort -@ ${task.cpus} - -o ${sampleId}.bwa.out.possorted.bam
        """
}


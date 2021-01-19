nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.bwamaptools

process SC__BWAMAPTOOLS__BWA_MEM_PE {

    container toolParams.container
    label 'compute_resources__bwa_mem'

    input:
        tuple path(bwa_fasta),
              path(bwa_index),
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
        set -euo pipefail
        bwa mem \
            -t ${task.cpus} \
            ${bwa_fasta} \
            ${fastq_PE1} \
            ${fastq_PE2} \
        | samtools sort -@ ${task.cpus} -n -O bam - \
        | samtools fixmate -@ ${task.cpus} -m -O bam - - \
        | samtools sort -@ ${task.cpus} -O bam - \
        | samtools markdup -@ ${task.cpus} -f ${sampleId}.markdup.log - ${sampleId}.bwa.out.possorted.bam
        """

}


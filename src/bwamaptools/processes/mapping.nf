nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.bwamaptools

process BWAMAPTOOLS__BWA_MEM_PE {

    container toolParams.container
    label 'compute_resources__bwa_mem'

    input:
        tuple path(bwa_fasta),
              path(bwa_index),
              val(unique_sampleId),
              val(sampleId),
              path(fastq_PE1),
              path(fastq_PE2)

    output:
        tuple val(sampleId),
              path("${sampleId}.bwa.out.fixmate.possorted.bam"),
              path("${sampleId}.bwa.out.fixmate.possorted.bam.bai")

    script:
        def sampleParams = params.parseConfig(unique_sampleId, params.global, toolParams)
        processParams = sampleParams.local
        """
        id=\$(zcat ${fastq_PE1} | head -n 1 | cut -f 1-4 -d':' | sed 's/@//')
        ${toolParams.bwa_version} mem \
            -t ${task.cpus} \
            -C \
            -R "@RG\\tID:\${id}\\tSM:${unique_sampleId}\\tLB:\${id}"__"${unique_sampleId}\\tPL:ILLUMINA" \
            ${bwa_fasta} \
            ${fastq_PE1} \
            ${fastq_PE2} \
        | samtools fixmate -u -m -O bam - - \
        | samtools sort -@ 2 -m 2G -O bam --write-index -T '${sampleId}.bwa.out.fixmate.possorted.TMP' -o '${sampleId}.bwa.out.fixmate.possorted.bam##idx##${sampleId}.bwa.out.fixmate.possorted.bam.bai' -
        """
}

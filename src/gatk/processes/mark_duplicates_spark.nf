nextflow.enable.dsl=2

//binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.gatk

process GATK__MARK_DUPLICATES_SPARK {

    container toolParams.container
    label 'compute_resources__cpu','compute_resources__24hqueue'

    input:
        tuple val(sampleId),
              path(bam)

    output:
        tuple val(sampleId),
              path("${sampleId}.bwa.out.fixmate.picard_markdup.possorted.bam"),
              path("${sampleId}.bwa.out.fixmate.picard_markdup.possorted.bam.bai")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        //processParams = sampleParams.local
        """
        gatk MarkDuplicatesSpark \
            -I ${bam} \
            -O ${sampleId}.bwa.out.fixmate.picard_markdup.possorted.bam \
            -- \
            --spark-runner LOCAL \
            --spark-master local[${task.cpus}]
        """
}


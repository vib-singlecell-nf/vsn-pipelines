nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.picard

process PICARD__MARK_DUPLICATES_AND_SORT {

    container toolParams.container
    label 'compute_resources__default','compute_resources__24hqueue'

    input:
        tuple val(sampleId),
              path(bam)

    output:
        tuple val(sampleId),
              path("${sampleId}.bwa.out.fixmate.picard_markdup.possorted.bam"),
              path("${sampleId}.bwa.out.fixmate.picard_markdup.possorted.bai")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        processParams = sampleParams.local
        """
        set -euo pipefail
        java -jar /picard.jar MarkDuplicates \
            I=${bam} \
            O=/dev/stdout \
            BARCODE_TAG=CB \
            COMPRESSION_LEVEL=0 \
            QUIET=true \
            ASSUME_SORT_ORDER=queryname \
        | java -jar /picard.jar SortSam \
            I=/dev/stdin \
            O=${sampleId}.bwa.out.fixmate.picard_markdup.possorted.bam \
            SORT_ORDER=coordinate \
            CREATE_INDEX=true
        """
}


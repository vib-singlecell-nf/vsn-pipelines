nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.gatk

process PICARD__MARK_DUPLICATES_AND_SORT {

    container toolParams.container
    label 'compute_resources__picard__mark_duplicates_and_sort'

    input:
        tuple val(sampleId),
              path(bam)

    output:
        tuple val(sampleId),
              path("${sampleId}.bwa.out.fixmate.picard_markdup.possorted.bam"),
              path("${sampleId}.bwa.out.fixmate.picard_markdup.possorted.bai"),
              path("${sampleId}.bwa.out.fixmate.picard_markdup.metrics.txt")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        processParams = sampleParams.local
        """
        set -euo pipefail
        gatk MarkDuplicates \
            -I ${bam} \
            -O /dev/stdout \
            --METRICS_FILE ${sampleId}.bwa.out.fixmate.picard_markdup.metrics.txt \
            --BARCODE_TAG CB \
            --COMPRESSION_LEVEL 0 \
            --QUIET true \
            --ASSUME_SORT_ORDER queryname \
        | gatk SortSam \
            -I /dev/stdin \
            -O ${sampleId}.bwa.out.fixmate.picard_markdup.possorted.bam \
            --SORT_ORDER coordinate \
            --CREATE_INDEX true
        """
}


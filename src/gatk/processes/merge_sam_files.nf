nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.gatk

process PICARD__MERGE_SAM_FILES_AND_SORT {

    container toolParams.container
    label 'compute_resources__default','compute_resources__24hqueue'

    input:
        tuple val(sampleId),
              path(bams)

    output:
        tuple val(sampleId),
              path("${sampleId}.bwa.out.fixmate.merged.bam")

    script:
        //def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        //processParams = sampleParams.local
        """
        gatk MergeSamFiles \
            ${"-I "+bams.join(" -I ")} \
            -O /dev/stdout \
        | gatk SortSam \
            -I /dev/stdin \
            -O ${sampleId}.bwa.out.fixmate.merged.bam \
            --SORT_ORDER queryname
        """
}


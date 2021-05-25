nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.picard

process PICARD__MERGE_SAM_FILES {

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
        java -jar /picard.jar MergeSamFiles \
            ${"I="+bams.join(" I=")} \
            O=/dev/stdout \
        | java -jar /picard.jar SortSam \
            I=/dev/stdin \
            O=${sampleId}.bwa.out.fixmate.merged.bam \
            SORT_ORDER=queryname
        """
}


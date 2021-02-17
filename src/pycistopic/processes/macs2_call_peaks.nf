nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.getToolParams("pycistopic")

process SC__PYCISTOPIC__MACS2_CALL_PEAKS {

    container toolParams.container
    publishDir "${params.global.outdir}/peaks/", mode: 'symlink'
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId),
              path(bam),
              path(bam_index),
              val(filetype)

    output:
        tuple val(sampleId),
              path("${sampleId}_peaks.narrowPeak"),
              path("${sampleId}_summits.bed")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        processParams = sampleParams.local
        """
        macs2 callpeak \
            --treatment ${bam} \
            --name ${sampleId} \
            --outdir . \
            --format BAMPE \
            --gsize ${processParams.gsize} \
            --qvalue ${processParams.qvalue} \
            --nomodel \
            --shift ${processParams.shift} \
            --extsize ${processParams.extsize} \
            --keep-dup ${processParams.keepdup} \
            --call-summits \
            --nolambda
        """
}


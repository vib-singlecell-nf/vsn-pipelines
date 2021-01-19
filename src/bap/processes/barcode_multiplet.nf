nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.bap

process SC__BAP__BARCODE_MULTIPLET_PIPELINE {

    container toolParams.container
    publishDir "${params.global.outdir}/bap", mode: 'symlink'
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId),
              path(bam),
              path(bai)

    output:
        tuple val(sampleId),
              path("${sampleId}.SC__TEMPLATE__PROCESS1.h5ad")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.barcode_multiplet)
		processParams = sampleParams.local
        """
        bap2 bam \
            -i ${bam} \
            -o bap_output \
            -r ${toolParams.genome} \
            -c ${task.cpus} \
            -bt ${toolParams.barcode}
        """
}


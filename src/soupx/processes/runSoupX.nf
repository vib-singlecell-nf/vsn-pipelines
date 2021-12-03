nextflow.enable.dsl=2

import java.nio.file.Paths

moduleName = "soupx"
toolParams = params.tools[moduleName]
binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/${moduleName}/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "${moduleName}/bin")


process SC__SOUPX {
    
    container toolParams.container
    publishDir "${params.global.outdir}/data/${moduleName}", mode: 'link'
    label 'compute_resources__default'

    input:
        tuple \
            val(sampleId), \
            path(f)

    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__SOUPX.Rds"), \
            emit: main
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__SOUPX.*.pdf"), \
            emit: other

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
        processParams = sampleParams.local

        def roundToIntAsArgument = ''
        if(processParams?.roundToInt) {
            roundToIntAsArgument = '--round-to-int '+ processParams.roundToInt
        }

        """
        export NXF_BIN_DIR=$binDir
        ${binDir}/run_soupx.R \
            --sample-id ${sampleId} \
            --seed ${params.global.seed} \
            ${roundToIntAsArgument} \
            --output-prefix "${sampleId}.SC__SOUPX" \
            $f
        """

}

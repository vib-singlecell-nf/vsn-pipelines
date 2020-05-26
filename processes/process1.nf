nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

process SC__TEMPLATE__PROCESS1 {

    container params.sc.template.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__TEMPLATE__PROCESS1.h5ad")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sc.template)
		processParams = sampleParams.local
        //varsUseAsArguments = processParams.varsUse.collect({ '--vars-use' + ' ' + it }).join(' ')
        """
        ${binDir}process1.py \
            --input ${f} \
            --output ${sampleId}.SC__TEMPLATE__PROCESS1.h5ad
        """
}


nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/harmony/bin/" : ""

process SC__HARMONY__HARMONY_MATRIX {
    
    container params.tools.harmony.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'
    label 'compute_resources__default'

    input:
        tuple \
            val(sampleId), \
            path(f)

    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__HARMONY__HARMONY_MATRIX.tsv")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.harmony)
		processParams = sampleParams.local
        varsUseAsArguments = processParams.varsUse.collect({ '--vars-use' + ' ' + it }).join(' ')
        """
        ${binDir}run_harmony.R \
            ${f} \
            --seed ${params.global.seed} \
            ${varsUseAsArguments} \
            ${processParams?.theta ? "--theta "+ processParams.theta : "" } \
            ${processParams?.lambda ? "--lambda "+ processParams.lambda : "" } \
            ${processParams?.epsilonHarmony ? "--epsilon-harmony "+ processParams.epsilonHarmony : "" } \
            --output-prefix "${sampleId}.SC__HARMONY__HARMONY_MATRIX"
        """

}

nextflow.preview.dsl=2

import java.nio.file.Paths

moduleName = "celda"
binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/${moduleName}/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "${moduleName}/bin")


process SC__CELDA__DECONTX {
    
    container params.sc.celda.container
    publishDir "${params.global.outdir}/data/${moduleName}", mode: 'link'
    label 'compute_resources__default'

    input:
        tuple \
            val(sampleId), \
            path(f)

    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.CELDA__DECONTX.Rds"), \
            emit: main
        tuple \
            val(sampleId), \
            path("${sampleId}.CELDA__DECONTX.{*.pdf,*.tsv}"), \
            emit: other

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.sc.celda.decontx)
        processParams = sampleParams.local
        
        def filterNumMadsThresholdsAsArguments = processParams.filters.containsKey('numMadsThresholds') ?
                processParams.filters.numMadsThresholds.collect({
                    '--num-mads-threshold ' + ' ' + it 
                }).join(' ') :
                ''
        def filterContaminationScoreThresholdsAsArguments = processParams.filters.containsKey('contaminationScoreThresholds') ?
                processParams.filters.contaminationScoreThresholds.collect({ 
                    '--custom-threshold ' + ' ' + it
                }).join(' ') :
                ''

        """
        ${binDir}/run_decontx.R \
            --sample-id ${sampleId} \
            --seed ${params.global.seed} \
            ${filterNumMadsThresholdsAsArguments} \
            ${filterContaminationScoreThresholdsAsArguments} \
            --output-prefix "${sampleId}.CELDA__DECONTX" \
            $f
        """

}

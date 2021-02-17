nextflow.enable.dsl=2

import java.nio.file.Paths

moduleName = "celda"
binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/${moduleName}/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "${moduleName}/bin")


process SC__CELDA__DECONTX {
    
    container params.tools.celda.container
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
            path("${sampleId}.CELDA__DECONTX.Contamination_Outlier_Table.tsv"), \
            emit: outlier_table
        tuple \
            val(sampleId), \
            path("${sampleId}.CELDA__DECONTX.{*.pdf,*.tsv}"), \
            emit: other

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.celda.decontx)
        processParams = sampleParams.local
        
        def filterNumMadsThresholdsAsArguments = ''
        def filterContaminationScoreThresholdsAsArguments = ''

        if(processParams?.filters) {
            filterNumMadsThresholdsAsArguments = processParams.filters.containsKey('numMadsThresholds') ?
                processParams.filters.numMadsThresholds.collect({
                    '--num-mads-threshold ' + ' ' + it 
                }).join(' ') :
                ''
            filterContaminationScoreThresholdsAsArguments = processParams.filters.containsKey('contaminationScoreThresholds') ?
                processParams.filters.contaminationScoreThresholds.collect({ 
                    '--custom-threshold ' + ' ' + it
                }).join(' ') :
                ''
        }

        def roundToIntAsArgument = ''
        if(processParams?.roundToInt) {
            roundToIntAsArgument = '--round-to-int '+ processParams.roundToInt
        }
        def filterEmptyCells = ''
        if(processParams?.filterEmptyCells) {
            filterEmptyCells = '--filter-empty-cells '+ processParams.filterEmptyCells
        }

        """
        ${binDir}/run_decontx.R \
            --sample-id ${sampleId} \
            --seed ${params.global.seed} \
            ${filterNumMadsThresholdsAsArguments} \
            ${filterContaminationScoreThresholdsAsArguments} \
            ${roundToIntAsArgument} \
            ${filterEmptyCells} \
            --output-prefix "${sampleId}.CELDA__DECONTX" \
            $f
        """

}
